// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Description:
//      fits2itk is a program that reads a FITS file into ITK and then uses ITK
//      to write the data back out in a form that can be understood by 3D
//      Slicer or another multidimensional imaging program.
//
// Author:
//      Douglas Alan <douglas_alan AT harvard.edu>
//                   <doug AT alum.mit.edu>
//      Initiative in Innovative Computing at Harvard University
//
// Copyright (c) 2006-2007 Harvard University
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details.
//=============================================================================

#include <libgen.h>			   // For basename()
#include <getopt.h>
#include <sys/param.h>                     // For MAXPATHLEN

#include <cassert>
#include <fstream>
using std::ifstream;
using std::istream;
using std::ostream;
using std::ios;
#include <memory>
using std::auto_ptr;
#include <string>
using std::string;

#include <itkFITSImageIOFactory.h>
#include <itkFITSImageIO.h>
#include <itkImage.h>
using itk::Image;

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkFlipImageFilter.h>
#include <itkBinomialBlurImageFilter.h>
// #include <itkDerivativeImageFilter.h>
// #include <itkMeanImageFilter.h>
// #include <itkBinaryMedianImageFilter.h>
// #include <itkGradientAnisotropicDiffusionImageFilter.h>

#include <pathToExecutable.h>
#include <da_util.h>
#include <da_sugar.h>

extern const char fits2itkVersion[];
const int c_dims = itk::FITSImageIO::c_dims;


//-----------------------------------------------------------------------------
// usage()
//-----------------------------------------------------------------------------

local proc void
usage1(bool exitWithFailureP, bool verboseUsageP)
{
//  if (daProgramName().size()) cerr << basename(daProgramName().c_str());
//  else cerr << "fits2itk";
  
  ostream& out = exitWithFailureP ? cerr : cout;
  out << "fits2itk " << fits2itkVersion << "\n\n";
  out << "usage:\n";

  const char* shortUsageFilepath = pteJoinPath(pathToExecutableDir(),
					       "shortUsageMessage.txt");

  auto_ptr<ifstream> shortUsageMessage =
    da::openFileForReading(shortUsageFilepath, da::dieOnError);
  da::copyStream(*shortUsageMessage, out, da::dieOnError);

  if (verboseUsageP) {
    out << "\nUse the \"-h\" option to get a terser usage message than this"
      " one.\n\n";

    const char* verboseUsageFilepath =
      pteJoinPath(pathToExecutableDir(), "verboseUsageMessage.txt");
    auto_ptr<ifstream> verboseUsageMessage = 
      da::openFileForReading(verboseUsageFilepath, da::dieOnError);
    da::copyStream(*verboseUsageMessage, out, da::dieOnError);
  } else {
    out <<  "\nUse the \"--help\" option to get a longer usage message.\n";
  }
  if (exitWithFailureP) exit(EXIT_FAILURE);
  else exit(0);
}


local proc void
usage(bool exitWithFailureP=true)
{
  ::usage1(exitWithFailureP, false);
}


local proc void
verboseUsage()
{
  ::usage1(false, true);
}


local proc void
checkEndptr(const char* endptr)
{
  if (*endptr != '\0') usage();
}


//-----------------------------------------------------------------------------
// debugPrint(): local macro
//-----------------------------------------------------------------------------

#define debugPrint(message) \
   { if (Cl::getDebugLevel()) { \
        cerr << message << endl; \
     } \
   }


//=============================================================================
// CommandLineParser: local class
//=============================================================================

namespace {
  class CommandLineParser {

    // Private class variables:
    static const char*    _cv_inputFilepath;
    static const char*    _cv_outputFilepath;
    static int		  _cv_debugLevel;
    static bool           _cv_coerceToShorts;
    static bool           _cv_coerceToUnsignedShorts;
    static bool		  _cv_flipvFlag;
    static bool		  _cv_binomialBlurFlag;
    static bool           _cv_derivativeImageFilterFlag;
    static bool           _cv_flipImageFlag;
    static bool		  _cv_identityFlipFlag;

    // Private non-virtual methods:
    static void        parseExtendedOption(const char* const option);

    // Deactivate constructor:
    CommandLineParser();

  public:
    
    // Class methods:
    static void        parseCommandLine(const int argc,
					const char* const argv[]);

    // Class accessors:
    static const char* getInputFilepath()  { return _cv_inputFilepath;}
    static const char* getOutputFilepath() { return _cv_outputFilepath; }
    static int         getDebugLevel() { return _cv_debugLevel; }
    static bool        getCoerceToShorts() { return _cv_coerceToShorts; }
    static bool        getCoerceToUnsignedShorts()
                          { return _cv_coerceToUnsignedShorts; }
    static bool	       getBinomialBlurFlag()
                          { return _cv_binomialBlurFlag; }
    static bool        getDerivativeImageFilterFlag()
                          { return _cv_derivativeImageFilterFlag; }
    static bool        getFlipImageFlag() { return _cv_flipImageFlag; }
    static bool	       getFlipvFlag() { return _cv_flipvFlag; }
    static bool	       getIdentityFlipFlag() { return _cv_identityFlipFlag; }
  };

  const char* CommandLineParser::_cv_inputFilepath = 0;
  const char* CommandLineParser::_cv_outputFilepath = 0;
  int  	      CommandLineParser::_cv_debugLevel = 0;
  bool	      CommandLineParser::_cv_coerceToShorts = false;
  bool        CommandLineParser::_cv_coerceToUnsignedShorts = false;
  bool	      CommandLineParser::_cv_flipvFlag = false;
  bool	      CommandLineParser::_cv_binomialBlurFlag = false;
  bool        CommandLineParser::_cv_derivativeImageFilterFlag = false;
  bool        CommandLineParser::_cv_flipImageFlag = false;
  bool	      CommandLineParser::_cv_identityFlipFlag = false;

  typedef     CommandLineParser  Cl;

} // namespace


//-----------------------------------------------------------------------------
// parseCommandLine(): class method
//-----------------------------------------------------------------------------

method void CommandLineParser::
parseCommandLine(const int argc, const char* const argv[])
{
  ::opterr = true;
  char* endptr;
  
  char optionChar;
  const int cBase10 = 10;
  string extendedOption;
  int verboseHelpFlag = false;
  int flipvFlag = false;
  int noWcsFlag = false;
  int ripOrientationFlag = false;
  int rotateSkyFlag = false;
  int scaleDecFlag = false;
  int typicalFlag = false;

  // Specify the allowed options:
  const char shortopts[] = "Aa:D:fhN:o:Rr:Ss:Uv:";
  struct option longopts[] = {
    { "help", no_argument, &verboseHelpFlag, true },
    { "flipv", no_argument, &flipvFlag, true },
    { "no-wcs", no_argument, &noWcsFlag, true},
    { "rotate-sky", required_argument, &rotateSkyFlag,
      true },
    { "RIP", no_argument, &ripOrientationFlag, true },
    { "scale-dec", required_argument, &scaleDecFlag, true },
    { "typical", no_argument, &typicalFlag, true },
    { null, 0, null, 0 }
  };


  // Do the parsing:
  while ((optionChar = ::getopt_long(argc,
				     const_cast<char**>(argv),
				     shortopts,
				     longopts,
				     null))
         != -1)
    {
      switch (optionChar) {

      case 'A':
	itk::FITSImageIO::SetAutoScaleVelocityAxis(true);
	break;

      case 'a':
	itk::FITSImageIO::SetScaleAllAxes(strtod(optarg, &endptr));
	::checkEndptr(endptr);
	break;

// 	// TODO: Implement this
//       case 'd':
// 	// itk::FITSImageIO::SetScaleDec(strtod(optarg, &endptr));
// 	// ::checkEndptr(endptr);
// 	break;

      case 'D':
	_cv_debugLevel = strtol(optarg, null, cBase10);
	itk::FITSImageIO::SetDebugLevel(_cv_debugLevel);
	break;

      case 'h': ::usage(false);

      case 'N':
	itk::FITSImageIO::SetNullValue(strtod(optarg, &endptr));
	::checkEndptr(endptr);
	break;

      case 'o':
	Cl::parseExtendedOption(optarg);
	break;

      case 'r':
	itk::FITSImageIO::SetScaleRA(strtod(optarg, &endptr));
	::checkEndptr(endptr);
	break;

      case 'S':
        _cv_coerceToShorts = true;
        break;

      case 's':
        itk::FITSImageIO::SetScaleVoxelValues(strtod(optarg, &endptr));
	::checkEndptr(endptr);
        break;

      case 'U':
        _cv_coerceToUnsignedShorts = true;
        break;

      case 'v':
        itk::FITSImageIO::SetScaleVelocity(strtod(optarg, &endptr));
	::checkEndptr(endptr);
        break;


      case '?': 
	// On BSD (e.g., OS X), we end up here for unknown long options.
	// In this case, it seems next to impossible to tell the difference
	// between this case and the case where the user specifies "-?" the
	// command line.  Consequently, we should always consider it to be
	// a parsing error.
	::usage();

      case ':': 
	// We might end up here for options that should have had arguments
	// but didn't, depending on the platform.
	::usage();

      case 0:
	// Parse long options:
	if (verboseHelpFlag) {
	  ::verboseUsage();
	} else if (flipvFlag) {
	  flipvFlag = false;
	  _cv_flipvFlag = true;
	} else if (noWcsFlag) {
	  noWcsFlag = false;
	  itk::FITSImageIO::SetSuppressWCS(true);
	} else if (ripOrientationFlag) {
	  ripOrientationFlag = false;
	  itk::FITSImageIO::SetRIPOrientation(true);
	  itk::FITSImageIO::SetSuppressWCS(true);
	} else if (rotateSkyFlag) {
	  rotateSkyFlag = false;
	  itk::FITSImageIO::SetRotateSky(strtod(optarg, &endptr));
	  ::checkEndptr(endptr);
	} else if (scaleDecFlag) {
	  scaleDecFlag = false;
	  itk::FITSImageIO::SetScaleDec(strtod(optarg, &endptr));
	  ::checkEndptr(endptr);
	} else if (typicalFlag) {
	  typicalFlag = false;
	  itk::FITSImageIO::SetAutoScaleVelocityAxis(true);
	  itk::FITSImageIO::SetScaleAllAxes(1000);
	  itk::FITSImageIO::SetScaleRA(-1);
	  itk::FITSImageIO::SetScaleVoxelValues(1000);
	} else {
	  ::usage();
	}
	break;

      default:
        cerr << "Bug!" << endl;
        assert(false);
      } // switch
    } // while

  // Parse the command line positional arguments:
  if (argc - ::optind != 2) ::usage();
  _cv_inputFilepath = argv[::optind];
  _cv_outputFilepath = argv[::optind + 1];

  // Do some sanity checking to make sure that the options specified are
  // consistent with each other:
  if (itk::FITSImageIO::GetSuppressWCS() and
      itk::FITSImageIO::GetAutoScaleVelocityAxis())
    {
      da::warning("Velocity axis auto-scaling does not work when WCS\n"
		  "     is suppressed.");
    }
}


//-----------------------------------------------------------------------------
// parseExtendedOption(): private non-virtual method
//-----------------------------------------------------------------------------

// TODO: Replace this stuff with long options, which I have started using since
// originally implementing this.

method void CommandLineParser::
parseExtendedOption(const char* const option)
{
  debugPrint("extendedOption=" << option);
  string optionStr = option;
  if (optionStr == "derivativeImageFilter") {
    _cv_derivativeImageFilterFlag = true;
  } else if (optionStr == "flipImage") {
    _cv_flipImageFlag = true;
  } else if (optionStr == "binomialBlur") {
    _cv_binomialBlurFlag = true;
  } else if (optionStr == "suppressMetaDataDictionary") {
    itk::FITSImageIO::SetSuppressMetaDataDictionary(true);
  } else if (optionStr == "identityFlip") {
    _cv_identityFlipFlag = true;
  } else ::usage();
}


//=============================================================================
// Local functions
//=============================================================================

namespace {


//-----------------------------------------------------------------------------
// applyMeanFilter(): local template function
//-----------------------------------------------------------------------------

// template <class PixelType>
// local proc void
// doMeanFilter(Image<PixelType, c_dims>& image)
// {
// //     typedef itk::MeanImageFilter<ImageType, ImageType> FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     typename ImageType::SizeType indexRadius;
// //     indexRadius[0] = 5;
// //     indexRadius[1] = 5;
// //     indexRadius[2] = 5;
// //     filter->SetRadius(indexRadius);


//-----------------------------------------------------------------------------
// flipImage(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc typename Image<PixelType, c_dims>::Pointer
flipImage(const typename Image<PixelType, c_dims>::Pointer& image)
{
  typedef Image<PixelType, c_dims> ImageType;
  typedef itk::FlipImageFilter<ImageType> FilterType;
  typedef typename FilterType::FlipAxesArrayType FlipAxesArrayType;
  typename FilterType::Pointer filter = FilterType::New();
  FlipAxesArrayType flipArray;
  flipArray[0] = 0;
  flipArray[1] = 1;
  flipArray[2] = 0;
  filter->SetFlipAxes(flipArray);
  filter->SetInput(image);
  filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// applyBinomialBlurFilter(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc typename Image<PixelType, c_dims>::Pointer
applyBinomialBlur(const typename Image<PixelType, c_dims>::Pointer& image)
{
  typedef Image<PixelType, c_dims> ImageType;
  typedef itk::BinomialBlurImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetRepetitions(1);
  filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// flipv(): local template function
//-----------------------------------------------------------------------------

// This functions flips the image around the velocity axis.  It does this by
// actually moving the pixels around, rather than by messing with the direction
// cosines.

template <class PixelType>
local proc void
flipv(const typename Image<PixelType, c_dims>::Pointer& image)
{
  // CAVEAT: This code will break if c_dims ever changes from 3.

  image->Update();

  typedef Image<PixelType, c_dims> ImageType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType SizeType;

  typename ImageType::RegionType allOfImage =
    image->GetLargestPossibleRegion();
  typename ImageType::SizeType imageSize = allOfImage.GetSize();
  IndexType imageOrigin = allOfImage.GetIndex();

  const size_t c_i    = itk::FITSImageIO::c_i;
  const size_t c_j    = itk::FITSImageIO::c_j;
  const size_t c_k    = itk::FITSImageIO::c_k;

  const size_t min_ra_index  = imageOrigin[c_i];
  const size_t min_dec_index = imageOrigin[c_j];
  const size_t min_vel_index = imageOrigin[c_k];

  const size_t max_ra_index = min_ra_index + imageSize[c_i];
  const size_t max_dec_index = min_dec_index + imageSize[c_j];
  const size_t max_vel_index = min_vel_index + imageSize[c_k];
  const size_t midway_vel_index = min_vel_index + (imageSize[c_k] / 2);

  for (size_t ra_i = min_ra_index; ra_i < max_ra_index; ++ra_i) {
    for (size_t dec_i= min_dec_index; dec_i < max_dec_index; ++dec_i) {
      for (size_t vel_i          = min_vel_index,
	          opposite_vel_i = max_vel_index - 1;
	   vel_i <= midway_vel_index;
	   ++vel_i, --opposite_vel_i)
	{
	  IndexType thisPixelIndex;
	  IndexType oppositePixelIndex;
	  thisPixelIndex[c_i] = oppositePixelIndex[c_i] = ra_i;
	  thisPixelIndex[c_j] = oppositePixelIndex[c_j] = dec_i;
	  thisPixelIndex[c_k] = vel_i;
	  oppositePixelIndex[c_k] = opposite_vel_i;
	  PixelType tmp = image->GetPixel(thisPixelIndex);
	  image->SetPixel(thisPixelIndex, image->GetPixel(oppositePixelIndex));
	  image->SetPixel(oppositePixelIndex, tmp);
	}
    }
  }
}


//-----------------------------------------------------------------------------
// convertInputFileToItkFile(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc int
convertInputFileToItkFile(const char* const inputFilepath,
			  const char* const outputFilepath)
{
  typedef itk::Image<PixelType, c_dims> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  typename WriterType::Pointer writer = WriterType::New();
  reader->SetFileName(inputFilepath);
  typename ImageType::Pointer image = reader->GetOutput();
  if (Cl::getFlipImageFlag()) image = ::flipImage<PixelType>(image);
  if (Cl::getBinomialBlurFlag()) image = ::applyBinomialBlur<PixelType>(image);
  if (Cl::getIdentityFlipFlag()) {
    image = ::flipImage<PixelType>(image);
    image = ::flipImage<PixelType>(image);
  }
  if (Cl::getFlipvFlag()) ::flipv<PixelType>(image);
  writer->SetInput(image);
  writer->SetFileName(outputFilepath);
  writer->Update();
  return EXIT_SUCCESS;

  // Note: The above call to Update() might raise an execption.  If you were to
  // want to catch this exception, here is now you might do it.  You don't
  // really need to, though, as the default exception handler will output a
  // message that is not particularly more cryptic than the following:
  //
  //   try {
  //     writer->Update(); 
  //   } catch(itk::ExceptionObject& err) { 
  //     cerr << "ExceptionObject caught !" << endl; 
  //     cerr << err << endl; 
  //     return EXIT_FAILURE;
  //   } 
}


//-----------------------------------------------------------------------------
} // END local functions
//-----------------------------------------------------------------------------


//=============================================================================
// main()
//=============================================================================

proc int
main(const int argc, const char* const argv[])
{
  ::setArgv(argc, argv);
  Cl::parseCommandLine(argc, argv);

  // Register FITS one factory with the ImageIOFactory.
  itk::FITSImageIOFactory::RegisterOneFactory();

  int status = -666;   // If the following code is correct, this value will
                       // always get overwritten.
  if (Cl::getCoerceToShorts()) {
    status = convertInputFileToItkFile<short>(
	        Cl::getInputFilepath(),
	        Cl::getOutputFilepath());
  } else if (Cl::getCoerceToUnsignedShorts()) {
    status = convertInputFileToItkFile<unsigned short>(
	        Cl::getInputFilepath(),
		Cl::getOutputFilepath());
  } else {
    status = convertInputFileToItkFile<float>(
                Cl::getInputFilepath(),
		Cl::getOutputFilepath());
  }
  return status;
}
