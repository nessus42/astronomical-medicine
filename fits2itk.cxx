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

// #include <itkDerivativeImageFilter.h>
// #include <itkMeanImageFilter.h>
// #include <itkBinaryMedianImageFilter.h>
// #include <itkGradientAnisotropicDiffusionImageFilter.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkFlipImageFilter.h>
#include <itkBinomialBlurImageFilter.h>

#include <libgen.h>			   // For basename()
#include <getopt.h>

#include <cassert>
#include <string>

#include <itkFITSImageIOFactory.h>
#include <itkFITSImageIO.h>
#include <da_sugar.h>

using std::string;
using itk::Image;

extern const char fits2itkVersion[];

//-----------------------------------------------------------------------------
// Usage messages
//-----------------------------------------------------------------------------

const char shortUsageMessage[] =

"  fits2itk [-ASU] [-a axes-scale] [-D debug-level] [-N null-value]\n"
"           [-r RA-scale] [-s pixel-scale] [-v velocity-scale]\n"
"           input-file output-file\n"
"\n"
"  A: auto-scale velocity axis\n"
"  S: coerce pixel values to shorts\n"
"  U: coerce pixel values to unsigned shorts\n"
"\n"
"\"-a\" scales all the axes, while \"-v\" scales only the velocity\n"
"axis.  If both are used, then both scaling factors will be applied.\n"
"\n"
"fits2itk supports CFITSIO's \"extended filename syntax\".\n"
;


const char verboseUsageMessage[] =

"Typical Usage\n"
"-------------\n"

"The typical usage of fits2itk for converting FITS files into \"nrrd\" files\n"
"for use with 3D Slicer is as follows:\n"
"\n"
"   $ fits2itk -A -a 1000 -r -1 -s 1000 inputfile.fits outputfile.nrrd\n"
"\n"
"The above command auto-scales the velocity axis, and then scales all of the\n"
"axes by 1000.  The right ascension axis is then flipped in order to present\n"
"it in the orientation to which astronomers are accustomed.  The pixel\n"
"values of the image are then all multiplied by 1000.\n"
"\n"
"To save you a bit of typing, this can be abbreviated as\n"
"\n"
"   $ fits2itk --typical inputfile.fits outputfile.nrrd\n"
"\n"
"\n"
"Auto Scaling\n"
"------------\n"
"\n"
"The \"-A\" (velocity auto-scale) option works by fitting the velocity axis\n"
"into a cube that is defined by the larger of the two positional axes.  If\n"
"this option is specified in conjunction with other scaling options (such as\n"
"-v, for instance), then -A is applied first.  The other scaling option are\n"
"then also applied, multiplicatively.\n"
"\n"
"\n"
"Dealing with NaN's\n"
"------------------\n"
"\n"
"\"NaN\" means \"not a number\".  This is a special floating point value\n"
"that represents undefined values.  FITS images often use NaN's to represent\n"
"undefined pixels, but unfortunately not all software can handle NaN's,\n"
"including some versions of 3D Slicer.  If you have images that contain\n"
"NaN's and want to use them with software that doesn't understand NaN's, you\n"
"can tell fits2itk to change all of the NaN's to a floating point value of\n"
"your choice using the \"-N\" option.  (Note: The value you specify is not\n"
"allowed to be 0, but it can be 0.000001, or somesuch.)\n"
"\n"
"In addition to supporting undefined floating point values, FITS also has a\n"
"notion of an undefined value for a FITS image with integral pixel values.\n"
"Unlike for floating point values, however, there is no special integer\n"
"value that represents an undefined values.  Consequently, FITS allows the\n"
"creator of a FITS image to specify whatever integer he or she would like to\n"
"represent undefined pixels.  If, as it so happens, that you don't like the\n"
"specific integer chosen by the creator of a FITS image for representing\n"
"undefined pixels, you can remap the undefined value to a different integer\n"
"using the -N option of fits2itk.\n"
"\n"
"\n"
"CFITSIO Extended Filename Syntax\n"
"--------------------------------\n"
"\n"
"fits2itk supports CFITSIO's \"extended filename syntax\", which allows all\n"
"sorts of interesting things.  For example, if you have a data cube with an\n"
"extra 1-pixel-thick fourth dimension, you can slice off the extra dimension\n"
"like so:\n"
"\n"
"   $ fits2itk \"input.fits[*,*,*,1:1][col #NAXIS=3]\" output.nrrd\n"
"\n"
"You can also specify a URL, rather than a filename, for the input file and\n"
"the input file will be automatically fetched via HTTP.  For the complete\n"
"manual on the extended filename syntax, see the CFITSIO User's Reference\n"
"Guide chapter on it here:\n"
"\n"
"   http://heasarc.nasa.gov/docs/software/fitsio/c/c_user/node79.html\n"
"\n"
"\n"
"Debugging Output\n"
"----------------\n"
"\n"
"Various bits of information that are useful to the developer of fits2itk,\n"
"but which probably aren't of much interest to you, can be output to stderr\n"
"during a file conversion by specifying \"-D1\" on the fits2itk command\n"
"line.\n"
"\n"
"\n"
"This Help Text\n"
"--------------\n"
"\n"
"If this help text scrolled by too quickly for you to read, you might try\n"
"piping it through \"more\", like so:\n"
"\n"
"  $ fits2itk --help | more\n"
"\n"
"If you would like a terser help message than this one, use the \"-h\"\n"
"option instead of \"--help\"\n"
;

//-----------------------------------------------------------------------------
// usage()
//-----------------------------------------------------------------------------

local proc void
usage1(bool exitWithFailure, bool verboseUsageFlag)
{
//  if (daProgramName().size()) cerr << basename(daProgramName().c_str());
//  else cerr << "fits2itk";
  
  std::ostream& out = exitWithFailure ? cerr : cout;
   out << "fits2itk " << fits2itkVersion << "\n\n";
  out << "usage:\n";
  out << ::shortUsageMessage;
  
  if (verboseUsageFlag) {
    out << "\nUse the \"-h\" option to get a terser usage message than this"
      " one.\n\n";
    out << ::verboseUsageMessage;
  } else {
    out <<  "\nUse the \"--help\" option to get a longer usage message.\n";
  }
  if (exitWithFailure) exit(EXIT_FAILURE);
  else exit(0);
}


local proc void
usage(bool exitWithFailure=true)
{
  ::usage1(exitWithFailure, false);
}


local proc void
verboseUsage()
{
  ::usage1(false, true);
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
    static bool	       getIdentityFlipFlag() { return _cv_identityFlipFlag; }
  };

  const char* CommandLineParser::_cv_inputFilepath = 0;
  const char* CommandLineParser::_cv_outputFilepath = 0;
  int  	      CommandLineParser::_cv_debugLevel = 0;
  bool	      CommandLineParser::_cv_coerceToShorts = false;
  bool        CommandLineParser::_cv_coerceToUnsignedShorts = false;
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
  char optionChar;
  const int cBase10 = 10;
  string extendedOption;
  int verboseHelpFlag = false;
  int typicalFlag = false;

  // Specify the allowed options:
  const char shortopts[] = "Aa:D:fhN:o:Rr:Ss:Uv:";
  struct option longopts[] = {
    { "help", no_argument, &verboseHelpFlag, true },
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
	itk::FITSImageIO::SetScaleAllAxes(strtod(optarg, null));
	break;

      case 'D':
	_cv_debugLevel = strtol(optarg, null, cBase10);
	itk::FITSImageIO::SetDebugLevel(_cv_debugLevel);
	break;

      case 'h': ::usage(false);

      case 'N':
	itk::FITSImageIO::SetNullValue(strtod(optarg, null));
	break;

      case 'o':
	Cl::parseExtendedOption(optarg);
	break;

      case 'R':
	// TODO: This isn't implemented yet.
	itk::FITSImageIO::SetRotateSky(strtod(optarg, null));
	break;

      case 'r':
	itk::FITSImageIO::SetScaleRA(strtod(optarg, null));
	break;

      case 'S':
        _cv_coerceToShorts = true;
        break;

      case 's':
        itk::FITSImageIO::SetScaleVoxelValues(strtod(optarg, null));
        break;

      case 'U':
        _cv_coerceToUnsignedShorts = true;
        break;

      case 'v':
        // TODO: Make the following (and other similar code here) more
        // robust:
        itk::FITSImageIO::SetScaleVelocityAxis(strtod(optarg, null));
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
	} else if (typicalFlag) {
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
}


//-----------------------------------------------------------------------------
// parseExtendedOption(): private non-virtual method
//-----------------------------------------------------------------------------

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
    itk::FITSImageIO::SuppressMetaDataDictionary();
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
// doMeanFilter(Image<PixelType, 3>& image)
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
local proc typename Image<PixelType, 3>::Pointer
flipImage(const typename Image<PixelType, 3>::Pointer& image)
{
  typedef Image<PixelType, 3> ImageType;
  typedef itk::FlipImageFilter<ImageType> FilterType;
  typedef typename FilterType::FlipAxesArrayType FlipAxesArrayType;
  typename FilterType::Pointer filter = FilterType::New();
  FlipAxesArrayType flipArray;
  flipArray[0] = 1;
  flipArray[1] = 0;
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
local proc typename Image<PixelType, 3>::Pointer
applyBinomialBlur(const typename Image<PixelType, 3>::Pointer& image)
{
  typedef Image<PixelType, 3> ImageType;
  typedef itk::BinomialBlurImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetRepetitions(1);
  filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// convertInputFileToItkFile(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc int
convertInputFileToItkFile(const char* const inputFilepath,
			  const char* const outputFilepath)
{
  const unsigned Dimensions = 3;
  typedef itk::Image<PixelType, Dimensions> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  typename WriterType::Pointer writer = WriterType::New();
  reader->SetFileName(inputFilepath);
  typename ImageType::Pointer image = reader->GetOutput();
  if (Cl::getFlipImageFlag()) {
    image = ::flipImage<PixelType>(image);
  }
  if (Cl::getBinomialBlurFlag()) {
    image = ::applyBinomialBlur<PixelType>(image);
  }
  if (Cl::getIdentityFlipFlag()) {
    image = ::flipImage<PixelType>(image);
    image = ::flipImage<PixelType>(image);
  }
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
  ::daSetProgramName(argv[0]);
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
