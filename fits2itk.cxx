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

#include <dlfcn.h>                         // For dlopen()
#include <libgen.h>			   // For basename()
#include <getopt.h>
#include <stdlib.h>			   // For getenv(), setenv()
#include <sys/param.h>                     // For MAXPATHLEN
#include <cassert>
#include <cmath>

#include <fstream>
using std::ifstream;
using std::istream;
using std::ostream;
using std::ios;

#include <memory>
using std::auto_ptr;

#include <string>
using std::string;

#include <itkImage.h>
using itk::Image;

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkFITSImageIO.h>
using itk::FITSImageIO;

#include <itkFITSWCSTransform.h>
using itk::FITSWCSTransform;

#include <itkFITSImageUtils.h>
using itk::fits::applyFlipImageFilter;
using itk::fits::applyBinomialBlurFilter;
using itk::fits::reflectPixels;
using itk::fits::initializeChangeOfBasis;
using itk::fits::transformToNorthOrientedEquiangular;
using itk::fits::transformToUnreorientedEquiangular;
using itk::fits::reorientNorth;

#include <pathToExecutable.h>
#include <da_util.h>
#include <da_sugar.h>

extern const char fits2itkVersion[];
const int c_dims = FITSImageIO::c_dims;

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
  usage1(exitWithFailureP, false);
}


local proc void
verboseUsage()
{
  usage1(false, true);
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

    static bool           _cv_coerceToShorts;
    static bool           _cv_coerceToUnsignedShorts;
    static int		  _cv_debugLevel;
    static bool		  _cv_dontWrite;
    static bool		  _cv_flipDecFlag;
    static bool		  _cv_flipRAFlag;
    static bool		  _cv_flipVFlag;
    static bool		  _cv_reorientNorth;
    static bool		  _cv_transformToEquiangular;

    static bool		  _cv_binomialBlurFlag;
    static bool           _cv_derivativeImageFilterFlag;
    static bool           _cv_flipImageFilterFlag;
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
    static bool        getCoerceToShorts() { return _cv_coerceToShorts; }
    static bool        getCoerceToUnsignedShorts()
                          { return _cv_coerceToUnsignedShorts; }
    static int         getDebugLevel() { return _cv_debugLevel; }
    static bool	       getDontWrite() { return _cv_dontWrite; }
    static bool	       getFlipDecFlag() { return _cv_flipDecFlag; }
    static bool        getFlipRAFlag() { return _cv_flipRAFlag; }
    static bool	       getFlipVFlag() { return _cv_flipVFlag; }
    static bool	       getReorientNorth() { return _cv_reorientNorth; }
    static bool	       getTransformToEquiangular()
                          { return _cv_transformToEquiangular; }

    static bool	       getBinomialBlurFlag()
                          { return _cv_binomialBlurFlag; }
    static bool        getDerivativeImageFilterFlag()
                          { return _cv_derivativeImageFilterFlag; }
    static bool        getFlipImageFilterFlag()
                          { return _cv_flipImageFilterFlag; }
    static bool	       getIdentityFlipFlag() { return _cv_identityFlipFlag; }
  };

  const char* CommandLineParser::_cv_inputFilepath = 0;
  const char* CommandLineParser::_cv_outputFilepath = 0;
  bool	      CommandLineParser::_cv_coerceToShorts = false;
  bool        CommandLineParser::_cv_coerceToUnsignedShorts = false;
  int  	      CommandLineParser::_cv_debugLevel = 0;
  bool	      CommandLineParser::_cv_dontWrite = false;
  bool	      CommandLineParser::_cv_flipDecFlag = false;
  bool	      CommandLineParser::_cv_flipRAFlag = false;
  bool	      CommandLineParser::_cv_flipVFlag = false;
  bool	      CommandLineParser::_cv_reorientNorth = false;
  bool	      CommandLineParser::_cv_transformToEquiangular = false;

  bool	      CommandLineParser::_cv_binomialBlurFlag = false;
  bool        CommandLineParser::_cv_derivativeImageFilterFlag = false;
  bool        CommandLineParser::_cv_flipImageFilterFlag = false;
  bool	      CommandLineParser::_cv_identityFlipFlag = false;

  typedef     CommandLineParser  Cl;

} // namespace


//-----------------------------------------------------------------------------
// parseCommandLine(): class method
//-----------------------------------------------------------------------------

method void CommandLineParser::
parseCommandLine(const int argc, const char* const argv[])
{
  opterr = true;
  char* endptr;
  
  char optionChar;
  const int cBase10 = 10;
  string extendedOption;
  int verboseHelpFlag = false;
  int equiangularFlag = false;
  int flipDecFlag = false;
  int flipRAFlag = false;
  int flipVFlag = false;
  int noWcsFlag = false;
  int reorientNorthFlag = false;
  int ripOrientationFlag = false;
  int rotateSkyFlag = false;
  int scaleDecFlag = false;
  int typicalFlag = false;
  int verboseFlag = false;

  // Specify the allowed options:
  const char shortopts[] = "Aa:D:fhnN:o:Rr:Ss:Uv:";

  struct option longopts[] = {
    { "equiangular", no_argument, &equiangularFlag, true },
    { "help", no_argument, &verboseHelpFlag, true },
    { "flip-dec", no_argument, &flipDecFlag, true },
    { "flip-ra", no_argument, &flipRAFlag, true },
    { "flip-v", no_argument, &flipVFlag, true },
    { "no-wcs", no_argument, &noWcsFlag, true},
    { "reorient-north", no_argument, &reorientNorthFlag, true },
    { "rotate-sky", required_argument, &rotateSkyFlag,
      true },
    { "RIP", no_argument, &ripOrientationFlag, true },
    { "scale-dec", required_argument, &scaleDecFlag, true },
    { "typical", no_argument, &typicalFlag, true },
    { "verbose", no_argument, &verboseFlag, true },
    { null, 0, null, 0 }
  };


  // Do the parsing:
  while ((optionChar = getopt_long(argc,
				   const_cast<char**>(argv),
				   shortopts,
				   longopts,
				   null))
         != -1)
    {
      switch (optionChar) {

      case 'A':
	//d FITSImageIO::SetAutoScaleVelocityAxis(true);
	break;

      case 'a':
	//d FITSImageIO::SetScaleAllAxes(strtod(optarg, &endptr));
	::checkEndptr(endptr);
	break;

// 	// TODO: Implement this
//       case 'd':
// 	// FITSImageIO::SetScaleDec(strtod(optarg, &endptr));
// 	// ::checkEndptr(endptr);
// 	break;

      case 'D':
	_cv_debugLevel = strtol(optarg, null, cBase10);
	//d FITSImageIO::SetDebugLevel(_cv_debugLevel);
	break;

      case 'h': usage(false);

      case 'N':
	//d FITSImageIO::SetNullValue(strtod(optarg, &endptr));
	checkEndptr(endptr);
	break;

      case 'n':
	_cv_dontWrite = true;
	break;

      case 'o':
	Cl::parseExtendedOption(optarg);
	break;

      case 'r':
	//d FITSImageIO::SetScaleRA(strtod(optarg, &endptr));
	checkEndptr(endptr);
	break;

      case 'S':
        _cv_coerceToShorts = true;
        break;

      case 's':
        //d FITSImageIO::SetScaleVoxelValues(strtod(optarg, &endptr));
	checkEndptr(endptr);
        break;

      case 'U':
        _cv_coerceToUnsignedShorts = true;
        break;

      case 'v':
        //d FITSImageIO::SetScaleVelocity(strtod(optarg, &endptr));
	checkEndptr(endptr);
        break;


      case '?': 
	// On BSD (e.g., OS X), we end up here for unknown long options.
	// In this case, it seems next to impossible to tell the difference
	// between this case and the case where the user specifies "-?" the
	// command line.  Consequently, we should always consider it to be
	// a parsing error.
	usage();

      case ':': 
	// We might end up here for options that should have had arguments
	// but didn't, depending on the platform.
	usage();

      case 0:
	// Parse long options:
	if (verboseHelpFlag) {
	  verboseUsage();
	} else if (equiangularFlag) {
	  equiangularFlag =  false;
	  _cv_transformToEquiangular = true;
	} else if (flipDecFlag) {
	  flipDecFlag = false;
	  _cv_flipDecFlag = true;
	} else if (flipRAFlag) {
	  flipRAFlag = false;
	  _cv_flipRAFlag = true;
	} else if (flipVFlag) {
	  flipVFlag = false;
	  _cv_flipVFlag = true;
	} else if (noWcsFlag) {
	  noWcsFlag = false;
	  //d FITSImageIO::SetSuppressWCS(true);
	} else if (reorientNorthFlag) {
	  reorientNorthFlag = false;
	  _cv_reorientNorth = true;
	} else if (ripOrientationFlag) {
	  ripOrientationFlag = false;
	  //d FITSImageIO::SetRIPOrientation(true);
	  //d FITSImageIO::SetSuppressWCS(true);
	} else if (rotateSkyFlag) {
	  rotateSkyFlag = false;
	  //d FITSImageIO::SetRotateSky(strtod(optarg, &endptr));
	  checkEndptr(endptr);
	} else if (scaleDecFlag) {
	  scaleDecFlag = false;
	  //d FITSImageIO::SetScaleDec(strtod(optarg, &endptr));
	  checkEndptr(endptr);
	} else if (typicalFlag) {
	  typicalFlag = false;
	  //d FITSImageIO::SetAutoScaleVelocityAxis(true);
	  //d FITSImageIO::SetScaleAllAxes(1000);
	  //d FITSImageIO::SetScaleRA(-1);
	  //d FITSImageIO::SetScaleVoxelValues(1000);
	} else if (verboseFlag) {
	  verboseFlag = false;
	  //d FITSImageIO::SetVerbose(true);
	} else {
	  usage();
	}
	break;

      default:
        cerr << "Bug!" << endl;
        assert(false);
      } // switch
    } // while

  // Parse the command line positional arguments:
  if (_cv_dontWrite) {
    if (argc - ::optind != 1) usage();
  } else {
    if (argc - ::optind != 2) usage();
    _cv_outputFilepath = argv[::optind + 1];
  }
  _cv_inputFilepath = argv[::optind];

  // Do some sanity checking to make sure that the options specified are
  // consistent with each other:
//d   if (FITSImageIO::GetSuppressWCS() and
//d       FITSImageIO::GetAutoScaleVelocityAxis())
//d     {
//d       da::warning("Velocity axis auto-scaling does not work when WCS\n"
//d 		  "     is suppressed.");
//d     }
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
  } else if (optionStr == "flipImageFilter") {
    _cv_flipImageFilterFlag = true;
  } else if (optionStr == "binomialBlur") {
    _cv_binomialBlurFlag = true;
  } else if (optionStr == "suppressMetaDataDictionary") {
    //d FITSImageIO::SetSuppressMetaDataDictionary(true);
  } else if (optionStr == "identityFlip") {
    _cv_identityFlipFlag = true;
  } else usage();
}


//=============================================================================
// Local functions
//=============================================================================

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

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFilepath);
  typename ImageType::Pointer image = reader->GetOutput();
  if (Cl::getFlipImageFilterFlag()) {
    image = applyFlipImageFilter<PixelType>(image);
  }
  if (Cl::getBinomialBlurFlag()) {
    image = applyBinomialBlurFilter<PixelType>(image);
  }
  if (Cl::getIdentityFlipFlag()) {
    image = applyFlipImageFilter<PixelType>(image);
    image = applyFlipImageFilter<PixelType>(image);
  }
  reflectPixels<PixelType>(*image,
			     Cl::getFlipRAFlag(),
			     Cl::getFlipDecFlag(),
			     Cl::getFlipVFlag());


  if (Cl::getReorientNorth() or Cl::getTransformToEquiangular()) {
    // TODO: Figure out how to do this without reading in the entire image.
    reader->Update();
    initializeChangeOfBasis(*image);

    if (Cl::getTransformToEquiangular()) {
      if (Cl::getReorientNorth()) {
	transformToNorthOrientedEquiangular(*image);
      } else {
	transformToUnreorientedEquiangular(*image);
      }
    } else if (Cl::getReorientNorth()) {
      reorientNorth(*image);
    }
  }

  if (outputFilepath) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(image);
    writer->SetFileName(outputFilepath);
    writer->Update();
  } else {
    reader->Update();
  }

  //d if (FITSImageIO::GetVerbose()) ::writeImageInfo(*image, cout);
  
  return EXIT_SUCCESS;

  // Note: The above call to writer->Update() might raise an execption.  If we
  // were to want to catch this exception, here is how we might do it.  You
  // don't really need to, though, as the default exception handler will output
  // a message that is not particularly more cryptic than the following:
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
// END local functions
//-----------------------------------------------------------------------------


//=============================================================================
// main()
//=============================================================================

proc int
main(const int argc, const char* const argv[])
{
  setArgv(argc, argv);

  // Set the environment variable $ITK_AUTOLOAD_PATH so that
  // ObjectFactoryBase::LoadLibrariesInPath() can find FITSImageIOFactory:
  {
#ifdef _WIN32
    const char pathSeparator = ';';
#else
    const char pathSeparator = ':';
#endif
  
    const char* oldAutoloadPath = getenv("ITK_AUTOLOAD_PATH");
    const string newAutoloadPath = oldAutoloadPath
      ? string(pathToExecutableDir()) + pathSeparator + oldAutoloadPath
      : string(pathToExecutableDir());
    setenv("ITK_AUTOLOAD_PATH", newAutoloadPath.c_str(), true);
  }

  Cl::parseCommandLine(argc, argv);

  // This is how we used to register FITSImageIOFactory, before we changed to
  //dynamic loading:
  //
  // FITSImageIOFactory::RegisterOneFactory();

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

