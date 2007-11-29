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
using itk::fits::initializeChangeOfBasis;
using itk::fits::reflectPixels;
using itk::fits::reorientNorth;
using itk::fits::setNullValue;
using itk::fits::transformToNorthOrientedEquiangular;
using itk::fits::transformToUnreorientedEquiangular;
using itk::fits::writeImageInfo;

#include <pathToExecutable.h>
#include <da_util.h>

#include <da_sugar.h>
using da::setDebugLevel;

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


//=============================================================================
// CommandLineParser: local class
//=============================================================================

// BEGIN local namespace
namespace {

  class CommandLineParser {

    // Instance variables:
    const char*    _inputFilepath;
    const char*    _outputFilepath;

    bool           _coerceToShorts;
    bool           _coerceToUnsignedShorts;
    int		   _debugLevel;
    bool	   _dontWrite;
    bool	   _flipDecFlag;
    bool	   _flipRAFlag;
    bool	   _flipVFlag;
    double	   _nullValue;
    bool	   _reorientNorth;
    bool	   _transformToEquiangular;

    bool	   _binomialBlurFlag;
    bool           _derivativeImageFilterFlag;
    bool           _flipImageFilterFlag;
    bool	   _identityFlipFlag;

    // Private non-virtual methods:
    void        parseExtendedOption(const char* const option);

  public:
    
    // Constructor:
    CommandLineParser(const int argc, const char* const argv[]);

    // Accessor methods:
    const char* getInputFilepath() const { return _inputFilepath;}
    const char* getOutputFilepath() const { return _outputFilepath; }
    bool        getCoerceToShorts() const { return _coerceToShorts; }
    bool        getCoerceToUnsignedShorts() const
                          { return _coerceToUnsignedShorts; }
    int         getDebugLevel() const { return _debugLevel; }
    bool	getDontWrite() const { return _dontWrite; }
    bool	getFlipDecFlag() const { return _flipDecFlag; }
    bool        getFlipRAFlag() const { return _flipRAFlag; }
    bool	getFlipVFlag() const { return _flipVFlag; }
    double	getNullValue() const { return _nullValue; }
    bool	getReorientNorth() const { return _reorientNorth; }
    bool	getTransformToEquiangular() const
                          { return _transformToEquiangular; }
    bool	getBinomialBlurFlag() const
                          { return _binomialBlurFlag; }
    bool        getDerivativeImageFilterFlag() const
                          { return _derivativeImageFilterFlag; }
    bool        getFlipImageFilterFlag() const
                          { return _flipImageFilterFlag; }
    bool	getIdentityFlipFlag() const { return _identityFlipFlag; }
  };


ctor
CommandLineParser::CommandLineParser(const int argc, const char* const argv[])
  : _inputFilepath(0),
    _outputFilepath(0),
    _coerceToShorts(false),
    _coerceToUnsignedShorts(false),
    _debugLevel(0),
    _dontWrite(false),
    _flipDecFlag(false),
    _flipRAFlag(false),
    _flipVFlag(false),
    _nullValue(0.0),
    _reorientNorth(false),
    _transformToEquiangular(false),
    _binomialBlurFlag(false),
    _derivativeImageFilterFlag(false),
    _flipImageFilterFlag(false),
    _identityFlipFlag(false)
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

  // Specify the allowed short options:
  const char shortopts[] = "a:D:fhnN:o:Rr:Ss:Uv:";

  // Specify the allowed long options:
  struct option longopts[] = {

    { "equiangular", no_argument, &equiangularFlag, true },

    { "help", no_argument, &verboseHelpFlag, true },

    { "flip-dec", no_argument, &flipDecFlag, true },

    { "flip-ra", no_argument, &flipRAFlag, true },

    { "flip-v", no_argument, &flipVFlag, true },

    { "no-wcs", no_argument, &noWcsFlag, true},
    { "nw", no_argument, &noWcsFlag, true},

    { "null-value", required_argument, 0, 'N'},

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

      case 'a':
	// XXX: You need to figure out a new way to scale all axes.

	//d FITSImageIO::SetScaleAllAxes(strtod(optarg, &endptr));
	::checkEndptr(endptr);
	break;

      case 'D':
	setDebugLevel(strtol(optarg, null, cBase10));
	break;

      case 'h': usage(false);

      case 'N':
	_nullValue = strtod(optarg, &endptr);
	checkEndptr(endptr);
	break;

      case 'n':
	_dontWrite = true;
	break;

      case 'o':
	this->parseExtendedOption(optarg);
	break;

      case 'r':
	//d FITSImageIO::SetScaleRA(strtod(optarg, &endptr));
	checkEndptr(endptr);
	break;

      case 'S':
        _coerceToShorts = true;
        break;

      case 's':
        //d FITSImageIO::SetScaleVoxelValues(strtod(optarg, &endptr));
	checkEndptr(endptr);
        break;

      case 'U':
        _coerceToUnsignedShorts = true;
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
	  _transformToEquiangular = true;
	} else if (flipDecFlag) {
	  flipDecFlag = false;
	  _flipDecFlag = true;
	} else if (flipRAFlag) {
	  flipRAFlag = false;
	  _flipRAFlag = true;
	} else if (flipVFlag) {
	  flipVFlag = false;
	  _flipVFlag = true;
	} else if (noWcsFlag) {
	  noWcsFlag = false;
	  //d FITSImageIO::SetSuppressWCS(true);
	} else if (reorientNorthFlag) {
	  reorientNorthFlag = false;
	  _reorientNorth = true;
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
	  da::setVerbosityLevel(1);
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
  if (_dontWrite) {
    if (argc - ::optind != 1) usage();
  } else {
    if (argc - ::optind != 2) usage();
    _outputFilepath = argv[::optind + 1];
  }
  _inputFilepath = argv[::optind];

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
    _derivativeImageFilterFlag = true;
  } else if (optionStr == "flipImageFilter") {
    _flipImageFilterFlag = true;
  } else if (optionStr == "binomialBlur") {
    _binomialBlurFlag = true;
  } else if (optionStr == "suppressMetaDataDictionary") {
    //d FITSImageIO::SetSuppressMetaDataDictionary(true);
  } else if (optionStr == "identityFlip") {
    _identityFlipFlag = true;
  } else usage();
}


} // END local namespace


//=============================================================================
// Local functions
//=============================================================================


//-----------------------------------------------------------------------------
// showFactoryClasses(): local function
//-----------------------------------------------------------------------------

// This function was used for debugging.

// proc void
// showFactoryClasses()
// {
//   cout << "Registered ITK factory classes: ";
//   typedef list<ObjectFactoryBase*> List;
//   typedef List::iterator Iter;
//   List factories = ObjectFactoryBase::GetRegisteredFactories();
//   for (Iter factoryIter = factories.begin();
//        factoryIter != factories.end();
//        ++factoryIter)
//     {
//       const char* const nameOfClass = (*factoryIter)->GetNameOfClass();
//       if (strcmp(nameOfClass, "FITSImageIOFactory") == 0) {
// 	cout << " YAAAAAAY!!!! ";
// 	const FITSImageIOFactory* factory =
// 	  (const FITSImageIOFactory*) *factoryIter;
// 	factory->SetTestValue(10);
// 	factory->GetTestValue();
//       } else {
// 	cout << nameOfClass << " ";
//       }
//     }
//   cout << endl;
//   cout << "Test Value1=" << FITSImageIOFactory::_cv_testValue << endl;
// }


//-----------------------------------------------------------------------------
// convertInputFileToItkFile(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc int
convertInputFileToItkFile(CommandLineParser& cl)
{
  typedef itk::Image<PixelType, c_dims> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(cl.getInputFilepath());
  typename ImageType::Pointer image = reader->GetOutput();
  if (cl.getFlipImageFilterFlag()) {
    image = applyFlipImageFilter<PixelType>(image);
  }
  if (cl.getBinomialBlurFlag()) {
    image = applyBinomialBlurFilter<PixelType>(image);
  }
  if (cl.getIdentityFlipFlag()) {
    image = applyFlipImageFilter<PixelType>(image);
    image = applyFlipImageFilter<PixelType>(image);
  }
  reflectPixels<PixelType>(*image,
			     cl.getFlipRAFlag(),
			     cl.getFlipDecFlag(),
			     cl.getFlipVFlag());


  if (cl.getReorientNorth() or cl.getTransformToEquiangular()) {
    // TODO: Figure out how to do this without reading in the entire image.
    reader->Update();
    initializeChangeOfBasis(*image);

    if (cl.getTransformToEquiangular()) {
      if (cl.getReorientNorth()) {
	transformToNorthOrientedEquiangular(*image);
      } else {
	transformToUnreorientedEquiangular(*image);
      }
    } else if (cl.getReorientNorth()) {
      reorientNorth(*image);
    }
  }

  if (cl.getOutputFilepath()) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(image);
    writer->SetFileName(cl.getOutputFilepath());
    writer->Update();
  } else {
    reader->Update();
  }

  if (da::getVerbosityLevel()) writeImageInfo(*image, cout);

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
// handleOptions(): local function
//-----------------------------------------------------------------------------

local proc void
handleOptions(const CommandLineParser& cl)
{
  setNullValue(cl.getNullValue());
  da::setDebugLevel(cl.getDebugLevel());
}


//-----------------------------------------------------------------------------
// setItkAutoloadPath(): local function
//-----------------------------------------------------------------------------

// Sets the environment variable $ITK_AUTOLOAD_PATH so that
// ObjectFactoryBase::LoadLibrariesInPath() can find FITSImageIOFactory:

local proc void
setItkAutoloadPath()
{
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
  setItkAutoloadPath();
  CommandLineParser cl(argc, argv);
  handleOptions(cl);

  // This is how we used to register FITSImageIOFactory, before we changed to
  // dynamic loading:
  //
  // FITSImageIOFactory::RegisterOneFactory();

  int status = -666;   // If the following code is correct, this value will
                       // always get overwritten.
  if (cl.getCoerceToShorts()) {
    status = convertInputFileToItkFile<short>(cl);
  } else if (cl.getCoerceToUnsignedShorts()) {
    status = convertInputFileToItkFile<unsigned short>(cl);
  } else {
    status = convertInputFileToItkFile<float>(cl);
  }
  return status;
}
