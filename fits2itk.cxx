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
using itk::fits::FITSImage;
using itk::fits::applyFlipImageFilter;
using itk::fits::applyBinomialBlurFilter;
using itk::fits::reflectPixels;
using itk::fits::setNullValue;
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
usage1(bool exitWithFailure, bool verboseUsage)
{
//  if (daProgramName().size()) cerr << basename(daProgramName().c_str());
//  else cerr << "fits2itk";
  
  ostream& out = exitWithFailure ? cerr : cout;
  out << "fits2itk " << fits2itkVersion << "\n\n";
  out << "usage:\n";

  const char* shortUsageFilepath = pteJoinPath(pathToExecutableDir(),
					       "shortUsageMessage.txt");

  auto_ptr<ifstream> shortUsageMessage =
    da::openFileForReading(shortUsageFilepath, da::dieOnError);
  da::copyStream(*shortUsageMessage, out, da::dieOnError);

  if (verboseUsage) {
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
  if (exitWithFailure) exit(EXIT_FAILURE);
  else exit(0);
}


local proc void
usage(bool exitWithFailure=true)
{
  usage1(exitWithFailure, false);
}


local proc void
verboseUsage()
{
  usage1(false, true);
}


//-----------------------------------------------------------------------------
// writeSlicerXmlModuleDescription()
//-----------------------------------------------------------------------------

local proc void
writeSlicerXmlModuleDescription(ostream& out)
{
  const char* moduleDescriptionFilepath =
    pteJoinPath(pathToExecutableDir(),
		"slicerModuleDescription.xml");
  auto_ptr<ifstream> moduleDescription =
    da::openFileForReading(moduleDescriptionFilepath, da::dieOnError);
  da::copyStream(*moduleDescription, out, da::dieOnError);
}


//=============================================================================
// CommandLineParser: local class
//=============================================================================

local proc void
checkEndptr(const char* endptr)
{
  if (*endptr != '\0') usage();
}


// BEGIN local namespace
namespace {

  class CommandLineParser {

    // TODO: Instead of calling usage() within this class upon a parsing error,
    // we should turn on a usage attribute that main() can test for and then
    // call usage().

    // Instance variables:
    const char*    _inputFilepath;
    const char*    _outputFilepath;

    bool           _coerceToShorts;
    bool           _coerceToUnsignedShorts;
    int		   _debugLevel;
    bool	   _dontWrite;
    bool	   _flipDecFlag;
    bool	   _flipRaFlag;
    bool	   _flipVFlag;
    bool	   _logoFlag;
    double	   _nullValue;
    bool	   _reorientNorth;
    // bool	   _transformToEquiangular;
    bool           _slicerXmlModuleDescriptionFlag;


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
    bool        getFlipRaFlag() const { return _flipRaFlag; }
    bool	getFlipVFlag() const { return _flipVFlag; }
    bool        getLogoFlag() const { return _logoFlag; }
    double	getNullValue() const { return _nullValue; }
    bool	getReorientNorth() const { return _reorientNorth; }
//     bool	getTransformToEquiangular() const
//                           { return _transformToEquiangular; }
    bool	getSlicerXmlModuleDescriptionFlag() const
                   { return _slicerXmlModuleDescriptionFlag; }


    bool	getBinomialBlurFlag() const
                   { return _binomialBlurFlag; }
    bool        getDerivativeImageFilterFlag() const
                   { return _derivativeImageFilterFlag; }
    bool        getFlipImageFilterFlag() const
                   { return _flipImageFilterFlag; }
    bool	getIdentityFlipFlag() const
                   { return _identityFlipFlag; }
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
    _flipRaFlag(false),
    _flipVFlag(false),
    _logoFlag(false),
    _nullValue(0.0),
    _reorientNorth(false),
    _slicerXmlModuleDescriptionFlag(false),
    // _transformToEquiangular(false),
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

  int verboseHelpFlag =                false;
  // int equiangularFlag =             false;
  int flipDecFlag =                    false;
  int flipRaFlag =                     false;
  int flipVFlag =                      false;
  int logoFlag  =                      false;
  int noWcsFlag =                      false;
  int reorientNorthFlag =              false;
  int ripOrientationFlag =             false;
  int rotateSkyFlag =                  false;
  int scaleDecFlag =                   false;
  int typicalFlag =                    false;
  int verboseFlag =                    false;
  int slicerXmlModuleDescriptionFlag = false;

  // Specify the allowed short options:
  const char shortopts[] = "a:D:fhnN:o:Rr:Ss:Uv:";

  // Specify the allowed long options:
  struct option longopts[] = {
    // { "equiangular", no_argument,       &equiangularFlag,   true },
    { "help",           no_argument,       &verboseHelpFlag,   true },
    { "flip-dec",       no_argument,       &flipDecFlag,       true },
    { "flip-ra",        no_argument,       &flipRaFlag,        true },
    { "flip-v",         no_argument,       &flipVFlag,         true },
    { "logo",           no_argument,       &logoFlag,          true },
    { "no-wcs",         no_argument,       &noWcsFlag,         true},
    { "nw",             no_argument,       &noWcsFlag,         true},
    { "null-value",     required_argument, 0,                  'N'},
    { "reorient-north", no_argument,       &reorientNorthFlag, true },
    { "rotate-sky",     required_argument, &rotateSkyFlag,     true },
    { "RIP",            no_argument,       &ripOrientationFlag, true },
    { "scale-dec",      required_argument, &scaleDecFlag,       true },
    { "typical",        no_argument,       &typicalFlag,        true },
    { "verbose",        no_argument,       &verboseFlag,        true },
    { "xml",            no_argument,       &slicerXmlModuleDescriptionFlag,
         true },
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

	// FITSImageIO::SetScaleAllAxes(strtod(optarg, &endptr));
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
	// FITSImageIO::SetScaleRA(strtod(optarg, &endptr));
	checkEndptr(endptr);
	break;

      case 'S':
        _coerceToShorts = true;
        break;

      case 's':
        // FITSImageIO::SetScaleVoxelValues(strtod(optarg, &endptr));
	checkEndptr(endptr);
        break;

      case 'U':
        _coerceToUnsignedShorts = true;
        break;

      case 'v':
        // FITSImageIO::SetScaleVelocity(strtod(optarg, &endptr));
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
// 	} else if (equiangularFlag) {
// 	  equiangularFlag =  false;
// 	  _transformToEquiangular = true;
	} else if (flipDecFlag) {
	  flipDecFlag = false;
	  _flipDecFlag = true;
	} else if (flipRaFlag) {
	  flipRaFlag = false;
	  _flipRaFlag = true;
	} else if (flipVFlag) {
	  flipVFlag = false;
	  _flipVFlag = true;
	} else if (logoFlag) {
	  logoFlag = false;
	  _logoFlag = true;
	} else if (noWcsFlag) {
	  noWcsFlag = false;
	  // FITSImageIO::SetSuppressWCS(true);
	} else if (reorientNorthFlag) {
	  reorientNorthFlag = false;
	  _reorientNorth = true;
	} else if (ripOrientationFlag) {
	  ripOrientationFlag = false;
	  // FITSImageIO::SetRIPOrientation(true);
	  // FITSImageIO::SetSuppressWCS(true);
	} else if (rotateSkyFlag) {
	  rotateSkyFlag = false;
	  // FITSImageIO::SetRotateSky(strtod(optarg, &endptr));
	  checkEndptr(endptr);
	} else if (scaleDecFlag) {
	  scaleDecFlag = false;
	  // FITSImageIO::SetScaleDec(strtod(optarg, &endptr));
	  checkEndptr(endptr);
	} else if (typicalFlag) {
	  typicalFlag = false;
	  // FITSImageIO::SetAutoScaleVelocityAxis(true);
	  // FITSImageIO::SetScaleAllAxes(1000);
	  // FITSImageIO::SetScaleRA(-1);
	  // FITSImageIO::SetScaleVoxelValues(1000);
	} else if (verboseFlag) {
	  verboseFlag = false;
	  da::setVerbosityLevel(1);
	} else if (slicerXmlModuleDescriptionFlag) {
	  slicerXmlModuleDescriptionFlag = false;
	  _slicerXmlModuleDescriptionFlag = true;
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
  const int nPositionalParams = argc - ::optind;
  if (_slicerXmlModuleDescriptionFlag or _logoFlag) {
    if (nPositionalParams != 0) usage();
  } else if (_dontWrite) {
    if (nPositionalParams != 1) usage();
  } else {
    if (nPositionalParams != 2) usage();
    _outputFilepath = argv[::optind + 1];
  }
  _inputFilepath = argv[::optind];

  // Do some sanity checking to make sure that the options specified are
  // consistent with each other:
//   if (FITSImageIO::GetSuppressWCS() and
//       FITSImageIO::GetAutoScaleVelocityAxis())
//     {
//       da::warning("Velocity axis auto-scaling does not work when WCS\n"
// 		  "     is suppressed.");
//     }
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
    // FITSImageIO::SetSuppressMetaDataDictionary(true);
  } else if (optionStr == "identityFlip") {
    _identityFlipFlag = true;
  } else usage();
}


} // END local namespace


//=============================================================================
// Local functions
//=============================================================================


//-----------------------------------------------------------------------------
// convertInputFileToItkFile(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc int
convertInputFileToItkFile(CommandLineParser& cl)
{
  printClock("Begin convertInputFileToItkFile()"); //d
  typedef itk::Image<PixelType, c_dims> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  printClock("Before ReaderType::New()"); //d  
  typename ReaderType::Pointer reader = ReaderType::New();
  printClock("Before reader->SetFileName()"); //d
  reader->SetFileName(cl.getInputFilepath());
  printClock("After reader->SetFileName()"); //d
  typename ImageType::Pointer image = reader->GetOutput();
  reflectPixels<PixelType>(*image,
			   cl.getFlipRaFlag(),
			   cl.getFlipDecFlag(),
			   cl.getFlipVFlag());
  printClock("After refect"); //d

  // TODO: Figure out how to do this without reading in the entire image.
  printClock("Before reader->Update()"); //d
  reader->Update();
  printClock("After reader->Update()"); //d
  typename FITSImage<ImageType>::Params params;
  params.itkImage = image;
  FITSImage<ImageType> fitsImage (params);

//   if (cl.getReorientNorth() or cl.getTransformToEquiangular()) {
//     // TODO: Figure out how to do this without reading in the entire image.
//     reader->Update();
//     initializeChangeOfBasis(*image);

//     if (cl.getTransformToEquiangular()) {
//       if (cl.getReorientNorth()) {
// 	transformToNorthOrientedEquiangular(*image);
//       } else {
// 	transformToUnreorientedEquiangular(*image);
//       }
//     } else if (cl.getReorientNorth()) {
//       reorientNorth(*image);
//     }
//   }

  printClock("Before filters"); //d
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
  printClock("After filters, before reflect"); //d

  if (cl.getOutputFilepath()) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(image);
    writer->SetFileName(cl.getOutputFilepath());

    printClock("Before writer->Update()"); //d
    writer->Update();
    printClock("After writer->Update()"); //d
  } else {
    reader->Update();
  }

  if (da::getVerbosityLevel()) writeImageInfo(fitsImage, cout);

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
    if (oldAutoloadPath) { //d
      cout << "oldAutoloadPath=" << oldAutoloadPath << endl; //d
    } else { //d
      cout << "No oldAutoloadPath." << endl; //d
    }
    const string newAutoloadPath =
      oldAutoloadPath
      ? string(pathToExecutableDir()) + pathSeparator + oldAutoloadPath
      : string(pathToExecutableDir());
//    const string newAutoloadPath = string(pathToExecutableDir()); //d
    cout << "newAutoloadPath=" << newAutoloadPath << endl; //d
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

  if (!cl.getSlicerXmlModuleDescriptionFlag()) { //d
    printClock("Main start"); //d
  } //d

  if (cl.getSlicerXmlModuleDescriptionFlag()) {
    writeSlicerXmlModuleDescription(cout);
    return 0;
  }

  if (cl.getLogoFlag()) usage();

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
    printClock("Before conversion"); //d
    status = convertInputFileToItkFile<float>(cl);
    printClock("After conversion"); //d
  }
  cout << "num of factories = " //d
       << itk::ObjectFactoryBase::GetRegisteredFactories().size() //d
       << endl; //d
  printClock("Main end"); //d
  return status;
}
