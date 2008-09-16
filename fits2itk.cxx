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

#include <libgen.h>                        // For basename()
#include <getopt.h>
#include <stdlib.h>                        // For getenv(), setenv()
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

#include <itkMetaDataObject.h>
using itk::MetaDataDictionary;
using itk::MetaDataObject;
using itk::MetaDataObjectBase;

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVector.h>

#include <itkFITSImageIO.h>
using itk::FITSImageIO;

#include <itkFITSWCSTransform.h>
using itk::FITSWCSTransform;

#include <itkFITSImageUtils.h>
using itk::fits::FITSImage;
using itk::fits::reflectPixels;
using itk::fits::scalePixelValues;
using itk::fits::setNullValue;
using itk::fits::setFITSImageIODebugLevel;
using itk::fits::writeImageInfo;
using itk::fits::writeFitsHeader;

#include <pathToExecutable.h>
#include <da_util.h>

#include <da_sugar.h>
using da::setDebugLevel;
using da::endMatchesP;

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


//=============================================================================
// CommandLineParser: local class
//=============================================================================

// TODO: Move this into its own file.  It should probably go into a "fits2itk"
// namespace.

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

//     bool	   _quietModeP;          // Set with -q, --quiet //?
    bool           _dontWriteP;          // Set with -n, --no-write //?
    bool	   _wcsP;                // Set with --wcs
    bool	   _equiangularP;        // Set with --equiangular
    bool	   _northUpP;            // Set with --north-up
//     bool	   _eastLeftP;           // Set with --east-left
    bool	   _autoscaleZAxisP;	 // Set with --autoscale-z-axis
    bool	   _flipxP;		 // Set with --flipx
    bool	   _flipyP;		 // Set with --flipy
    bool	   _flipzP;		 // Set with --flipz
    bool	   _verboseP;		 // Set with --verbose

    //? I don't think that we need both --quiet and --verbose, do we?

    bool	   _showFitsHeaderP;     // Set with --show-fits-header //?
    bool	   _coerceToShortsP;	 // Set with --coerce-to-shorts
    bool	   _coerceToUnsignedShortsP;
    bool	   _outputFileIsInAnalyzeFormat;
//     bool	   _lpsP;		 // Set with --lps //?
    
    double	   _pixelScale;	         // Set with --pixel-scale
    double         _xAxisScale;		 // Set with --x-scale
    double 	   _yAxisScale;		 // Set with --y-scale
    double	   _zAxisScale; 	 // Set with --z-scale
    double	   _nullValue;		 // Set with --null-value //?
    double	   _debugLevel;		 // Set with --debug-level
    double	   _rotateSky;		 // Set with --rotate-sky




//     unsigned       _wcsGridStride;	 // Set with --wcs-grid-stride //?

  public:
    
    // Constructor:
    CommandLineParser(const int argc, const char* const argv[]);

    // Accessor methods:
    const char* inputFilepath() const { return _inputFilepath;}
    const char* outputFilepath() const { return _outputFilepath; }

//     bool quietModeP() const { return _quietModeP; }
    bool dontWriteP() const { return _dontWriteP; }
    bool wcsP() const { return _wcsP; }
    bool equiangularP() const { return _equiangularP; }
    bool northUpP() const { return _northUpP; }
//     bool eastLeftP() const { return _eastLeftP; }
    bool autoscaleZAxisP() const { return _autoscaleZAxisP; }
    bool flipxP() const { return _flipxP; }
    bool flipyP() const { return _flipyP; }
    bool flipzP() const { return _flipzP; }
    bool verboseP() const { return _verboseP; }
    bool showFitsHeaderP() const { return _showFitsHeaderP; }
    bool coerceToShortsP() const { return _coerceToShortsP; }
    bool coerceToUnsignedShortsP() const { return _coerceToUnsignedShortsP; }
    bool outputFileIsInAnalyzeFormat() const
            { return _outputFileIsInAnalyzeFormat; }
//     bool lpsP() const { return _lpsP	 }


    double pixelScale() const { return _pixelScale; }
    double xAxisScale() const { return _xAxisScale; }
    double yAxisScale() const { return _yAxisScale; }
    double zAxisScale() const { return _zAxisScale; }
    double nullValue() const { return _nullValue; }
    double debugLevel() const { return _debugLevel; }
    double rotateSky() const { return _rotateSky; }

//     unsigned wcsGridStride() const { return _wcsGridStride; }
  };

ctor
CommandLineParser::CommandLineParser(const int argc, const char* const argv[])
  : _inputFilepath(0),
    _outputFilepath(0),
//     _quietModeP(false),
    _dontWriteP(false),
    _wcsP(false),
    _equiangularP(false),
    _northUpP(false),
//     _eastLeftP(false),
    _autoscaleZAxisP(false),
    _flipxP(false),
    _flipyP(false),
    _flipzP(false),
    _verboseP(false),
    _showFitsHeaderP(false),
    _coerceToShortsP(false),
    _coerceToUnsignedShortsP(false),
    _outputFileIsInAnalyzeFormat(false),
//     _lpsP(false),
    _pixelScale(1),
    _xAxisScale(1),
    _yAxisScale(1),
    _zAxisScale(1),
    _nullValue(0),  // 0 means "leave alone", not set to 0
    _debugLevel(0),
    _rotateSky(0)
//     _wcsGridStride(0)
{
  opterr = true;
  char* endptr;
  
  char optionChar;
  const int c_base10 = 10;

  int verboseHelpParsingFlag = false;
//   int quietModeParsingFlag = false;
  int wcsParsingFlag = false;
  int equiangularParsingFlag = false;
  int northUpParsingFlag = false;
//int eastLeftParsingFlag  = false;
  int autoscaleZAxisParsingFlag = false;
  int flipxParsingFlag = false;
  int flipyParsingFlag = false;
  int flipzParsingFlag  = false;
  int verboseParsingFlag  = false;
  int showFitsHeaderParsingFlag  = false;
  int coerceToShortsParsingFlag  = false;
  int coerceToUnsignedShortsParsingFlag  = false;
//int lpsParsingFlag  = false;
  int pixelScaleParsingFlag  = false;
  int xAxisScaleParsingFlag  = false;
  int yAxisScaleParsingFlag  = false;
  int zAxisScaleParsingFlag  = false;
  int skyScaleParsingFlag  = false;
  int nullValueParsingFlag  = false;
  int debugLevelParsingFlag  = false;
  int rotateSkyParsingFlag  = false;
//int wcsGridStrideParsingFlag  = false;

  // Specify the allowed short options:
  const char shortopts[] = "hnq?";

  // Specify the allowed long options:
  struct option longopts[] = {
    { "help",             no_argument,       &verboseHelpParsingFlag,    true },
//     { "quiet",		  no_argument,       &quietModeParsingFlag,      'q' },
    { "no-write",	  no_argument,       null,      		 'n' },
    { "wcs",              no_argument,       &wcsParsingFlag,            true },
    { "equiangular",      no_argument,       &equiangularParsingFlag,    true },
    { "north-up",         no_argument,       &northUpParsingFlag,        true },
//  { "east-left",        no_argument,       &eastLeftParsingFlag,       true },
    { "autoscale-z-axis", no_argument,       &autoscaleZAxisParsingFlag, true },
    { "flipx",            no_argument,       &flipxParsingFlag,          true },
    { "flipy",            no_argument,       &flipyParsingFlag,          true },
    { "flipz",            no_argument,       &flipzParsingFlag,          true },
    { "verbose",          no_argument,       &verboseParsingFlag,        true },
    { "show-fits-header", no_argument,       &showFitsHeaderParsingFlag, true },
    { "coerce-to-shorts", no_argument,       &coerceToShortsParsingFlag, true },
    { "coerce-to-unsigned-shorts",
                          no_argument,       &coerceToUnsignedShortsParsingFlag,
      									 true },
//  { "lps",              no_argument,       &lpsParsingFlag,            true },
    { "pixel-scale",      required_argument, &pixelScaleParsingFlag,     true },
    { "x-scale",          required_argument, &xAxisScaleParsingFlag,     true },
    { "y-scale",          required_argument, &yAxisScaleParsingFlag,     true },
    { "z-scale",          required_argument, &zAxisScaleParsingFlag,     true },
    { "sky-scale",	  required_argument, &skyScaleParsingFlag, 	 true },
    { "null-value",       required_argument, &nullValueParsingFlag,      true },
    { "debug-level",      required_argument, &debugLevelParsingFlag,     true },
    { "rotate-sky",       required_argument, &rotateSkyParsingFlag,      true },
//  { "wcs-grid-stride",  required_argument, &wcsGridStrideParsingFlag,  true },
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

      case 'h':
        usage();

      case 'n':
        _dontWriteP = true;
        break;

//       case 'q':
//         _quietModeP = true;
//         break;

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
        if (verboseHelpParsingFlag) {
          verboseUsage();
        } else if (wcsParsingFlag) {
          wcsParsingFlag = false;
          _wcsP = true;
        } else if (equiangularParsingFlag) {
          equiangularParsingFlag = false;
          _equiangularP = true;
        } else if (northUpParsingFlag) {
          northUpParsingFlag = false;
          _northUpP = true;
//         } else if (eastLeftParsingFlag) {
//           eastLeftParsingFlag = false;
//           _eastLeftP = true;
        } else if (autoscaleZAxisParsingFlag) {
          autoscaleZAxisParsingFlag = false;
          _autoscaleZAxisP = true;
        } else if (flipxParsingFlag) {
          flipxParsingFlag = false;
          _flipxP = true;
        } else if (flipyParsingFlag) {
          flipyParsingFlag = false;
          _flipyP = true;
        } else if (flipzParsingFlag) {
          flipzParsingFlag = false;
          _flipzP = true;
        } else if (verboseParsingFlag) {
          verboseParsingFlag = false;
          _verboseP = true;
        } else if (showFitsHeaderParsingFlag) {
          showFitsHeaderParsingFlag = false;
          _showFitsHeaderP = true;
        } else if (coerceToShortsParsingFlag) {
          coerceToShortsParsingFlag = false;
          _coerceToShortsP = true;
        } else if (coerceToUnsignedShortsParsingFlag) {
          coerceToUnsignedShortsParsingFlag = false;
          _coerceToUnsignedShortsP = true;
//         } else if (lpsParsingFlag) {
//           lpsParsingFlag = false;
//           _lpsP = true;
        } else if (pixelScaleParsingFlag) {
          pixelScaleParsingFlag = false;
          _pixelScale *= strtod(optarg, &endptr);
          checkEndptr(endptr);
        } else if (xAxisScaleParsingFlag) {
          xAxisScaleParsingFlag = false;
          _xAxisScale *= strtod(optarg, &endptr);
          checkEndptr(endptr);
        } else if (yAxisScaleParsingFlag) {
          yAxisScaleParsingFlag = false;
          _yAxisScale *= strtod(optarg, &endptr);
          checkEndptr(endptr);
        } else if (zAxisScaleParsingFlag) {
          zAxisScaleParsingFlag = false;
          _zAxisScale *= strtod(optarg, &endptr);
          checkEndptr(endptr);
        } else if (skyScaleParsingFlag) {
          skyScaleParsingFlag = false;
          _xAxisScale *= strtod(optarg, &endptr);
          checkEndptr(endptr);
          _yAxisScale *= strtod(optarg, &endptr);
          checkEndptr(endptr);
        } else if (nullValueParsingFlag) {
          nullValueParsingFlag = false;
          _nullValue = strtod(optarg, &endptr);
          checkEndptr(endptr);
          if (_nullValue == 0) runTimeError("Cannot set the null value to 0.");
        } else if (debugLevelParsingFlag) {
          debugLevelParsingFlag = false;
          _debugLevel = strtol(optarg, &endptr, c_base10);
          checkEndptr(endptr);
        } else if (rotateSkyParsingFlag) {
          rotateSkyParsingFlag = false;
          _rotateSky = strtod(optarg, &endptr);
          checkEndptr(endptr);
//         } else if (wcsGridStrideParsingFlag) {
//           wcsGridStrideParsingFlag = false;
//           _debugLevel = strtol(optarg, &endptr, c_base10);
//           checkEndptr(endptr);
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

  if (_dontWriteP) {
    if (nPositionalParams != 1) usage();
  } else {
    if (nPositionalParams != 2) usage();
    _outputFilepath = argv[::optind + 1];
  }
  _inputFilepath = argv[::optind];

  if (da::endMatchesP(_outputFilepath, ".hdr") or
      da::endMatchesP(_outputFilepath, ".img"))
    {
      _outputFileIsInAnalyzeFormat = true;
    } 
    
  // Do some sanity checking to make sure that the options specified are
  // consistent with each other:

  if (_wcsP and (_flipxP or _flipyP or _flipzP)) {
    runTimeError("You cannot use any of the flip options if WCS is turned on.");
  }

  if (_wcsP and _outputFileIsInAnalyzeFormat) {
    runTimeError("You cannot turn on WCS when writing file in Analyze format.");
  }

  if (_autoscaleZAxisP and !_wcsP) {
    runTimeError("Velocity axis auto-scaling does not work without WCS on.\n");
  }
  if (_northUpP and _equiangularP) {
    runTimeError("\"--equiangular causes north to be oriented up, so "
                 "please don't specify both.");
  }
  if (_northUpP and _wcsP) {
    runTimeError("\"--wcs\" causes north to be oriented up, so please don't "
                 "specify both.");
  }
}


} // END local namespace


//=============================================================================
// Local functions
//=============================================================================

//-----------------------------------------------------------------------------
// convertInputFileToItkFile(): local template function
//-----------------------------------------------------------------------------

template <class PixelT>
local proc int
convertInputFileToItkFile(CommandLineParser& cl)
{
  typedef itk::Image<PixelT, c_dims> ImageT;
  typedef itk::ImageFileReader<ImageT> ReaderT;

  typename ReaderT::Pointer reader = ReaderT::New();
  reader->SetFileName(cl.inputFilepath());
  typename ImageT::Pointer image = reader->GetOutput();
  reflectPixels(*image,
                cl.flipxP(),
                cl.flipyP(),
                cl.flipzP());

  // Note: scalePixelValues() uses a pixel iterator to scale all the pixels.
  // I'm guessing that iterating the pixels can't be done in a streaming
  // manner.  If we ever get to the point where we can do streaming i/o, then
  // I'd probably want to implement the pixel value scaling using a unary
  // functor filter instead:

  scalePixelValues(*image, cl.pixelScale());

  // The call to Update() immediately below causes the image to be read in to
  // ram.  This is not ideal as it may cause more ram to be used that would
  // otherwise be required.  I think that the only thing requiring that
  // Update() be called here is so that the MDD will exist when we get to the
  // construction of `fitsImage` immediately following.  TODO: Figure out how
  // how to avoid having to call Update() here:

  reader->Update();

  typename FITSImage<ImageT>::Params params;
  params.itkImage = image;
  params.wcsP = cl.wcsP();
  params.equiangularP = cl.equiangularP();
  params.northUpP = cl.northUpP();
//   params.eastLeftP = cl.eastLeftP();
  params.autoscaleZAxisP = cl.autoscaleZAxisP();
//   params.lpsP = cl.lpsP();
  params.xAxisScale = cl.xAxisScale();
  params.yAxisScale = cl.yAxisScale();
  params.zAxisScale = cl.zAxisScale();
  params.rotateSky = cl.rotateSky();

  if (cl.outputFileIsInAnalyzeFormat()) params.zAxisScale *= -1;

  FITSImage<ImageT> fitsImage (params);

  // We actually only write an output file if we were given the name of an
  // output file.  If we weren't, then we just read in the image for the
  // purposes of outputing any information to stdout that has been requested:

  if (cl.outputFilepath()) {
    typedef itk::ImageFileWriter<ImageT> WriterT;
    typename WriterT::Pointer writer = WriterT::New();
    writer->SetInput(image);
    writer->SetFileName(cl.outputFilepath());
    writer->Update();
  } else {
    reader->Update();
  }

  if (da::getVerbosityLevel()) writeImageInfo<ImageT>(fitsImage, cout);
  if (cl.showFitsHeaderP()) writeFitsHeader<ImageT>(*image, cout);

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
// mapRange(): local inline function
//-----------------------------------------------------------------------------

// Implements this function:
//
// http://processing.org/reference/map_.html

local inline double
mapRange(const double value,
         const double sourceLow, double sourceHigh,
         const double destLow, double destHigh)
{
  return (value - sourceLow) / (sourceHigh - sourceLow) * (destHigh - destLow)
    + destLow;
}


local inline double
mapRange(const int value,
         const int sourceLow, int sourceHigh,
         const double destLow, double destHigh)
{
  if (value == sourceLow) return destLow;
  if (value == sourceHigh) return destHigh;
  return mapRange(
            static_cast<double>(value),
            static_cast<double>(sourceLow), static_cast<double>(sourceHigh),
            destLow, destHigh);
}


//-----------------------------------------------------------------------------
// convertInputFileToWcsGridImage(): local function
//-----------------------------------------------------------------------------

// TODO: This function probably won't work right in the current version of
// fits2itk, as it doesn't support all the options that we now have, but I'm
// not sure that it should be removed, as we may still want it someday.  I
// think I'll just comment it out for now, but leave it here for future
// reference.

// local proc int
// convertInputFileToWcsGridImage(CommandLineParser& cl)
// {
//   typedef itk::Image<short, c_dims> InputImage;
//   typedef itk::ImageFileReader<InputImage> Reader;

//   Reader::Pointer reader = Reader::New();
//   reader->SetFileName(cl.getInputFilepath());
//   InputImage::Pointer inputImage = reader->GetOutput();

//   // TODO: Hopefully we can get rid of this call to Update() at some point.
//   // See note about this in convertInputFileToItkFile():
//   reader->Update();

//   FITSImage<InputImage>::Params params;
//   params.itkImage = inputImage;
//   FITSImage<InputImage> fitsImage (params);

//   // Create a new image to hold the WCS information:
//   typedef itk::Vector<double, 3> WcsMapPixel;
//   typedef itk::Image<WcsMapPixel, c_dims> WcsMapImage;
//   WcsMapImage::Pointer wcsMapImage = WcsMapImage::New();

//   // Extract the velocity information from the MDD:
//   double velocityAtIndexOrigin;
//   double velocityDelta;
//   {
//     const MetaDataDictionary& mdd = inputImage->GetMetaDataDictionary();
//     const MetaDataObjectBase* const velocityAtIndexOriginMdob
//       = mdd["fits2itk.velocityAtIndexOrigin"];
//     const MetaDataObjectBase* const velocityDeltaMdob
//       = mdd["fits2itk.velocityDelta"];
    
//     // TODO: Add error checking.  If these items aren't in the MDD, we'll dump
//     // core.
    
//     const MetaDataObject<double>* const velocityAtIndexOriginMdo
//       = dynamic_cast<const MetaDataObject<double>*>(velocityAtIndexOriginMdob);
//     const MetaDataObject<double>* const velocityDeltaMdo
//       = dynamic_cast<const MetaDataObject<double>*>(velocityDeltaMdob);

//     velocityAtIndexOrigin = velocityAtIndexOriginMdo->GetMetaDataObjectValue();
//     velocityDelta = velocityDeltaMdo->GetMetaDataObjectValue();
    
//     debugPrint("velocityAtIndexOrigin=" << velocityAtIndexOrigin);
//     debugPrint("velocityDelta=" << velocityDelta);
//   }


//   // Write the WCS info into WCS map image:
//   {
//     typedef WcsMapImage::IndexType WcsMapIndex;

//     enum { c_i = FITSImageIO::c_i,
//            c_j = FITSImageIO::c_j,
//            c_k = FITSImageIO::c_k };

//     InputImage::RegionType allOfInputImage =
//       inputImage->GetLargestPossibleRegion();
//     InputImage::SizeType inputImageSize = allOfInputImage.GetSize();
//     InputImage::IndexType inputImageOrigin = allOfInputImage.GetIndex();
    
//     const size_t minRaIndex   = inputImageOrigin[c_i];
//     const size_t minDecIndex = inputImageOrigin[c_j];
//     const size_t minVelIndex  = inputImageOrigin[c_k];

//     const size_t maxRaIndex  = minRaIndex  + inputImageSize[c_i] - 1;
//     const size_t maxDecIndex = minDecIndex + inputImageSize[c_j] - 1;
//     const size_t maxVelIndex = minVelIndex + inputImageSize[c_k] - 1;

//     debugPrint("wcsImageStride="
//                << cl.getWcsImageStride());

//     const size_t numOfRaStrides =
//       size_t(ceil((maxRaIndex - minRaIndex) / 
//                   cl.getWcsImageStride()) + 1);
//     const size_t numOfDecStrides =
//       size_t(ceil((maxDecIndex - minDecIndex) /
//                   cl.getWcsImageStride()) + 1);
//     const size_t velIndexCount = 2;

//     FITSImage<InputImage>::WcsTransformConstPtr wcsTransform =
//       fitsImage.wcsTransform();
//     WcsMapImage::IndexType pixelIndex;
//     WcsMapImage::PointType pixelPoint;
//     WcsMapImage::PixelType wcsPixel;
//     WcsMapImage::PointType wcsPoint;

//     // Size the WCS map image:
//     {
//       WcsMapImage::RegionType wcsMapDims;
//       wcsMapDims.SetSize(c_i, numOfRaStrides);
//       wcsMapDims.SetSize(c_j, numOfDecStrides);
//       wcsMapDims.SetSize(c_k, velIndexCount);
//       wcsMapImage->SetRegions(wcsMapDims);
//       wcsMapImage->Allocate();
//     }

//     for (unsigned slice = 0; slice < velIndexCount; ++slice) {
//       pixelIndex[c_k] = slice;
//       size_t vel_i = (slice == 0) ? minVelIndex : maxVelIndex;
//       double velocity = velocityAtIndexOrigin + vel_i * velocityDelta;
//       // pixelPoint[c_k] = vel_i;

//       pixelIndex[c_j] = 0;
//       for (size_t decStrideNum = 1;
//            decStrideNum <= numOfDecStrides;
//            ++decStrideNum, ++pixelIndex[c_j])
//         {
//           pixelPoint[c_j] = mapRange(decStrideNum,
//                                      1, numOfDecStrides,
//                                      minDecIndex, maxDecIndex);
//           pixelIndex[c_i] = 0;
//           for (size_t raStrideNum = 1;
//                raStrideNum <= numOfRaStrides;
//                ++raStrideNum, ++pixelIndex[c_i])
//             {
//               pixelPoint[c_i] = mapRange(raStrideNum, 
//                                          1, numOfRaStrides,
//                                          minRaIndex, maxRaIndex);
//               wcsPoint = wcsTransform->TransformPoint(pixelPoint);
//               wcsPixel[c_i] = wcsPoint[c_i];
//               wcsPixel[c_j] = wcsPoint[c_j];
//               wcsPixel[c_k] = velocity;
//               wcsMapImage->SetPixel(pixelIndex, wcsPixel);
//             }
//         }
//     }
//   }

//   // Write `wcsMapImage` to the output file:
//   typedef itk::ImageFileWriter<WcsMapImage> Writer;
//   Writer::Pointer writer = Writer::New();
//   writer->SetInput(wcsMapImage);
//   writer->SetFileName(cl.getOutputFilepath());
//   writer->Update();

//   return EXIT_SUCCESS;
// }

//-----------------------------------------------------------------------------
// handleOptions(): local function
//-----------------------------------------------------------------------------

local proc void
handleOptions(const CommandLineParser& cl)
{
  setNullValue(cl.nullValue());
  da::setDebugLevel(cl.debugLevel());
  setFITSImageIODebugLevel(cl.debugLevel());
  da::setVerbosityLevel(cl.verboseP());
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
    const string newAutoloadPath =
      oldAutoloadPath
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

// TODO: Make this work again someday, if need be:

//   if (cl.getWcsImageStride() > 0) {
//     status = convertInputFileToWcsGridImage(cl);
//   } else 

  if (cl.coerceToShortsP()) {
    status = convertInputFileToItkFile<short>(cl);
  } else if (cl.coerceToUnsignedShortsP()) {
    status = convertInputFileToItkFile<unsigned short>(cl);
  } else {
    status = convertInputFileToItkFile<float>(cl);
  }
  return status;
}
