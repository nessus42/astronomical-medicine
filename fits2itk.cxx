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
    static bool           _cv_coerceToShorts;
    static bool           _cv_coerceToUnsignedShorts;
    static int		  _cv_debugLevel;
    static bool		  _cv_dontWrite;
    static bool		  _cv_flipDecFlag;
    static bool		  _cv_flipRAFlag;
    static bool		  _cv_flipVFlag;

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
  ::opterr = true;
  char* endptr;
  
  char optionChar;
  const int cBase10 = 10;
  string extendedOption;
  int verboseHelpFlag = false;
  int flipDecFlag = false;
  int flipRAFlag = false;
  int flipVFlag = false;
  int noWcsFlag = false;
  int ripOrientationFlag = false;
  int rotateSkyFlag = false;
  int scaleDecFlag = false;
  int typicalFlag = false;
  int verboseFlag = false;

  // Specify the allowed options:
  const char shortopts[] = "Aa:D:fhnN:o:Rr:Ss:Uv:";

  struct option longopts[] = {
    { "help", no_argument, &verboseHelpFlag, true },
    { "flip-dec", no_argument, &flipDecFlag, true },
    { "flip-ra", no_argument, &flipRAFlag, true },
    { "flip-v", no_argument, &flipVFlag, true },
    { "no-wcs", no_argument, &noWcsFlag, true},
    { "rotate-sky", required_argument, &rotateSkyFlag,
      true },
    { "RIP", no_argument, &ripOrientationFlag, true },
    { "scale-dec", required_argument, &scaleDecFlag, true },
    { "typical", no_argument, &typicalFlag, true },
    { "verbose", no_argument, &verboseFlag, true },
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

      case 'n':
	_cv_dontWrite = true;
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
	} else if (verboseFlag) {
	  verboseFlag = false;
	  itk::FITSImageIO::SetVerbose(true);
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
  if (_cv_dontWrite) {
    if (argc - ::optind != 1) ::usage();
  } else {
    if (argc - ::optind != 2) ::usage();
    _cv_outputFilepath = argv[::optind + 1];
  }
  _cv_inputFilepath = argv[::optind];

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
  } else if (optionStr == "flipImageFilter") {
    _cv_flipImageFilterFlag = true;
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
// applyFlipImageFilter(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc typename Image<PixelType, c_dims>::Pointer
applyFlipImageFilter(const typename Image<PixelType, c_dims>::Pointer& image)
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
// isOdd(): local inline function
//-----------------------------------------------------------------------------

local proc inline bool
isOdd(size_t num)
{
  return num & 1;
}


//-----------------------------------------------------------------------------
// degreesToRadians(): local inline function
//-----------------------------------------------------------------------------

local proc inline double
degreesToRadians(double degrees)
{
  return degrees * M_PI/180;
}


//-----------------------------------------------------------------------------
// radiansToDegrees(): local inline function
//-----------------------------------------------------------------------------

local proc inline double
radiansToDegrees(double radians)
{
  return radians * 180/M_PI;
}


//-----------------------------------------------------------------------------
// reflectPixels(): local template function
//-----------------------------------------------------------------------------

// This functions flips the image around the specified axes.  It does this by
// actually moving the pixels around, rather than by messing with the direction
// cosines.  It is, unfortunately, quite hairy.

// CAVEAT: This function assumes that RA is aligned with the i-axis, that Dec
// is aligned with the j-axis, and that V is aligned with the k-axis.  This is
// not always the case, however.  Sometimes a FITS image is rotated in RA/Dec
// space.  Furthermore, even then, this is only an approximation that works for
// small areas of the sky that are not near a pole.  As this function is only
// designed to be used with OsiriX, we can live with this caveat.

template <class PixelType>
local proc void
reflectPixels(Image<PixelType, c_dims>& image,
	      bool flipRAFlag, bool flipDecFlag, bool flipVFlag)
{
  // CAVEAT: This code will break if c_dims ever changes from 3.

  // Return immediately if a NOP is requested:
  if (!flipRAFlag and !flipDecFlag and !flipVFlag) return;

  // Make sure that we actually have the pixels loaded into the image:
  image.Update();

  typedef Image<PixelType, c_dims> ImageType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType SizeType;

  typename ImageType::RegionType allOfImage = image.GetLargestPossibleRegion();
  typename ImageType::SizeType imageSize = allOfImage.GetSize();
  IndexType imageOrigin = allOfImage.GetIndex();

  const size_t c_i    = itk::FITSImageIO::c_i;
  const size_t c_j    = itk::FITSImageIO::c_j;
  const size_t c_k    = itk::FITSImageIO::c_k;

  const size_t minRaIndex  = imageOrigin[c_i];
  const size_t minDecIndex = imageOrigin[c_j];
  const size_t minVelIndex = imageOrigin[c_k];

  const size_t maxRaIndex  = minRaIndex  + imageSize[c_i] - 1;
  const size_t maxDecIndex = minDecIndex + imageSize[c_j] - 1;
  const size_t maxVelIndex = minVelIndex + imageSize[c_k] - 1;

  const size_t middleRaIndex  = minRaIndex  + ((imageSize[c_i] + 1) / 2) - 1;
  const size_t middleDecIndex = minDecIndex + ((imageSize[c_j] + 1) / 2) - 1;
  const size_t middleVelIndex = minVelIndex + ((imageSize[c_k] + 1) / 2) - 1;
  
  const bool lastFlipIsV = flipVFlag;
  const bool lastFlipIsDec = !lastFlipIsV and flipDecFlag;
  const bool lastFlipIsRA = !lastFlipIsV and !lastFlipIsDec and flipRAFlag;

  const size_t raStopIndex = lastFlipIsRA   ? middleRaIndex
                                            : maxRaIndex;
  const size_t decStopIndex = lastFlipIsDec ? middleDecIndex
                                            : maxDecIndex;
  const size_t velStopIndex = lastFlipIsV   ? middleVelIndex
                                            : maxVelIndex;

  const bool needToWorryAboutMiddleVel = flipVFlag and isOdd(imageSize[c_k]);
  const bool needToWorryAboutMiddleDec = flipDecFlag and isOdd(imageSize[c_j]);

  for (size_t vel_i = minVelIndex, velReverse_i = maxVelIndex;
       vel_i <= velStopIndex;
       ++vel_i, --velReverse_i)
    {
      for (size_t dec_i= minDecIndex, decReverse_i = maxDecIndex;
	   dec_i <= decStopIndex;
	   ++dec_i, --decReverse_i)
	{
	  for (size_t ra_i = minRaIndex, raReverse_i = maxRaIndex;
	       ra_i <= raStopIndex;
	       ++ra_i, --raReverse_i)
	    {
	      // Break out of the loop in situations in which we'd be swapping
	      // a pair of pixels back to their original locations due to
	      // swapping the pixels more than once:
	      if (needToWorryAboutMiddleVel) {
		if (vel_i == middleVelIndex) {
		  if (flipDecFlag) {
		    if (dec_i > middleDecIndex) break;
		    else if (needToWorryAboutMiddleDec and
			     flipRAFlag and
			     dec_i == middleDecIndex and
			     ra_i > middleRaIndex) {
		      break;
		    }
		  } else if (flipRAFlag and ra_i > middleRaIndex) break;
		}
	      }

	      IndexType thisPixelIndex;
	      IndexType oppositePixelIndex;
	      thisPixelIndex[c_i] = ra_i;
	      thisPixelIndex[c_j] = dec_i;
	      thisPixelIndex[c_k] = vel_i;
	      oppositePixelIndex[c_i] = flipRAFlag  ? raReverse_i  : ra_i;
	      oppositePixelIndex[c_j] = flipDecFlag ? decReverse_i : dec_i;
	      oppositePixelIndex[c_k] = flipVFlag   ? velReverse_i : vel_i;
	      PixelType tmp = image.GetPixel(thisPixelIndex);
	      image.SetPixel(thisPixelIndex,
			      image.GetPixel(oppositePixelIndex));
	      image.SetPixel(oppositePixelIndex, tmp);
	    }
	}
    }
}


//-----------------------------------------------------------------------------
// writeImageInfo(): local template function
//-----------------------------------------------------------------------------

template <class PixelType>
local proc void
writeImageInfo(const Image<PixelType, c_dims>& image,
	       ostream& out)
{
  typedef Image<PixelType, c_dims> ImageType;
  typename ImageType::RegionType allOfImage = image.GetLargestPossibleRegion();
  typename ImageType::SizeType imageSize = allOfImage.GetSize();
  typename ImageType::IndexType imageOrigin = allOfImage.GetIndex();

  const size_t c_i    = itk::FITSImageIO::c_i;
  const size_t c_j    = itk::FITSImageIO::c_j;
  const size_t c_k    = itk::FITSImageIO::c_k;

  float middleIIndex = imageOrigin[c_i] + imageSize[c_i]/2.0 - 0.5;
  float middleJIndex = imageOrigin[c_j] + imageSize[c_j]/2.0 - 0.5;
  float middleKIndex = imageOrigin[c_k] + imageSize[c_k]/2.0 - 0.5;

  out << "Image center in IJK space with (0,0,0) origin: ("
      << middleIIndex << ", "
      << middleJIndex << ", "
      << middleKIndex << ")\n";

  if (itk::g_theFITSWCSTransform) {
    typedef itk::FITSWCSTransform<double, c_dims> WCS;

    // URGENT TODO: Fix this attrocity!
    const WCS::ConstPointer wcs = (WCS*) itk::g_theFITSWCSTransform;

    typedef WCS::InputPointType IjkPoint;
    IjkPoint ijkCenter;
    ijkCenter[c_i] = middleIIndex;
    ijkCenter[c_j] = middleJIndex;
    ijkCenter[c_k] = middleKIndex;

    typedef WCS::OutputPointType WcsPoint;
    WcsPoint wcsCenter = wcs->TransformPoint(ijkCenter);

    const size_t raIndex  = 0;
    const size_t decIndex = 1;

    const double raAtCenter  = wcsCenter[raIndex];
    const double decAtCenter = wcsCenter[decIndex];


    out << "Image center in RA/Dec coordinates: ("
	<< raAtCenter << ", " << decAtCenter << ")\n";

    // Calculate lengths of unit i and j vectors in RA/Dec space:
    { 
      IjkPoint ijkLeftHalfAPixel  = ijkCenter;
      IjkPoint ijkRightHalfAPixel = ijkCenter;
      IjkPoint ijkDownHalfAPixel  = ijkCenter;
      IjkPoint ijkUpHalfAPixel    = ijkCenter;
      ijkLeftHalfAPixel[c_i]  -= .5;
      ijkRightHalfAPixel[c_i] += .5;
      ijkDownHalfAPixel[c_j]  -= .5;
      ijkUpHalfAPixel[c_j]    += .5;
      
      WcsPoint wcsLeftHalfAPixel  = wcs->TransformPoint(ijkLeftHalfAPixel);
      WcsPoint wcsRightHalfAPixel = wcs->TransformPoint(ijkRightHalfAPixel);
      WcsPoint wcsDownHalfAPixel  = wcs->TransformPoint(ijkDownHalfAPixel);
      WcsPoint wcsUpHalfAPixel    = wcs->TransformPoint(ijkUpHalfAPixel);
      
      // TODO: It's confusing that in some situations V is 1 and dec is 2, and
      // in others, dec is 1 and V is 2.  We need a better way to denote this.

      const double raPerI = 
	wcsRightHalfAPixel[raIndex]  - wcsLeftHalfAPixel[raIndex];
      const double decPerI =
	wcsRightHalfAPixel[decIndex] - wcsLeftHalfAPixel[decIndex];
      const double raPerJ =
	wcsUpHalfAPixel[raIndex]     - wcsDownHalfAPixel[raIndex];
      const double decPerJ =
	wcsUpHalfAPixel[decIndex]    - wcsDownHalfAPixel[decIndex];

//       itk::Vector<double> unitJInWcs = wcsUpHalfAPixel - wcsDownHalfAPixel;
//       cout << "unitJInWcs=" << unitJInWcs << "\n";

//       cout << "raPerI=" << raPerI << "\n";
//       cout << "decPerI=" << decPerI << "\n";
//       cout << "raPerJ=" << raPerJ << "\n";
//       cout << "decPerJ=" << decPerJ << "\n";

//       cout << "ijkLeftHalfAPixel=" << ijkLeftHalfAPixel << "\n";
//       cout << "ijkRightHalfAPixel=" << ijkRightHalfAPixel << "\n";
//       cout << "ijkDownHalfAPixel=" << ijkDownHalfAPixel << "\n";
//       cout << "ijkUpHalfAPixel=" << ijkUpHalfAPixel << "\n";

//       cout << "wcsLeftHalfAPixel=" << wcsLeftHalfAPixel << "\n";
//       cout << "wcsRightHalfAPixel=" << wcsRightHalfAPixel << "\n";
//       cout << "wcsDownHalfAPixel=" << wcsDownHalfAPixel << "\n";
//       cout << "wcsUpHalfAPixel=" << wcsUpHalfAPixel << "\n";

      out << "Unit i vector in RA/Dec space: ("
	  << raPerI << ", " << decPerI << ")\n";
      out << "Unit j vector in RA/Dec space: ("
	  << raPerJ << ", " << decPerJ << ")\n";
      out << "Unit i vector in approximate angular space: ("
	  << raPerI * cos(degreesToRadians(decAtCenter)) << ", "
	  << decPerI << ")\n";
      out << "Unit j vector in approximate angular space: ("
	  << raPerJ * cos(degreesToRadians(decAtCenter)) << ", "
	  << decPerJ << ")\n";
      out << "Image rotation in degrees clockwise:  "
	  << radiansToDegrees(atan(-raPerJ/decPerJ)) << "\n";
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

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFilepath);
  typename ImageType::Pointer image = reader->GetOutput();
  if (Cl::getFlipImageFilterFlag()) {
    image = ::applyFlipImageFilter<PixelType>(image);
  }
  if (Cl::getBinomialBlurFlag()) {
    image = ::applyBinomialBlur<PixelType>(image);
  }
  if (Cl::getIdentityFlipFlag()) {
    image = ::applyFlipImageFilter<PixelType>(image);
    image = ::applyFlipImageFilter<PixelType>(image);
  }
  ::reflectPixels<PixelType>(*image,
			     Cl::getFlipRAFlag(),
			     Cl::getFlipDecFlag(),
			     Cl::getFlipVFlag());
  if (outputFilepath) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(image);
    writer->SetFileName(outputFilepath);
    writer->Update();
  } else {
    reader->Update();
  }

  if (itk::FITSImageIO::GetVerbose()) ::writeImageInfo(*image, cout);
  
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
