// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
//   Program:   fits2itk
//   Module:    fits2itk.cxx
//   Package: 	FITS IO
//   Author:    Douglas Alan <doug AT alum.mit.edu>
//              Initiative in Innovative Computing at Harvard University
//
//   Copyright (c) 2006 Douglas Alan
//
//   This software is freely distributable under the open source MIT X11
//   License.
//
//   See
//
//      http://www.opensource.org/licenses/mite-license
//
//   for details.
//
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
#include <unistd.h>                        // For getopt()

#include <cassert>
#include <string>

#include <itkFITSImageIOFactory.h>
#include <itkFITSImageIO.h>
#include <da_sugar.h>

using std::string;
using itk::DataObject;
using itk::Image;

static const char fits2itkVersion[] = "Version 0.3dev.1";

//-----------------------------------------------------------------------------
// usage()
//-----------------------------------------------------------------------------

// Terminates the program with a return status of EXIT_FAILURE after
// outputting a usage message to 'cerr'.

local proc void
usage()
{
  cerr << fits2itkVersion << "\n\n";
  cerr << "usage:\n";

//  if (daProgramName().size()) cerr << basename(daProgramName().c_str());
//  else cerr << "fits2itk";

  cerr << 

"  fits2itk [-ASU] [-a axes-scale] [-D debug-level] [-N null-value]\n"
"           [-r RA-scale] [-s pixel-scale] [-v velocity-scale]\n"
"           input-file output-file\n"
"\n"
"  A: auto-scale velocity axis\n"
"  S: coerce pixel values to shorts\n"
"  U: coerce pixel values to unsigned shorts\n"
"\n"
"-a scales all the axes, while -v scales only the velocity axis.  If both\n"
"are used, then both scaling factors will be applied.\n"
"\n"
"fits2itk supports CFITSIO's \"extended filename syntax\", which allows\n"
"all sorts of interesting things.  For example, if you have a data cube\n"
"with an extra 1-voxel-thick fourth dimension, you can slice off the\n"
"extra dimension like so:\n"
"\n"
"   fits2itk \"input.fits[*,*,*,1:1][col #NAXIS=3]\" output.nrrd\n"
"\n"
"You can also specify a URL, rather than a filename, for the input file\n"
"and the input file will be automatically fetched via HTTP.\n"
"For the complete manual on the extended filename syntax, see\n"
"the CFITSIO User's Reference Guide chapter on it here:\n"
"\n"
"   http://heasarc.nasa.gov/docs/software/fitsio/c/c_user/node79.html\n"

    ;
  exit(EXIT_FAILURE);
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
  // Parse command line options:
  ::opterr = true;
  char optionChar;
  const char options[] = "Aa:D:fN:o:Rr:Ss:Uv:";
  char** const null = 0;
  const int cBase10 = 10;
  string extendedOption;

  while ((optionChar = ::getopt(argc, const_cast<char**>(argv), options))
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

      case '?': usage();

      default:
        cerr << "Bug!" << endl;
        assert(false);
      } // switch
    } // while

  // Parse the command line positional arguments:
  if (argc - ::optind != 2) usage();
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
  } else usage();
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

// template <class PixelType>
// local proc typename DataObject::Pointer
// flipImage(const typename DataObject::Pointer image)
// {
//   typedef itk::Image<PixelType, 3> ImageType;
//   typedef itk::FlipImageFilter<ImageType> FilterType;
//   typedef typename FilterType::FlipAxesArrayType FlipAxesArrayType;
//   typename FilterType::Pointer filter = FilterType::New();
//   FlipAxesArrayType flipArray;
//   flipArray[0] = 1;
//   flipArray[1] = 1;
//   flipArray[2] = 0;
//   filter->SetFlipAxes(flipArray);
//   filter->SetInput(image);
//   filter->UpdateLargestPossibleRegion();
//   return filter->GetOutput();
// }


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
  // typename DataObject::Pointer image = reader->GetOutput();
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


// TODO: The following version uses the ITK recommended way of doing
// pipelining, while the above code does everything in RAM.  We might want to
// switch back to the "right" way, which is how things are done below.
// Unfortunately, the code below isn't working quite right, and also it is
// requiring seemingly unecessary calls to update().

// //-----------------------------------------------------------------------------
// // convertInputFileToItkFile(): local template function
// //-----------------------------------------------------------------------------

// template <class PixelType>
// local proc int
// convertInputFileToItkFile(const char* const inputFilepath,
//                          const char* const outputFilepath)
// {
//   const unsigned int Dimensions = 3;
//   typedef Image<PixelType, Dimensions> ImageType;
//   typedef itk::ImageFileReader<ImageType> ReaderType;
//   typedef itk::ImageFileWriter<ImageType> WriterType;

//   typename ReaderType::Pointer reader = ReaderType::New();
//   typename WriterType::Pointer writer = WriterType::New();
//   reader->SetFileName(inputFilepath);

//   // This option to run the image through a derivative filter, exists purely
//   // for debugging purposes at the moment:
//   if (::derivativeImageFilterFlag) {
//     cerr << "YO YO!" << endl;
// //     typedef itk::DerivativeImageFilter<ImageType, ImageType> FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     filter->SetOrder(1);
// //     filter->SetDirection(0);

// //     typedef itk::MeanImageFilter<ImageType, ImageType> FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     typename ImageType::SizeType indexRadius;
// //     indexRadius[0] = 5;
// //     indexRadius[1] = 5;
// //     indexRadius[2] = 5;
// //     filter->SetRadius(indexRadius);

// //     typedef itk::BinaryMedianImageFilter<ImageType, ImageType> FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     typename ImageType::SizeType indexRadius;
// //     indexRadius[0] = 1;
// //     indexRadius[1] = 1;
// //     indexRadius[2] = 1;
// //     filter->SetRadius(indexRadius);

//     typedef itk::FlipImageFilter<ImageType> FilterType;
//     typedef typename FilterType::FlipAxesArrayType FlipAxesArrayType;
//     typename FilterType::Pointer filter = FilterType::New();
//     FlipAxesArrayType flipArray;
//     flipArray[0] = 1;
//     flipArray[1] = 1;
//     flipArray[2] = 1;
//     filter->SetFlipAxes(flipArray);

// //     typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType>
// //       FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     filter->SetNumberOfIterations(1);
// //     filter->SetTimeStep(0.125);
// //     filter->SetConductanceParameter(10);

//     filter->SetInput(reader->GetOutput());
//     filter->Update();
//     writer->SetInput(filter->GetOutput());

//   } else {
//     writer->SetInput(reader->GetOutput());
//   }
//   writer->SetFileName(outputFilepath);
//   writer->Update();
//   return EXIT_SUCCESS;

//   // Note: The above call to Update() might raise an execption.  If you were to
//   // want to catch this exception, here is now you might do it.  You don't
//   // really need to, though, as the default exception handler will output a
//   // message that is not particularly more cryptic than the following:
//   //
//   //   try {
//   //     writer->Update(); 
//   //   } catch(itk::ExceptionObject& err) { 
//   //     cerr << "ExceptionObject caught !" << endl; 
//   //     cerr << err << endl; 
//   //     return EXIT_FAILURE;
//   //   } 
// }

} // namespace

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


//=============================================================================
// Change Log
//=============================================================================

// See http://wiki.iic.harvard.edu/AstroMed/Software_Version_Numbers for
// documentation on how these version numbers work.


//---------------------------
// Version 0.1dev.0
//---------------------------

// *** Tue Dec  5, 2006 ***

// Added version number to usage output and a change log comments section.

// Changed autoscaling of velocity axis to be average of raPerI and decPerI,
// rather than just raPerI.

// RA/DEC/V were getting mapped to L/P/S rather than to L/P/S.  Fixed this.

// Added option to turn on debugging messages and turned them off by default.
// They were all getting output all the time prior to this change.

// Note this version not checked into subversion due to repository being down.

//---------------------------
// Version 0.1dev.1
//---------------------------

// *** Thu Dec  7, 2006 ***

// Added option to scale all axes and removed hack that automatically scaled
// all axes by 1000.

// Implemented a better method of autoscaling the velocity axis by
// fitting it maximally into cube.

// Completely reworked the change-of-basis matrix stuff to be clearer
// and to actually work correctly.

// *** Fri Dec  8, 2006 ***

// Added option to scale RA axis.  The most important purpose for this option
// is the scale RA by -1 so that it will display properly in Slicer.

//---------------------------
// Version 0.1
//---------------------------

// *** Tue Dec 12, 2006 ***

// Changed version to "0.1".  Could not check it into the NA-MIC sandbox due to
// the corruption of their repository.  So started my own Mercurial repository.

//---------------------------
// Version 0.2dev.0
//---------------------------

// *** Wed Dec 31, 2007 ***

// itkFITSImageIO.cxx: Fixed bug that caused core dump when filename extention
// test was longer than actual filename.

// cfitsio/CMakeLists.txt: Modified to set the correct "-D" options to gcc.
// (At least for OS X.)

// libwcs/fitshead.h: Added a #define for "dec2str" and "str2dec" because
// functions were conflicting with functions of the same name somewhere in
// Carbon.

// Recompiled with ITK 3.0.1 so as to get turn off preemptive denial of reading
// non-existant files.  (This is needed to support cfitsio's extended file
// syntax, which allows you to specify "filenames" that indicate more than just
// the filename.

// *** Thu Feb  1, 2007 ***

// itkFITSImageIO.cxx: Modified checkExtension() to use a complicated regular
// expression, rather than just doing a string compare on the end of the
// string.

//---------------------------
// Version 0.2
//---------------------------

// *** Mon Feb  5, 2007 ***

// Added a command line option to set the FITS null value for reading FITS
// files.

// *** Thu Feb  8, 2007 ***

// Added support for debug level -1, which is for testing out my
// itk::Transform.  The transform has been copied from
// itk::TranslationTransform, but I haven't done anything to it yet.  I've just
// made sure that it compiles.

// *** Tue Feb 13, 2007 ***

// Made itkFITSWCSTransform.h and itkFITSWCSTransform.txx and got them
// to compile.  They don't do anything interesting yet.

// Prettied up all the FITS Reader source files.

// *** Thu Feb 15, 2007 ***

// Added note to usage string telling how to slice off a useless forth
// dimension.

//---------------------------
// Version 0.3dev.0
//---------------------------

// *** Tue Feb 20, 2007 ***

// Added WorldCoor to FITSWCSTransform.

// *** Fri Feb 23, 2007

// Got both forward and reverse transforms to work for FITSWCSTransform, and
// made the debugging output prove this fact.

//---------------------------
// Version 0.3dev.1pending
//---------------------------

// *** Mon Feb 26, 2007 ***

// Modified to optionally feed the image through an ITK FlipImageFilter and see
// if the metadata dictionary is preserved.  Alas, it wasn't.  Also, the
// resulting nrrd output crashes Slicer, and if I try to output an Analyze file
// instead, I get an error about the image not being in RPI, PIR, or RIP
// orientation.

//---------------------------
// Version 0.3dev.1
//---------------------------

// *** Mon Mar  5, 2007 ***

// Added an option to apply the ITK BinomialBlurImageFilter.

// Reworked the code that applies filters so that you can apply multiple
// filters.

// Moved command line parsing into a class.

//----------------------------------------------------------------------
// *** Changes described above this line are checked in to Mercurial ***
//----------------------------------------------------------------------
