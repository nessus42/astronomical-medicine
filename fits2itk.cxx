//===============================================================================
//
//   Program:   fits2itk
//   Module:    fits2itk.cxx
//   Language:  C++
//   Author:    Douglas Alan <doug AT alum.mit.edu>
//
//   Copyright (c) 2006 Douglas Alan
//
//   This software is freely distributable under the open source MIT X11 License.
//   See
//
//      http://www.opensource.org/licenses/mit-license.php
//
//   for details.
//
// ==============================================================================


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <unistd.h>                        // For getopt()
#include <libgen.h>			   // For basename()
#include <assert.h>

#include <itkFITSImageIOFactory.h>
#include <itkFITSImageIO.h>
#include <da_sugar.h>

using std::cerr;
using std::endl;


//-----------------------------------------------------------------------------
// Local error procedures
//-----------------------------------------------------------------------------

// Terminates the program with a return status of EXIT_FAILURE after
// outputting a usage message to 'cerr'.

local proc void
usage()
{
  cerr << "\nVersion 0.2dev.1pending\n\n";
  cerr << "usage: ";
  if (daProgramName().size()) cerr << basename(daProgramName().c_str());
  else cerr << "fits2itk";
  cerr << 
    " [-ASU] [-a axes-scale] [-D debug-level] [-N null-value] [-r RA-scale]\n"
    "    [-s pixel-scale] [-v velocity-scale] input-file output-file\n"
    "\n"
    "  A: auto-scale velocity axis\n"
    "  S: coerce pixel values to shorts\n"
    "  U: coerce pixel values to unsigned shorts\n"
    "\n"
    "-a scales all the axes, while -v scales only the velocity axis.  If both\n"
    "are used, then both scaling factors will be applied.\n"
    ;
  exit(EXIT_FAILURE);
}


//-----------------------------------------------------------------------------
// Local procedures
//-----------------------------------------------------------------------------

template <class PixelType>
local proc int
convertFitsFileToItkFile(const char* const inputFilepath,
                         const char* const outputFilepath)
{
  const unsigned int Dimensions = 3;
  typedef itk::Image<PixelType, Dimensions> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  typename WriterType::Pointer writer = WriterType::New();
  reader->SetFileName(inputFilepath);
  writer->SetFileName(outputFilepath);
  writer->SetInput(reader->GetOutput());

  try {
    writer->Update(); 
  } catch(itk::ExceptionObject& err) { 
    cerr << "ExceptionObject caught !" << endl; 
    cerr << err << endl; 
    return EXIT_FAILURE;
  } 
  return EXIT_SUCCESS;
}


//-----------------------------------------------------------------------------
// main()
//-----------------------------------------------------------------------------

proc int
main(const int argc, const char* const argv[])
{
  
  ::daSetProgramName(argv[0]);

  // Flags controllable from command line:
  bool coerceToShorts = false;
  bool coerceToUnsignedShorts = false;

  // Parse command line options:
  ::opterr = true;
  char optionChar;
  const char options[] = "Aa:D:fN:Rr:Ss:Uv:";
  char** const null = 0;
  const int cBase10 = 10;

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
	itk::FITSImageIO::SetDebugLevel(strtol(optarg, null, cBase10));
	break;

      case 'N':
	itk::FITSImageIO::SetNullValue(strtod(optarg, null));;
	break;

      case 'R':
	// TODO: This isn't implemented yet.
	itk::FITSImageIO::SetRotateSky(strtod(optarg, null));
	break;

      case 'r':
	itk::FITSImageIO::SetScaleRA(strtod(optarg, null));
	break;

      case 'S':
        coerceToShorts = true;
        break;

      case 's':
        itk::FITSImageIO::SetScaleVoxelValues(strtod(optarg, null));
        break;

      case 'U':
        coerceToUnsignedShorts = true;
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
      }
    }


  // Parse the command line positional arguments:
  if (argc - ::optind != 2) usage();
  const char* inputFilepath = argv[::optind];
  const char* outputFilepath = argv[::optind + 1];

  // Register FITS one factory with the ImageIOFactory.
  itk::FITSImageIOFactory::RegisterOneFactory();

  int status = -99;
  if (coerceToShorts) {
    status = convertFitsFileToItkFile<short>(inputFilepath, outputFilepath);
  } else if (coerceToUnsignedShorts) {
    status = convertFitsFileToItkFile<unsigned short>(inputFilepath,
                                                      outputFilepath);
  } else {
    status = convertFitsFileToItkFile<float>(inputFilepath, outputFilepath);
  }
  assert(status != -99);
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

//----------------------------------------------------------------------
// *** Changes described above this line are checked in to Mercurial ***
//----------------------------------------------------------------------

//---------------------------
// Version 0.2dev.1pending
//---------------------------

// *** Mon Feb  5, 2007 ***

// Added a command line option to set the FITS null value for reading FITS
// files.
