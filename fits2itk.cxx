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
  cerr << "\nVersion 0.1dev.0\n\n";
  cerr << "usage: ";
  if (daProgramName().size()) cerr << basename(daProgramName().c_str());
  else cerr << "fits2itk";
  cerr << 
    " [-ASU] [-s pixel-scale] [-v velocity-scale] \n"
    "         input-file output-file\n"
    "\n"
    "  A: auto-scale velocity axis\n"
    "  D: debug level\n"
    "  S: coerce pixel values to shorts\n"
    "  U: coerce pixel values to unsigned shorts\n"
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
  const char options[] = "AD:RSUv:s:";
  while ((optionChar = ::getopt(argc, const_cast<char**>(argv), options))
         != -1)
    {
      switch (optionChar) {

      case 'A':
        itk::FITSImageIO::SetAutoScaleVelocityAxis(true);
        break;

      case 'D':
	itk::FITSImageIO::SetDebugLevel(strtol(optarg, NULL, 10));
	break;

      case 'R':
	// TODO: This isn't implemented yet.
	itk::FITSImageIO::SetRotateSky(strtod(optarg, NULL));
	break;

      case 'S':
        coerceToShorts = true;
        break;

      case 's':
        itk::FITSImageIO::SetScaleVoxelValues(strtod(optarg, NULL));
        break;

      case 'U':
        coerceToUnsignedShorts = true;
        break;

      case 'v':
        // TODO: Make the following (and other similar code here) more
        // robust:
        itk::FITSImageIO::SetScaleVelocityAxis(strtod(optarg, NULL));
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


//-----------------------------------------------------------------------------
// Change Log
//-----------------------------------------------------------------------------

// See http://wiki.iic.harvard.edu/AstroMed/Software_Version_Numbers for
// documentation on how these version numbers work.

//------------------------
// Version 0.1dev.0
//------------------------

// *** Tue Dec  5, 2006 ***

// Added version number to usage output and a change log comments section.

// Changed autoscaling of velocity axis to be average of raPerI and decPerI,
// rather than just raPerI.

// RA/DEC/V were getting mapped to L/P/S rather than to L/P/S.  Fixed this.

// Added option to turn on debugging messages and turned them off by default.
// They were all getting output all the time prior to this change.

//-------------------------------------------------
// Changes described above this line are checked in
//-------------------------------------------------

