// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
//   Program:   FITS Reader for ITK
//   Module:    itkFITSImageIO.cxx
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
//      http://www.opensource.org/licenses/mit-license
//
//   for details.
//
//=============================================================================


#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <math.h>
#include <zlib.h>

#include <itkExceptionObject.h>
#include <itkByteSwapper.h>
#include <itkMetaDataObject.h>
#include <itksys/RegularExpression.hxx>

#include <wcs.h>
#include <da_usual.h>

#include "itkFITSImageIO.h"
#include <grparser.h> // for FITS NGP_MAX_ARRAY_DIM

extern "C" { char* GetFITShead(const char* filepath, bool); }

using std::string;
using std::stringstream;
using std::ostream;
using std::cerr;
using std::endl;
using std::vector;


//*****************************************************************************
//*****              Local procedures and classes                         *****
//*****************************************************************************

//-----------------------------------------------------------------------------
// debugPrint(): local macro
//-----------------------------------------------------------------------------

#define debugPrint(message) \
   { if (_cv_debugLevel) { \
        cerr << message << endl; \
     } \
   }

//-----------------------------------------------------------------------------
// max(): local proc
//-----------------------------------------------------------------------------

/*local proc*/ template <class T> T
max(T a, T b) {
  if (a > b) return a;
  else return b;
}


// //--------------------------------------------------------------------------
// // toLower(): local proc
// //--------------------------------------------------------------------------

// // Converts \a s to lowercase.

// /*local proc*/ void
// toLower(string& s)
// {
//    for (string::iterator p = s.begin(); p != s.end( ); ++p) {
//      *p = tolower(*p);
//    }

// }
// //--------------------------------------------------------------------------
// // endMatchesP(): local proc
// //--------------------------------------------------------------------------

// // Returns true iff \a extension is on the end of \a filepath.

// /*local proc*/ static bool
// endMatchesP(const string& filepath, const string& extension)
// {
//   typedef string::size_type StrLen;
//   const StrLen extensionLength = extension.length();
//   const StrLen filepathLength = filepath.length();
//   if (extensionLength >= filepathLength) return false;
//   string filepathEnd = filepath.substr(filepathLength - extensionLength);
//   if (filepathEnd == extension) return true;
//   else return false;
// }


//-----------------------------------------------------------------------------
// checkExtension(): local proc
//-----------------------------------------------------------------------------

// Returns true iff \a filepath ends with an extension that indicates
// that it is a FITS file.

/*local proc*/ static bool
checkExtension(const string& filepath)
{
  static itksys::RegularExpression fitsRE =
    "^.+\\.(fits?|imh)(.(gz|Z))?"
    "(\\([^()]+\\))?((\\[[^][]+\\])*|\\+[0-9]+)$";

  static itksys::RegularExpression rawRE = 
    "^.+\\.dat\\[[^][]+\\]$";

  return fitsRE.find(filepath) or rawRE.find(filepath);

// TODO: Remove the following commented-out code:

//   // TODO: Extend this to deal CFITSIO's "Extended File Name Syntax".

//   if (endMatchesP(loweredFilepath, ".fits")    or
//       endMatchesP(loweredFilepath, ".fits.gz") or
//       endMatchesP(loweredFilepath, ".fit")     or
//       endMatchesP(loweredFilepath, ".fit.gz")  or
//       endMatchesP(loweredFilepath, ".fits.Z")  or
//       endMatchesP(loweredFilepath, ".fit.Z"))
//     {
//       return true;
//     }
//   else return false;
}


//-----------------------------------------------------------------------------
// getAllFitsErrorMessages(): local proc
//-----------------------------------------------------------------------------

/*local proc*/ string
getAllFitsErrorMessages(const int status)
{
  string retval;
  char fitsError[FLEN_ERRMSG];
//   while (::fits_read_errmsg(fitsError)) {
//     retval += fitsError;
//   }
  retval = "FITSIO status = ";
  stringstream ss;
  ss << status;
  retval += ss.str();
  retval += ": ";
  ::fits_get_errstatus(status, fitsError);
  retval += fitsError;
  return retval;
}


//-----------------------------------------------------------------------------
// square(): local inline proc
//-----------------------------------------------------------------------------

/*local proc*/ inline double
square(double x)
{
  return x * x;
}


//*****************************************************************************
//*****                                                                   *****
//*****         FITSImageIO: leaf subclass of ImageIOBase                 *****
//*****                                                                   *****
//*****************************************************************************

namespace itk {

//-----------------------------------------------------------------------------
// Private class variables
//-----------------------------------------------------------------------------

int FITSImageIO::_cv_debugLevel = 0;
double FITSImageIO::_cv_nullValue = 0.0; // The default of 0.0 causes NaN's
                                         // to be left as NaN's, rather than
                                         // converted to 0.0, as one would
                                         // naively expect.
double FITSImageIO::_cv_rotateSky = 0;     
double FITSImageIO::_cv_rotateDecIntoVelocityAxis = 0;
double FITSImageIO::_cv_rotateRAIntoVelocityAxis = 0;
double FITSImageIO::_cv_scaleVoxelValues = 1;
double FITSImageIO::_cv_scaleRA = 1;
double FITSImageIO::_cv_scaleDec = 1;
bool   FITSImageIO::_cv_autoScaleVelocityAxis = false;
double FITSImageIO::_cv_scaleVelocityAxis = 1;
double FITSImageIO::_cv_scaleAllAxes = 1;


//=============================================================================
//                    Private methods                                   
//=============================================================================

//-----------------------------------------------------------------------------
// getFitsHeader(): private method of FITSImageIO
//-----------------------------------------------------------------------------

/*private method*/ string
FITSImageIO::getFitsHeader()
{
  int status = 0;
  const bool noComments = true;
  char* headersToExclude[0] = { };
  const int nHeadersToExclude = 0;
  char* retval1;
  int nKeysDummy;
  ::fits_hdr2str(m_fitsFile, !noComments, headersToExclude, nHeadersToExclude,
                 &retval1, &nKeysDummy, &status);
  DaFreer freer (retval1);
  if (status) {
    itkExceptionMacro("FITSImageIO could not get header from Primary Array of"
                      "FITS file: \""
                      << this->GetFileName() << "\": "
                      << ::getAllFitsErrorMessages(status) << '.');
  }
  string retval (retval1);
  return retval;
}


//=============================================================================
//                    Public methods                                   
//=============================================================================

//-----------------------------------------------------------------------------
// CanReadFile(): inherited virtual method
//-----------------------------------------------------------------------------

/*method*/ bool
FITSImageIO::CanReadFile(const char* const filepath) 
{ 
  debugPrint("Entering FITSImageIO::CanReadFile().")

  if (!filepath or *filepath == 0) {
    itkDebugMacro(<< "No filename specified.");
    return false;
  } else if (::checkExtension(filepath)) {
    debugPrint("Exiting FITSImageIO::CanReadFile() with true.");
    return true;
  } else {
    itkDebugMacro(<<"The filename extension is not recognized");
    return false;
  }
}


//-----------------------------------------------------------------------------
// CanWriteFile(): inherited virtual method
//-----------------------------------------------------------------------------

/*method*/ bool
FITSImageIO::CanWriteFile(const char* const name)
{
  
  // Return false because we currently don't implement writing FITS files:
  return false;

// TODO: More fully implement the following if we ever need to be able
// to write FITS files:

//   const string filepath = name;
//   if (filepath == "") {
//     itkDebugMacro(<< "No filename specified.");
//   }

//   bool extensionFound = ::checkExtension(name);

//   if(!extensionFound) {
//     itkDebugMacro(<<"The filename extension is not recognized");
//     return false;
//   }

//   if (extensionFound) {
//     return true;
//   }
//   return false;

}


//-----------------------------------------------------------------------------
// ReadImageInformation(): inherited virtual method
//-----------------------------------------------------------------------------

/*method*/ void
FITSImageIO::ReadImageInformation()
{
  // CFITSIO writes to this to indicate errors:
  int status = 0;   

  // Open the FITS file:
  ::fits_open_file(&m_fitsFile, this->GetFileName(), READONLY, &status);
  if (status) {
    itkExceptionMacro("FITSImageIO could not open FITS file: \""
                      << this->GetFileName() << "\" for reading: "
                      << ::getAllFitsErrorMessages(status) << '.');
  }

  debugPrint("Entering FITSImageIO::ReadImageInformation().");

  // Get the dimensions and type of the FITS Primary Array:
  int numOfAxes;
  long lengthOfAxisInPixels[NGP_MAX_ARRAY_DIM];
  int bitsPerPixel;
  ::fits_get_img_param(m_fitsFile, NGP_MAX_ARRAY_DIM, &bitsPerPixel,
                       &numOfAxes, lengthOfAxisInPixels, &status);
  if (status) {
    itkExceptionMacro("FITSImageIO could not read Primary Array parameters "
                      "from FITS file \""
                      << this->GetFileName() << "\":"
                      << ::getAllFitsErrorMessages(status) << '.');
  }

  debugPrint("NAXIS=" << numOfAxes);
  if (_cv_debugLevel) {
    cerr << "DIMS=";
    for (long i = 0; i < numOfAxes; ++i) {
      cerr << " " << lengthOfAxisInPixels[i];
    }
    cerr << endl;
  }
  debugPrint("BPP= " << bitsPerPixel);

  if (numOfAxes != 3) {
    itkExceptionMacro("FITSImageIO cannot handle FITS Primary Arrays that are"
                      " anything other than three dimensional.");
    
    // TODO: Remove this restriction?
  }

  // BEGIN getting RA and Dec of border pixels.

  // TODO: Check to make sure that wcsinit copies the string handed to it by
  // fitsHeader.c_str(), as I believe that the string returned by
  // fitsHeader.c_str() will go stale as soon as fitsHeader goes out of scope.

  string fitsHeader = getFitsHeader();
  WorldCoor* wcs = wcsinit(fitsHeader.c_str());
  const ConstRcMallocPointer<WorldCoor> wcsRcPtr = wcs;
  typedef itk::FITSWCSTransform<double, 3> WCSTransform;
  m_transform = WCSTransform::New();
  m_transform->SetWCS(wcsRcPtr);
  m_transform->ApplySettings();

  // TODO: Delete this commented out code:
//   Freer freer (wcs);

  double lowerLeftRA, lowerLeftDec;
  pix2wcs(wcs, 1, 1, &lowerLeftRA, &lowerLeftDec);
  double lowerRightRA, lowerRightDec;
  pix2wcs(wcs, lengthOfAxisInPixels[0], 1, &lowerRightRA, &lowerRightDec);
  double upperLeftRA, upperLeftDec;
  pix2wcs(wcs, 1, lengthOfAxisInPixels[1], &upperLeftRA, &upperLeftDec);

  debugPrint("Before correcting for equinox crossing:");
  debugPrint("LL: RA=" << lowerLeftRA << ' ' << "Dec=" << lowerLeftDec);
  debugPrint("UL: RA=" << upperLeftRA << ' ' << "Dec=" << upperLeftDec);
  debugPrint("LR: RA=" << lowerRightRA << ' ' << "Dec=" << lowerRightDec);

  // If we are in debug output mode, then test out the FITSWCSTransform object:
  if (_cv_debugLevel) {
    WCSTransform::Pointer inverseTransform = WCSTransform::New();
    m_transform->GetInverse(inverseTransform);
    WCSTransform::InputPointType ijkPoint;
    WCSTransform::OutputPointType wcsPoint;
    ijkPoint[0] = 1;
    ijkPoint[1] = lengthOfAxisInPixels[1];
    wcsPoint = m_transform->TransformPoint(ijkPoint);
    cerr << "Value of transformed UL point: " << wcsPoint << endl;
    ijkPoint[0] = -666;
    ijkPoint[1] = -666;
    ijkPoint[2] = -666;
    ijkPoint = inverseTransform->TransformPoint(wcsPoint);
    cerr << "Value of UL point tranformed back to ijk coordinates: "
	 << ijkPoint << endl;    
      
  }

  // Make some RAs negative (by substracting 360 degrees) if the image crosses
  // the equinox.  Otherwise, we will incorrectly think that we are
  // representing a huge area of the sky, instead of a small area:
  const double xRaDiff = lowerRightRA - lowerLeftRA;
  if (xRaDiff > 180.0) lowerRightRA -= 360.0;
  else if (xRaDiff < -180.0) lowerLeftRA -= 360.0;

  const double yRaDiff = upperLeftRA - lowerLeftRA;
  if (yRaDiff > 180.0) upperLeftRA -= 360.0;
  else if (yRaDiff < -180.0) lowerLeftRA -= 360.0;

  debugPrint("After correcting for equinox crossing:");
  debugPrint("LL: RA=" << lowerLeftRA << ' ' << "Dec=" << lowerLeftDec);
  debugPrint("UL: RA=" << upperLeftRA << ' ' << "Dec=" << upperLeftDec);
  debugPrint("LR: RA=" << lowerRightRA << ' ' << "Dec=" << lowerRightDec);

  // TODO: The above method for removing an RA discontinuity should be improved
  // for large areas of the sky or for areas near the poles (in the unlikely
  // case that we were to continue to take this sort of approach at all), as we
  // would have to determine the direction that RA is moving along each axis by
  // looking at a small increment along the axis to make sure that the total
  // change in RA from one side of the image to the other is not more than 180
  // degrees.  If we want to handle large areas of the sky, we will also have
  // to do something similar for Dec.

  // END getting RA and Dec of border pixels.


  // The following is an explanation of the math involved in the coordinate
  // transformation that is to come immediately below: One way of transforming
  // between two different Cartesian coordinate systems (let's call the source
  // coordinate system A and the destination coordinate system B) for the same
  // Euclidean space involves knowing two things:

  // 1. The coordinates in B corresponding to the origin in A.  Let's call
  // these coordinates in B, "o".

  // 2. A change-of-basis matrix that transforms the unit basis vectors for A
  // into (unnormalized) basis vectors for B.  Let's call this change-of-basis
  // matrix, "C".

  // Knowing C and o, one can then transform the coordinates in A for any given
  // point (let's call the coordinates in A for this given point, "c") into the
  // corresponding coordinates in B by calculating C * c + o.

  // As it turns out, ITK uses just such a method for transforming coordinate
  // systems, but rather than accepting o and C as the parameters for this
  // transformation, it wants us to first factor C into a "direction cosine
  // matrix" (let's call this matrix "D") and a spacing vector (let's call the
  // spacing vector "s").  After factoring C into D and s, D consists of
  // normalized basis vectors (i.e., unit vectors) for B and s consists of the
  // length of each unnormalized vector in C.  The transformation from A to B
  // for a given point is then calculated as (D * s) * c + o.

  // TODO: RA and Dec represent a polar coordinate system, not a Cartesian
  // coordinate system, and hence this method sucks for fields of view that are
  // not small or are near the poles.  Unfortunately, for the time being, we
  // have to work with the facilities that 3D Slicer provides us, and so this
  // is the best that we can do for right now.

  // BEGIN Calculate the change-of-basis matrix, "C":
  const double raPerI =
    (lowerRightRA - lowerLeftRA) / (lengthOfAxisInPixels[0] - 1);
  const double decPerI =
    (lowerRightDec - lowerLeftDec)/(lengthOfAxisInPixels[0] - 1);
  const double raPerJ =
    (upperLeftRA - lowerLeftRA) / (lengthOfAxisInPixels[1] - 1);
  const double decPerJ = 
    (upperLeftDec - lowerLeftDec)/(lengthOfAxisInPixels[1] - 1);

  debugPrint("raPerI=" << raPerI);
  debugPrint("decPerI=" <<  decPerI);
  debugPrint("raPerJ=" << raPerJ);
  debugPrint("decPerJ=" << decPerJ);

  double velocityPerK;

  if (_cv_autoScaleVelocityAxis) {

    // Calculate length of each side of the sky rectangle in RA/Dec
    // coordinates:
    const double rectWidth = sqrt(square(lowerRightRA - lowerLeftRA) +
				  square(lowerRightDec - lowerLeftDec));
    const double rectHeight = sqrt(square(upperLeftRA - lowerRightRA) +
				   square(upperLeftDec - lowerRightDec));
    debugPrint("rectWidth= " << rectWidth);
    debugPrint("rectHeight= " << rectHeight);

    velocityPerK = max(rectWidth, rectHeight) / lengthOfAxisInPixels[2]
      * _cv_scaleVelocityAxis;

  } else {

    // TODO: This is not the right value at all!  It really needs to be
    // determinded from the FITS headers (and then scaled by
    // _cv_scaleVelocityAxis):
    velocityPerK = _cv_scaleVelocityAxis;
  }

  debugPrint("velocityPerK=" << velocityPerK);

  // TODO: Extend the following to more than three dimensions.

  // TODO: We need to extract velocity information using a WCS library of some
  // sort.  Is there one that does this?  At the moment, I just set velocity to
  // nonsensical things, like 0 for the image origin, immediately below.

  // Create the origin vector (what we called "o" in the exposition above):
  const double origin[3] = { lowerLeftRA, 0, lowerLeftDec };

  // Create the change-of-basis matrix (what we called "C" in the exposition
  // above):
  vector<double> changeOfBasisMatrix[3];
  for (int axis = 0; axis < 3; ++axis) {
    changeOfBasisMatrix[axis].resize(3);
  }

  // Each column of the matrix represents a vector that transform a
  // unit vector in i, j, k space to RA, V, DEC (a.k.a. LPS) space:

  const int ra  = 0;
  const int vel = 1;
  const int dec = 2;

  const int i   = 0;
  const int j   = 1;
  const int k   = 2;

  changeOfBasisMatrix[ra] [i] = raPerI * _cv_scaleRA;
  changeOfBasisMatrix[vel][i] = 0;
  changeOfBasisMatrix[dec][i] = decPerI;

  changeOfBasisMatrix[ra] [j] = raPerJ * _cv_scaleRA;
  changeOfBasisMatrix[vel][j] = 0;
  changeOfBasisMatrix[dec][j] = decPerJ;

  changeOfBasisMatrix[ra] [k] = 0;
  changeOfBasisMatrix[vel][k] = velocityPerK;
  changeOfBasisMatrix[dec][k] = 0;


  // END Calculate the change-of-basis matrix.

  if (_cv_debugLevel) {
    cerr << "Change-of-basis matrix:\n";
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
	cerr << "    " << changeOfBasisMatrix[row][col];
      }
      cerr << endl;
    }
  }

  // TODO: Extract the dimensions for the velocity spectrum from the FITS 
  // file rather than punting and just setting it to 1 as we have done above.

  // Create a direction cosine matrix and spacing vector from the
  // change-of-basis matrix.  The direction cosines are calculated by
  // normalizing the direction vectors contained in the change-of-basis matrix
  // (i.e., dividing their values by , while also populating the spacing
  // vector:
  double spacing[3];
  vector<double> directionCosines[3];

//   for (int i = 0; i < 3; ++i) {
//     directionCosines[i].resize(3);
//     spacing[i] = sqrt(square(changeOfBasisMatrix[i][0]) + 
//                       square(changeOfBasisMatrix[i][1]) +
//                       square(changeOfBasisMatrix[i][2]));
//     for (int j = 0; j < 3; ++j) {
//       directionCosines[i][j] = changeOfBasisMatrix[i][j] / spacing[i];
//     }
//   }


  for (int indexAxis = 0; indexAxis < 3; ++indexAxis) {
    directionCosines[indexAxis].resize(3);
    spacing[indexAxis] = sqrt(square(changeOfBasisMatrix[0][indexAxis]) + 
			      square(changeOfBasisMatrix[1][indexAxis]) +
			      square(changeOfBasisMatrix[2][indexAxis]));
  }
  for (int indexAxis = 0; indexAxis < 3; ++indexAxis) {
    for (int physicalAxis = 0; physicalAxis < 3; ++physicalAxis) {
      directionCosines[indexAxis][physicalAxis] =
	changeOfBasisMatrix[physicalAxis][indexAxis] / spacing[indexAxis];
    }
  }

  // TODO: Adapt the following code to do rotations, using the -R
  // option, which you already partially implemented the command-line
  // parsing for.

//   For testing purposes the following code rotates the image by 30 degrees,
//   to make sure that the direction cosines are being interpreted properly:

//    double angle = PI/6; //d
//    directionCosines[0][0] = cos(angle); //d
//    directionCosines[0][1] = sin(angle); //d
//    directionCosines[1][0] = -sin(angle); //d
//    directionCosines[1][1] = cos(angle); //d
//    spacing[0] = 1; //d
//    spacing[1] = 1; //d
//    spacing[2] = 1; //d


  for (int indexAxis = 0; indexAxis < 3; ++indexAxis) {
    spacing[indexAxis] *= _cv_scaleAllAxes;
  }

  // Set up the ITK image:
  this->SetNumberOfComponents(1);
  this->SetPixelType(SCALAR);
  this->SetComponentType(FLOAT);
  this->SetNumberOfDimensions(numOfAxes);
  for (int indexAxis = 0; indexAxis < numOfAxes; ++indexAxis) {
    this->SetDimensions(indexAxis, lengthOfAxisInPixels[indexAxis]);
    this->SetOrigin(indexAxis, origin[indexAxis]);
    this->SetSpacing(indexAxis, spacing[indexAxis]);
    this->SetDirection(indexAxis, directionCosines[indexAxis]);

    // Write debugging output if debug output option is set:
    if (_cv_debugLevel) {
      cerr << "spacing=" << indexAxis << "," << spacing[indexAxis] << " ";
      cerr << "directions=";
      for (int physicalAxis = 0; physicalAxis < 3; ++physicalAxis) {
	cerr << directionCosines[indexAxis][physicalAxis] << " ";
      }
      cerr << endl;
    }
  }


  // Put the FITS Primary Array Header into the ITK MetaDataDictionary as one
  // big string:
  MetaDataDictionary& dict = this->GetMetaDataDictionary();

  itk::EncapsulateMetaData<string>(dict, "FITS Header", fitsHeader);

  // Also break up the aforementioned FITS header into individual entries and
  // add each entry into the MetaDataDictionary as a separate entry:
  int nKeys = 0;
  int dummy;
  ::fits_get_hdrspace(m_fitsFile, &nKeys, &dummy, &status);
  for (int keyIndex = 1; keyIndex <= nKeys; ++keyIndex) {
    char keyName[FLEN_KEYWORD];
    char keyValue[FLEN_VALUE];
    char keyComment[FLEN_COMMENT];
    ::fits_read_keyn(m_fitsFile, keyIndex, keyName, keyValue, keyComment,
                     &status);
    itk::EncapsulateMetaData<string>(dict,
				     string("FITS.") + keyName,
				     keyValue);

    // TODO: You need to put the comments and units somewhere too.
  }

  // TODO: Check status and raise exception.
  

  // TODO: Above we put the uninterpreted FITS Primary Array header into the
  // ITK image as single huge value.  But, in addition to doing that, we should
  // also parse all the entries and make a separate ITK MetaDataDictionary
  // entry for each FITS header entry.  To do so, wade through the comments and
  // commented-out code below, and turn it into working code:

  // Documentation regarding CFITSIO to parse the FITS headers:
  //
  //
  //    CFITSIO Quick Start Guide, p. 4: How to get the list of tags in the
  //    header.
  //
  //    CFITSIO User's Guide, p. 34 
   
// typedef string defaultValueType;  

//    // Get the header
//    int numberOfHeaderEntries = 0;
//    ::fits_get_hdrspace(fptr, &numberOfHeaderEntries, NULL, &status);

//    // TODO: check this number by looking at cfits code...
//    const unsigned int maxRecordLength = 1024;

//    char charKey[ maxRecordLength ];
//    char charValue[ maxRecordLength ];
//    char* dontNeedTheComment = NULL;

//    // Get the records and populate the MetaDataDictionary 
//    MetaDataDictionary& thisDic = this->GetMetaDataDictionary();

//    for(unsigned int keyentry=0; keyentry < numberOfHeaderEntries; keyentry++)
//      { 
//        unsigned int keyentryFortran = keyentry + 1;

//        // Get the record from cfits
//        ::fits_read_keyn(fptr, keyentry, charKey, charValue,
//                      dontNeedTheComment, &status);
                        

//        if (status != 0) {
//         break;  // TODO: report an error
//        }
 
//       string stringKey   = charKey;
//       string stringValue = charValue;

//       // Put the record in the MetaDataDictionary
//       itk::EncapsulateMetaData<defaultValueType>(thisDic, stringKey,
//                                               stringValue);
//     }

}


//-----------------------------------------------------------------------------
// Read(): inherited virtual method
//-----------------------------------------------------------------------------

/*method*/ void
FITSImageIO::Read(void* const buffer)
{
  // TODO: At some point we might want to preserve the native pixel
  // type, rather than converting everything into a float.  Achieving
  // this is complicated, however, by the fact that ITK image types
  // are parameterized by the pixel type.

  float* bufferAsFloats = static_cast<float*>(buffer);

  // Make a pixel location array that represents the origin pixel; i.e., [1,
  // 1, 1, ...]:
  const unsigned nDimensions = this->GetNumberOfDimensions();
  long origin[nDimensions];
  for (int i = 0; i < nDimensions; ++i) {
    origin[i] = 1;
  }

  // Calculate the total number of pixels in the image:
  long nPixels = 1;
  for (int i = 0; i < nDimensions; ++i) {
    nPixels *= GetDimensions(i);
  }

  // Read the FITS image into the ITK image buffer:
  int status = 0;
  const float nullValue = _cv_nullValue;
  const float* const nullValuePtr = &nullValue;
  int anyNull = false;
  ::fits_read_pix(m_fitsFile, TFLOAT, origin, nPixels, (void*) nullValuePtr,
                  bufferAsFloats, &anyNull, &status);
  if (nullValue != 0) {
    debugPrint("Any null values? " << anyNull);
    debugPrint("Null value will be fetched as: " << nullValue);
  }
  if (status) {
    itkExceptionMacro("FITSImageIO::Read() could not read the primary data "
                      "array from FITS file \""
                      << this->GetFileName() << "\": "
                      << ::getAllFitsErrorMessages(status) << ".");
  }
                

  // Scale the voxel values by the voxel value scaling factor:
  if (_cv_scaleVoxelValues != 1) {
    const float* const lastPixel = bufferAsFloats + nPixels;
    for (float* pixelPtr = bufferAsFloats; pixelPtr < lastPixel; ++pixelPtr) {
      *pixelPtr = _cv_scaleVoxelValues * *pixelPtr;
    }
  }

  // Close the FITS file:
  ::fits_close_file(m_fitsFile, &status);
  if (status) {
    itkExceptionMacro("FITSImageIO::Read() could not close FITS file \""
                      << this->GetFileName()
                      << "\" after reading the primary data array: "
                      << ::getAllFitsErrorMessages(status) << ".");
  }
}


//-----------------------------------------------------------------------------
// WriteImageInformation(): inherited virtual method
//-----------------------------------------------------------------------------

//! Not yet implemented.

/*method*/ void 
FITSImageIO::WriteImageInformation()
{
  // TODO: Implement with help of cfitsio
  itkExceptionMacro("FITSImageIO::WriteImageInformation() not implemented "
		    "yet.");
}


//-----------------------------------------------------------------------------
// Write(): inherited virtual method
//-----------------------------------------------------------------------------

//! Not yet implemented.

/*method*/ void 
FITSImageIO::Write(const void* const buffer) 
{
  // TODO: Maybe implement this someday.

  itkExceptionMacro("FITSImageIO::Write() not implemented.");
}


//-----------------------------------------------------------------------------
// PrintSelf(): inherited virtual method
//-----------------------------------------------------------------------------

/*method*/ void
FITSImageIO::PrintSelf(ostream& os, const Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PixelType " << m_PixelType << "\n";
}

} // namespace itk
