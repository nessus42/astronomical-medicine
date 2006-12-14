/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFITSImageIO.cxx,v $
  Language:  C++
  Date:      $Date: 2005/10/10 19:33:50 $
  Version:   $Revision: 1.19 $  

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// TODO: Fix the above copyright, which maybe should be IIC, and the RCS stuff,
// which is not right for Subversion.


#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <math.h>
#include <zlib.h>

#include <itkExceptionObject.h>
#include <itkByteSwapper.h>
#include <itkMetaDataObject.h>
#include <itksys/SystemTools.hxx>

#include <wcs.h>

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
//***              Local procedures and classes                             ***
//*****************************************************************************

//-----------------------------------------------------------------------------
// debugOutput(): local macro
//-----------------------------------------------------------------------------

#define debugOutput(message) \
   { if (_cv_debugLevel) { \
        cerr << message << endl; \
     } \
   }

//-----------------------------------------------------------------------------
// Freer: local class
//-----------------------------------------------------------------------------

// Frees a malloced object via RAII.

class Freer {
  void* _malloced_object;

  // Disable copying:
  Freer(const Freer&);           
  void operator=(const Freer&);

public:
  Freer(void* malloced_object) : _malloced_object(malloced_object) {}
  ~Freer() { free((char *)_malloced_object); }
};


//-----------------------------------------------------------------------------
// toLower(): local proc
//-----------------------------------------------------------------------------

// Converts \a s to lowercase.

/*local proc*/ void
toLower(string& s)
{
   for (string::iterator p = s.begin(); p != s.end( ); ++p) {
     *p = tolower(*p);
   }
}


//-----------------------------------------------------------------------------
// endMatchesP(): local proc
//-----------------------------------------------------------------------------

// Returns true iff \a extension is on the end of \a filepath.

/*local proc*/ static bool
endMatchesP(const string& filepath, const string& extension)
{
  typedef string::size_type StrLen;
  const StrLen extensionLength = extension.length();
  const StrLen filepathLength = filepath.length();
  string filepathEnd = filepath.substr(filepathLength - extensionLength);
  if (filepathEnd == extension) {
    return true;
  }
  else return false;
}


//-----------------------------------------------------------------------------
// checkExtension(): local proc
//-----------------------------------------------------------------------------

// Returns true iff \a filepath ends with an extension that indicates
// that it is a FITS file.

/*local proc*/ static bool
checkExtension(const string& filepath)
{
  // TODO: Extend this to deal CFITSIO's "Extended File Name Syntax".

  string loweredFilepath = filepath;
  toLower(loweredFilepath);
  if (endMatchesP(loweredFilepath, ".fits")    or
      endMatchesP(loweredFilepath, ".fits.gz") or
      endMatchesP(loweredFilepath, ".fit")     or
      endMatchesP(loweredFilepath, ".fit.gz")  or
      endMatchesP(loweredFilepath, ".fits.Z")  or
      endMatchesP(loweredFilepath, ".fit.Z"))
    {
      return true;
    }
  else return false;
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
// square() local inline proc
//-----------------------------------------------------------------------------

/*local proc*/ inline double
square(double x)
{
  return x * x;
}


//*****************************************************************************
//***                 Private Class Variables                               ***
//*****************************************************************************

namespace itk {

int FITSImageIO::_cv_debugLevel = 0;
double FITSImageIO::_cv_scaleVoxelValues = 1;
double FITSImageIO::_cv_scaleRA = 1;
double FITSImageIO::_cv_scaleDec = 1;
double FITSImageIO::_cv_scaleVelocityAxis = 1;
bool   FITSImageIO::_cv_autoScaleVelocityAxis = false;                  
double FITSImageIO::_cv_rotateSky = 0;
double FITSImageIO::_cv_rotateDecIntoVelocityAxis = 0;
double FITSImageIO::_cv_rotateRAIntoVelocityAxis = 0;


//*****************************************************************************
//***                     Private Methods                                   ***
//*****************************************************************************

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
  Freer freer (retval1);
  if (status) {
    itkExceptionMacro("FITSImageIO could not get header from Primary Array of"
                      "FITS file: \""
                      << this->GetFileName() << "\": "
                      << ::getAllFitsErrorMessages(status) << '.');
  }
  string retval (retval1);
  return retval;
}


//*****************************************************************************
//***                     Public Methods                                    ***
//*****************************************************************************

//-----------------------------------------------------------------------------
// CanReadFile(): virtual method of FITSImageIO
//   implements pure virtual method of ImageIOBase
//-----------------------------------------------------------------------------

/*method*/ bool
FITSImageIO::CanReadFile(const char* const filepath) 
{ 
  // This is set this way merely for the purpose of trying to see what
  // ITK will and won't allow us to read:
  return true;  //d

  debugOutput("Entering FITSImageIO::CanReadFile().")

  if (!filepath or *filepath == 0) {
    itkDebugMacro(<< "No filename specified.");
    return false;
  } else if (::checkExtension(filepath)) {
    debugOutput("Exiting FITSImageIO::CanReadFile() with true.");
    return true;
  } else {
    itkDebugMacro(<<"The filename extension is not recognized");
    return false;
  }
}


//-----------------------------------------------------------------------------
// CanWriteFile(): virtual method of FITSImageIO
//   implements pure virtual method of ImageIOBase
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
// ReadImageInformation(): virtual method of FITSImageIO
//   implements pure virtual method of ImageIOBase
//-----------------------------------------------------------------------------

//! Read information about the FITS file and put the cursor of the
//! stream just before the first data pixel.

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

  debugOutput("Entering FITSImageIO::ReadImageInformation().");

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

  debugOutput("NAXIS=" << numOfAxes);
  if (_cv_debugLevel) {
    cerr << "DIMS=";
    for (long i = 0; i < numOfAxes; ++i) {
      cerr << " " << lengthOfAxisInPixels[i];
    }
    cerr << endl;
  }
  debugOutput("BPP= " << bitsPerPixel);

  if (numOfAxes != 3) {
    itkExceptionMacro("FITSImageIO cannot handle FITS Primary Arrays that are"
                      " anything other than three dimensional.");
    
    // TODO: Remove this restriction?
  }

  // BEGIN getting RA and Dec of border pixels.

  string fitsHeader = getFitsHeader();
  WorldCoor* const wcs = wcsinit(fitsHeader.c_str());
  Freer freer (wcs);
  double lowerLeftRA, lowerLeftDec;
  pix2wcs(wcs, 1, 1, &lowerLeftRA, &lowerLeftDec);
  double lowerRightRA, lowerRightDec;
  pix2wcs(wcs, lengthOfAxisInPixels[0], 1, &lowerRightRA, &lowerRightDec);
  double upperLeftRA, upperLeftDec;
  pix2wcs(wcs, 1, lengthOfAxisInPixels[1], &upperLeftRA, &upperLeftDec);

  debugOutput("LL: RA=" << lowerLeftRA << ' ' << "Dec=" << lowerLeftDec);
  debugOutput("UL: RA=" << upperLeftRA << ' ' << "Dec=" << upperLeftDec);
  debugOutput("LR: RA=" << lowerRightRA << ' ' << "Dec=" << lowerRightDec);

  // Make some RAs negative (by substracting 360 degrees) if the image crosses
  // the equinox.  Otherwise, we will incorrectly think that we are
  // representing a huge area of the sky, instead of a small area:
  const double xRaDiff = lowerRightRA - lowerLeftRA;
  if (xRaDiff > 180.0) lowerRightRA -= 360.0;
  else if (xRaDiff < -180.0) lowerLeftRA -= 360.0;

  const double yRaDiff = upperLeftRA - lowerLeftRA;
  if (yRaDiff > 180.0) upperLeftRA -= 360.0;
  else if (yRaDiff < -180.0) lowerLeftRA -= 360.0;

  debugOutput("LL: RA=" << lowerLeftRA << ' ' << "Dec=" << lowerLeftDec);
  debugOutput("UL: RA=" << upperLeftRA << ' ' << "Dec=" << upperLeftDec);
  debugOutput("LR: RA=" << lowerRightRA << ' ' << "Dec=" << lowerRightDec);

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

  // Calculate the change-of-basis matrix:
  const double raPerI = (lowerRightRA - lowerLeftRA)/(lengthOfAxisInPixels[0]-1);
  const double decPerI =
    (lowerRightDec - lowerLeftDec)/(lengthOfAxisInPixels[0]-1);
  const double raPerJ = (upperLeftRA - lowerLeftRA)/(lengthOfAxisInPixels[1]-1);
  const double decPerJ = 
    (upperLeftDec - lowerLeftDec)/(lengthOfAxisInPixels[1]-1);

  debugOutput("raPerI=" << raPerI);
  debugOutput("decPerI=" <<  decPerI);
  debugOutput("raPerJ=" << raPerJ);
  debugOutput("decPerJ=" << decPerJ);

  // TODO: Extend the following to more than three dimensions.

  // YOU ARE HERE: trying to figure out how RA/DEC/V should map onto LPS space.

  // TODO: We need to extract velocity information using a WCS library of some
  // sort.  Is there one that does this?  At the moment, I just set velocity to
  // nonsensical things, like 0 for the image origin, immediately below.

  // Create the origin vector (what we called "o" in the exposition above):
  const double origin[3] = { lowerLeftRA, 0, lowerLeftDec };

  // Create the change-of-basis matrix (what we called "C" in the exposition
  // above):
  vector<double> changeOfBasisMatrix[3];
  for (int axis = 0; axis < 3; ++axis) changeOfBasisMatrix[axis].resize(3);

  changeOfBasisMatrix[0][0] = raPerI;
  changeOfBasisMatrix[0][1] = decPerI;
  changeOfBasisMatrix[0][2] = 0;

  changeOfBasisMatrix[1][0] = 0;
  changeOfBasisMatrix[1][1] = 0;
  changeOfBasisMatrix[1][2] = _cv_autoScaleVelocityAxis
    ? (sqrt(square(raPerI) + square(raPerJ)) +
       sqrt(square(decPerI) + square(decPerJ))) / 2
    : _cv_scaleVelocityAxis;

  changeOfBasisMatrix[2][0] = raPerJ;
  changeOfBasisMatrix[2][1] = decPerJ;
  changeOfBasisMatrix[2][2] = 0;


  // TODO: Extract the dimensions for the velocity spectrum from the FITS 
  // file rather than punting and just setting it to 1 as we have done above.

  // Create a direction cosine matrix and spacing vector from the
  // change-of-basis matrix.  The direction cosines are calculated by
  // normalizing the direction vectors contained in the change-of-basis matrix
  // (i.e., dividing their values by , while also populating the spacing
  // vector:
  double spacing[3];
  vector<double> directionCosines[3];
  for (int i = 0; i < 3; ++i) {
    directionCosines[i].resize(3);
    spacing[i] = sqrt(square(changeOfBasisMatrix[i][0]) + 
                      square(changeOfBasisMatrix[i][1]) +
                      square(changeOfBasisMatrix[i][2]));
    for (int j = 0; j < 3; ++j) {
      directionCosines[i][j] = changeOfBasisMatrix[i][j] / spacing[i];
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


  // TODO: The following is a gross hack.  The problem with not doing this hack
  // is that the velocity spectrum dimension is likely to have a completely
  // different scale than the other two dimensions, and Slicer isn't clever
  // enough yet to handle such an issue in any reasonable way.  It will either
  // show the dimension as completely scrunched up, or it will be elongated
  // (and then clipped) to the point of uselessness.  So for now, we just set
  // the spacing to be the same as it is for Right Ascension:
  // spacing[2] = spacing[0];  //d

  // TODO: Remove this gross hack.  (Actually, we should put this in as a
  // command line option.) We scale the physical values by 1,000 so that they
  // are easier to use with Slicer:
  for (int i = 0; i < 3; ++i) spacing[i] *= 1000; //d

  // Set up the ITK image:
  this->SetNumberOfComponents(1);
  this->SetPixelType(SCALAR);
  this->SetComponentType(FLOAT);
  this->SetNumberOfDimensions(numOfAxes);
  for (int axis = 0; axis < numOfAxes; ++axis) {
    this->SetDimensions(axis, lengthOfAxisInPixels[axis]);
    this->SetOrigin(axis, origin[axis]);
    this->SetSpacing(axis, spacing[axis]);

    if (_cv_debugLevel) {
      cerr << "spacing=" << axis << "," << spacing[axis] << " ";
      cerr << "directions=";
      for (int j = 0; j < 3; ++j) {
	cerr << directionCosines[axis][j] << " ";
      }
      cerr << endl;
    }
    this->SetDirection(axis, directionCosines[axis]);
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
    itk::EncapsulateMetaData<string>(dict, string("FITS.") + keyName, keyValue);

    // YOU ARE HERE: You need to put the comments and units somewhere too.
  }

  // YOU ARE HERE: Check status and raise exception.
  

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
// Read(): virtual method of FITSImageIO
//   implements pure virtual method of ImageIOBase
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
  ::fits_read_pix(m_fitsFile, TFLOAT, origin, nPixels, NULL,
                  bufferAsFloats, NULL, &status);
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
// WriteImageInformation(): virtual method of FITSImageIO
//   implements pure virtual method of ImageIOBase
//-----------------------------------------------------------------------------

/*method*/ void 
FITSImageIO::WriteImageInformation()
{
  // TODO: Implement with help of cfitsio
}


//-----------------------------------------------------------------------------
// Write(): virtual method of FITSImageIO
//   implements pure virtual method of ImageIOBase
//-----------------------------------------------------------------------------

//! The Write() function is not yet implemented.

/*method*/ void 
FITSImageIO::Write(const void* const buffer) 
{
  // TODO: Maybe implement this someday.

  itkExceptionMacro("FITSImageIO::Write() not implemented.");
}


//-----------------------------------------------------------------------------
// PrintSelf(): virtual method of FITSImageIO
//   overrides partially implement method of ImageIOBase
//-----------------------------------------------------------------------------

/*method*/ void
FITSImageIO::PrintSelf(ostream& os, const Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PixelType " << m_PixelType << "\n";
}


//-----------------------------------------------------------------------------
} // END namespace itk
//-----------------------------------------------------------------------------
