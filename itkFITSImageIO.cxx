// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  FITS Reader for ITK
// Module:   itkFITSImageIO.cxx
// Package:  FITS IO
// Author:   Douglas Alan <douglas_alan AT harvard.edu>
//                        <doug AT alum.mit.edu>
//           Initiative in Innovative Computing at Harvard University
//
// Copyright (c) 2006-2007 Harvard University
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details.
//=============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itksys/RegularExpression.hxx>

#include <itkFITSImageIO.h>
#include <grparser.h> // for FITS NGP_MAX_ARRAY_DIM

#include <da_sugar.h>

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
// checkExtension(): local proc
//-----------------------------------------------------------------------------

// Returns true iff \a filepath ends with an extension that indicates
// that it is a FITS file.

local proc bool
checkExtension(const string& filepath)
{
  static itksys::RegularExpression fitsRE =
    "^.+\\.(fits?|imh)(.(gz|Z))?"
    "(\\([^()]+\\))?((\\[[^][]+\\])*|\\+[0-9]+)$";

  static itksys::RegularExpression rawRE = 
    "^.+\\.dat\\[[^][]+\\]$";

  return fitsRE.find(filepath) or rawRE.find(filepath);
}


//-----------------------------------------------------------------------------
// getAllFitsErrorMessages(): local proc
//-----------------------------------------------------------------------------

local proc string
getAllFitsErrorMessages(const int status)
{
  string retval;
  char fitsError[FLEN_ERRMSG];
  retval = "FITSIO status = ";
  stringstream ss;
  ss << status;
  retval += ss.str();
  retval += ": ";
  fits_get_errstatus(status, fitsError);
  retval += fitsError;
  return retval;
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

double FITSImageIO::_cv_nullValue = 0; // The default of 0 causes NaN's
                                       // to be left as NaN's, rather than
                                       // converted to 0, as one would
                                       // naively expect.


//=============================================================================
//                    Private methods                                   
//=============================================================================

//-----------------------------------------------------------------------------
// getFitsHeader(): private method of FITSImageIO
//-----------------------------------------------------------------------------

private_method string
FITSImageIO::getFitsHeader()
{
  int status = 0;
  const bool noComments = true;
  char* headersToExclude[0] = { };
  const int nHeadersToExclude = 0;
  char* retval1;
  int nKeysDummy;
  fits_hdr2str(m_fitsFile, !noComments, headersToExclude, nHeadersToExclude,
	       &retval1, &nKeysDummy, &status);
  da::Freer freer (retval1);
  if (status) {
    itkExceptionMacro("FITSImageIO could not get header from Primary Array of"
                      "FITS file: \""
                      << this->GetFileName() << "\": "
                      << getAllFitsErrorMessages(status) << '.');
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

method bool
FITSImageIO::CanReadFile(const char* const filepath) 
{ 
  if (!filepath or *filepath == 0) {
    itkDebugMacro(<< "No filename specified.");
    return false;
  } else if (checkExtension(filepath)) {
    return true;
  } else {
    itkDebugMacro(<<"The filename extension is not recognized");
    return false;
  }
}


//-----------------------------------------------------------------------------
// CanWriteFile(): inherited virtual method
//-----------------------------------------------------------------------------

method bool
FITSImageIO::CanWriteFile(const char* const name)
{
  // Return false because we currently don't implement writing FITS files:
  return false;
}


//-----------------------------------------------------------------------------
// ReadImageInformation(): inherited virtual method
//-----------------------------------------------------------------------------

method void
FITSImageIO::ReadImageInformation()
{
  // CFITSIO writes to this to indicate errors:
  int status = 0;   

  // Open the FITS file:
  ::fits_open_file(&m_fitsFile, this->GetFileName(), READONLY, &status);
  if (status) {
    itkExceptionMacro("FITSImageIO could not open FITS file: \""
                      << this->GetFileName() << "\" for reading: "
                      << getAllFitsErrorMessages(status) << '.');
  }

  // Get the dimensions and type of the FITS Primary Array:
  int numOfAxes;
  long lengthsOfAxesInPixels[NGP_MAX_ARRAY_DIM];
  int bitsPerPixel;
  fits_get_img_param(m_fitsFile, NGP_MAX_ARRAY_DIM, &bitsPerPixel,
		     &numOfAxes, lengthsOfAxesInPixels, &status);
  if (status) {
    itkExceptionMacro("FITSImageIO could not read Primary Array parameters "
                      "from FITS file \""
                      << this->GetFileName() << "\":"
                      << getAllFitsErrorMessages(status) << '.');
  }

  debugPrint("NAXIS=" << numOfAxes);
  if (da::getDebugLevel()) {
    cerr << "DIMS=";
    for (long i = 0; i < numOfAxes; ++i) {
      cerr << " " << lengthsOfAxesInPixels[i];
    }
    cerr << endl;
  }
  debugPrint("BPP= " << bitsPerPixel);

  if (numOfAxes != c_dims) {
    itkExceptionMacro("FITSImageIO cannot handle FITS Primary Arrays that are"
                      " anything other than three dimensional.");
    
    // TODO: Remove this restriction?
  }

  // BEGIN getting RA and Dec of border pixels.

  // TODO: Check to make sure that wcsinit copies the string handed to it by
  // fitsHeader.c_str(), as I believe that the string returned by
  // fitsHeader.c_str() will go stale as soon as fitsHeader goes out of scope.

  string fitsHeader = getFitsHeader();

  // TODO: The code below for getting the velocity information is fragile, and
  // depends on quite a few assumptions about the format of the FITS file.  It
  // would be quite nice of this could be made more robust, but it's not clear
  // at all how to do so.  A more fully-featured WCS library might be able to
  // navigate this minefield.  But then again, maybe it can't.

  // Get the coordinate frame information regarding the velocity axis from the
  // FITS file:
  double velocityAtIndexOrigin;
  double velocityDelta;
  {
    double referenceVelocity;
    double referenceVelocityIndex;
    char dummyComment[81];
    
    fits_read_key(m_fitsFile, TDOUBLE, "CRVAL3", &referenceVelocity,
		  dummyComment, &status);
    fits_read_key(m_fitsFile, TDOUBLE, "CRPIX3", &referenceVelocityIndex,
		  dummyComment, &status);
    fits_read_key(m_fitsFile, TDOUBLE, "CDELT3", &velocityDelta,
		  dummyComment, &status);
//     if (status) {
//       itkExceptionMacro("FITSImageIO could not read velocity WCS info from"
// 			" FITS file \""
// 			<< this->GetFileName() << "\": "
// 			<< ::getAllFitsErrorMessages(status) << '.');
//      }

    if (!status) {

      // Convert from m/s to km/s:
      referenceVelocity /= 1000; 
      velocityDelta /= 1000; 
      velocityAtIndexOrigin = 
        referenceVelocity - velocityDelta * ( 1 - referenceVelocityIndex);
      debugPrint(
                 "velocityAtIndexOrigin=" << velocityAtIndexOrigin << "\n"
                 "indexOfZeroVelocity="
                 << referenceVelocityIndex - (referenceVelocity / velocityDelta)
      );
    }
  }

  // Set up the ITK image:
  { 
    // TODO: Note that we set the "component type" to float immediately below,
    // and yet in fits2itk, we allow the creation of short and unsigned short
    // images.  This means that when an image of shorts is being read in, it
    // gets converted to floats here, and then gets converted back to shorts
    // later (by ITK, not by us).  It might be better to do this in a less
    // round-about fashion.  Note: If we do, we also have to change the call to
    // cfitsio's fits_read_pix() to allow types other than floats.

    this->SetNumberOfComponents(1);
    this->SetPixelType(SCALAR);
    this->SetComponentType(FLOAT);
    this->SetNumberOfDimensions(numOfAxes);
    for (int indexAxis = 0; indexAxis < numOfAxes; ++indexAxis) {
      this->SetDimensions(indexAxis, lengthsOfAxesInPixels[indexAxis]);
    }
  }

  // Add FITS information to the ITK MetaDataDictionary:
  {
    // Put the FITS Primary Array Header into the MDD as one big string:
    MetaDataDictionary& dict = this->GetMetaDataDictionary();
    EncapsulateMetaData(dict, "FITS Header", fitsHeader);

    // Also break up the aforementioned FITS header into individual entries
    // and add each entry into the MDD as a separate entry:
    int nKeys = 0;
    int dummy;
    fits_get_hdrspace(m_fitsFile, &nKeys, &dummy, &status);
    for (int keyIndex = 1; keyIndex <= nKeys; ++keyIndex) {
      char keyName[FLEN_KEYWORD];
      char keyValue[FLEN_VALUE];
      char keyComment[FLEN_COMMENT];
      fits_read_keyn(m_fitsFile, keyIndex, keyName, keyValue, keyComment,
		     &status);
      EncapsulateMetaData(dict, string("FITS.") + keyName, string(keyValue));

      // TODO: You need to put the comments and units somewhere too.
    }

    // Put the velocity of the first slice and the last slice into the MDD so
    // that I can extract it easily in fits2itk.c:
    {
      EncapsulateMetaData(dict, "fits2itk.velocityAtIndexOrigin",
			  velocityAtIndexOrigin);
      EncapsulateMetaData(dict, "fits2itk.velocityDelta", velocityDelta);
    }
    // TODO: Making separate MDD entries for the above velocites is a bit
    // gross, and perhaps it would just be better to grab them out of the ITK
    // Image spacing and origin information.  The problem with doing that,
    // however, is that if velocity autoscaling is set, then the velocity
    // information will not be accurate.  I think, however, that in the future
    // taking that approach may be okay, as I think that autoscaling will not
    // be done here, but rather more directly in fits2itk.c (or one of its
    // helper files), and it can do the autoscaling after it has extracted the
    // real velocity information out of the Image object.

    // TODO: Check status and raise exception.
  }
} // end namespace fitsio


//-----------------------------------------------------------------------------
// Read(): inherited virtual method
//-----------------------------------------------------------------------------

method void
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
  int foundNull = false;
  fits_read_pix(m_fitsFile, TFLOAT, origin, nPixels, (void*) nullValuePtr,
		bufferAsFloats, &foundNull, &status);
  if (nullValue != 0) {
    debugPrint("Any null values? " << foundNull);
    debugPrint("Null value will be fetched as: " << nullValue);
  }
  if (status) {
    itkExceptionMacro("FITSImageIO::Read() could not read the primary data "
                      "array from FITS file \""
                      << this->GetFileName() << "\": "
                      << getAllFitsErrorMessages(status) << ".");
  }
                
  // Close the FITS file:
  fits_close_file(m_fitsFile, &status);
  if (status) {
    itkExceptionMacro("FITSImageIO::Read() could not close FITS file \""
                      << this->GetFileName()
                      << "\" after reading the primary data array: "
                      << getAllFitsErrorMessages(status) << ".");
  }
}


//-----------------------------------------------------------------------------
// WriteImageInformation(): inherited virtual method
//-----------------------------------------------------------------------------

method void 
FITSImageIO::WriteImageInformation()
{
  itkExceptionMacro("FITSImageIO::WriteImageInformation() not implemented "
		    "yet.");
}


//-----------------------------------------------------------------------------
// Write(): inherited virtual method
//-----------------------------------------------------------------------------

method void 
FITSImageIO::Write(const void* const buffer) 
{
  itkExceptionMacro("FITSImageIO::Write() not implemented.");
}


//-----------------------------------------------------------------------------
// PrintSelf(): inherited virtual method
//-----------------------------------------------------------------------------

method void
FITSImageIO::PrintSelf(ostream& os, const Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PixelType " << m_PixelType << "\n";
}

//-----------------------------------------------------------------------------
// itkFITSImageIO_setNullValue(): function exported for dynamic loading
//-----------------------------------------------------------------------------

proc extern "C" void 
itkFITSImageIO_setNullValue(double nullValue)
{ 
  FITSImageIO::SetNullValue(nullValue);
}

//-----------------------------------------------------------------------------
// itkFITSImageIO_setNullValue(): function exported for dynamic loading
//-----------------------------------------------------------------------------

proc extern "C" void 
itkFITSImageIO_setDebugLevel(int debugLevel)
{ 
  da::setDebugLevel(debugLevel);
}


} // namespace itk
