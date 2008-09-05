// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  FITS Reader for ITK
// Module:   itkFITSImageIO.h
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

#ifndef __itkFITSImageIO_h
#define __itkFITSImageIO_h

#include <fstream>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <fitsio.h>
#include <itkImageIOBase.h>

/*BEGIN*/ namespace itk
{


//*****************************************************************************
//*****                                                                   *****
//*****         FITSImageIO: leaf subclass of ImageIOBase                 *****
//*****                                                                   *****
//*****************************************************************************

//! @class FITSImageIO
//!
//! @brief Read FITS file format. 
//!
//! TODO: Add here a link to documentation...
//!
//! @ingroup IOFilters

class ITK_EXPORT FITSImageIO : public ImageIOBase
{
public:

  // The number of dimensions in an image:
  enum { c_dims = 3 };
  enum CelestialCoordinateAxis { c_ra, c_vel, c_dec };
  enum FitsImageArrayAxis      { c_i,  c_j,   c_k };

  // Standard ITK typedefs:
  typedef FITSImageIO         Self;
  typedef ImageIOBase         Superclass;
  typedef SmartPointer<Self>  Pointer;

  //! Method for creation through the object factory.
  itkNewMacro(Self);

  //! Run-time type information (and related methods).
  itkTypeMacro(FITSImageIO, ImageIOBase);

  // Static setter methods:
  static void SetNullValue(double nullValue) { _cv_nullValue = nullValue; }
  static void deprecated_SetDebugLevel(int debugLevel)
                 { _cv_deprecated_debugLevel = debugLevel; }
                 
  // Static getter methods:
  static double GetNullValue() { return _cv_nullValue; }
                   

  // Virtual methods implementing pure virtual methods of ImageIOBase:
  virtual bool CanReadFile(const char* filename);
  virtual void ReadImageInformation();
  virtual void Read(void* buffer);
  virtual bool CanWriteFile(const char*);
  virtual void WriteImageInformation();
  virtual void Write(const void* buffer);

  // Functions exported for dynamic loading:
  typedef void (*NullValueSetter)(double nullValue);

protected:

  // Constructors, etc:
  FITSImageIO() {};

  // Virtual methods overriding partially implemented methods of ImageIOBase:
  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  
private:

  // Instance variables:
  fitsfile*                  m_fitsFile;

                             // Note: the memory for 'm_fitsFile' is managed
			     // by CFITSIO.  I.e., its memory is freed when
			     // when m_fitsFile is closed.


  // Private class variables:   // TODO
  static double _cv_nullValue;
  static int    _cv_deprecated_debugLevel;

  // Deactivate object copying:
  FITSImageIO(const Self&);      // Intentionally not implemented.
  void operator=(const Self&);   // Intentionally not implemented.

  // Private methods:
  std::string getFitsHeader();
};

} // END namespace itk

#endif // __itkFITSImageIO_h
