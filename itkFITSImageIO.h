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
// #include <itkFITSWCSTransform.h>

// BEGIN
namespace itk
{

// URGENT TODO: Fix this attrocity!  Pass the information in the
// meta-data dictionary, or something, instead.  This won't work right if the
// program uses more than one FITS Image!

// extern void* g_theFITSWCSTransform;

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
  enum FitsImageArrayAxis { c_i, c_j, c_k };

  // Standard ITK typedefs:
  typedef FITSImageIO         Self;
  typedef ImageIOBase         Superclass;
  typedef SmartPointer<Self>  Pointer;

  // Class-specific typedefs:
  // typedef FITSWCSTransform<double, c_dims> WCSTransform;
  
  //! Method for creation through the object factory.
  itkNewMacro(Self);

  //! Run-time type information (and related methods).
  itkTypeMacro(FITSImageIO, ImageIOBase);

  // Static setter methods:
  static void SetNullValue(double nullValue)
                 { _cv_nullValue = nullValue; }
  static void deprecated_SetDebugLevel(int debugLevel)
                 { _cv_deprecated_debugLevel = debugLevel; }

//   static void deprecated_SetSuppressWCS(bool flag)
//                  { _cv_deprecated_suppressWCS = flag; }
//   static void deprecated_SetRIPOrientation(bool flag)
//                  { _cv_deprecated_RIPOrientationFlag = flag; }
//   static void deprecated_SetRotateSky(double degrees)
//                  { _cv_deprecated_rotateSky = degrees; }
//   static void deprecated_SetScaleVelocity(double scalingFactor)
//                  { _cv_deprecated_scaleVelocity = scalingFactor; }
//   static void deprecated_SetAutoScaleVelocityAxis(bool flag)
//                  { _cv_deprecated_autoScaleVelocityAxis = flag; }
//   static void deprecated_SetScaleVoxelValues(double scalingFactor)
//                  { _cv_deprecated_scaleVoxelValues = scalingFactor; }
//   static void deprecated_SetScaleAllAxes(double scalingFactor)
//                  { _cv_deprecated_scaleAllAxes = scalingFactor; }
//   static void deprecated_SetScaleDec(double scalingFactor)
//                  { _cv_deprecated_scaleDec = scalingFactor; }
//   static void deprecated_SetScaleRA(double scalingFactor)
//                  { _cv_deprecated_scaleRA = scalingFactor; }
//  static void deprecated_SetSuppressMetaDataDictionary(bool flag)
//               { _cv_deprecated_suppressMetaDataDictionary = flag; }
//   static void deprecated_SetVerbose(bool flag)
//                  { _cv_deprecated_verbose = flag; }
					
  // Static getter methods:
  static double GetNullValue()
                   { return _cv_nullValue; }

//   static bool   deprecated_GetSuppressWCS()
//                   { return _cv_deprecated_suppressWCS; }
//   static int    deprecated_GetDebugLevel()
//                   { return _cv_deprecated_debugLevel; }
//   static bool   deprecated_GetRIPOrientation()
//                   { return _cv_deprecated_RIPOrientationFlag; }
//   static double deprecated_GetRotateSky()
//                   { return _cv_deprecated_rotateSky; }
//   static double deprecated_GetScaleVelocity()
//                   { return _cv_deprecated_scaleVelocity; }
//   static bool   deprecated_GetAutoScaleVelocityAxis()
//                   { return _cv_deprecated_autoScaleVelocityAxis; }
//   static double deprecated_GetScaleVoxelValues()
//                   { return _cv_deprecated_scaleVoxelValues; }
//   static double deprecated_GetScaleAllAxes()
//                   { return _cv_deprecated_scaleAllAxes; }
//   static double deprecated_GetScaleDec()
//                   { return _cv_deprecated_scaleDec; }
//   static double deprecated_GetScaleRA()
//                   { return _cv_deprecated_scaleRA; }
//   static bool   deprecated_GetSuppressMetaDataDictionary()
//                   { return _cv_deprecated_suppressMetaDataDictionary; }
//   static bool   deprecated_GetVerbose()
//                   { return _cv_deprecated_verbose; }

  // Virtual methods implementing pure virtual methods of ImageIOBase:
  virtual bool CanReadFile(const char* filename);
  virtual void ReadImageInformation();
  virtual void Read(void* buffer);
  virtual bool CanWriteFile(const char*);
  virtual void WriteImageInformation();
  virtual void Write(const void* buffer);

  // Functions exported for dynamic loading:
  typedef void (*NullValueSetter)(double nullValue);
  // typedef void* (*WCSTransformGetter)();

protected:

  // Constructors, etc:
  FITSImageIO() {};

  // Virtual methods overriding partially implemented methods of ImageIOBase:
  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  
private:

  // Instance variables:
  fitsfile*                  m_fitsFile;
  // WCSTransform::Pointer      m_WCSTransform;

			     // Note: the memory for 'm_fitsFile' is managed
			     // by CFITSIO.  I.e., its memory is freed when
			     // when m_fitsFile is closed.


  // Private class variables:   // TODO
  static double _cv_nullValue;
  static int    _cv_deprecated_debugLevel;
//   static bool   _cv_deprecated_suppressWCS;
//   static bool   _cv_deprecated_RIPOrientationFlag;
//   static double _cv_deprecated_rotateSky;
//   static double _cv_deprecated_rollRA;
//   static double _cv_deprecated_rollDec;
//   static double _cv_deprecated_scaleVoxelValues;
//   static double _cv_deprecated_scaleRA;
//   static double _cv_deprecated_scaleDec;
//   static double _cv_deprecated_scaleVelocity;
//   static bool   _cv_deprecated_autoScaleVelocityAxis;
//   static double _cv_deprecated_scaleAllAxes;
//   static bool   _cv_deprecated_suppressMetaDataDictionary;
//   static bool	_cv_deprecated_verbose;

  // Deactivate object copying:
  FITSImageIO(const Self&);      // Intentionally not implemented.
  void operator=(const Self&);   // Intentionally not implemented.

  // Private methods:
  std::string getFitsHeader();
};

} // END namespace itk

#endif // __itkFITSImageIO_h
