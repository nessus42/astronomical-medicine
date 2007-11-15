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

#include <itkFITSWCSTransform.h>

// BEGIN
namespace itk
{

// URGENT TODO: Fix this attrocity!  Pass the information in the
// meta-data dictionary, or something, instead.  This won't work right if the
// program uses more than one FITS Image!

extern void* g_theFITSWCSTransform;

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

  // Standard class typedefs:
  typedef FITSImageIO         Self;
  typedef ImageIOBase         Superclass;
  typedef SmartPointer<Self>  Pointer;
  
  //! Method for creation through the object factory.
  itkNewMacro(Self);

  //! Run-time type information (and related methods).
  itkTypeMacro(FITSImageIO, ImageIOBase);

  // Static setter methods:
  static void SetSuppressWCS(bool flag)
                 { _cv_suppressWCS = flag; }
  static void SetDebugLevel(int debugLevel)
                 { _cv_debugLevel = debugLevel; }
  static void SetRIPOrientation(bool flag)
                 { _cv_RIPOrientationFlag = flag; }
  static void SetRotateSky(double degrees)
                 { _cv_rotateSky = degrees; }
  static void SetScaleVelocity(double scalingFactor)
                 { _cv_scaleVelocity = scalingFactor; }
  static void SetAutoScaleVelocityAxis(bool flag)
                 { _cv_autoScaleVelocityAxis = flag; }
  static void SetScaleVoxelValues(double scalingFactor)
                 { _cv_scaleVoxelValues = scalingFactor; }
  static void SetScaleAllAxes(double scalingFactor)
                 { _cv_scaleAllAxes = scalingFactor; }
  static void SetScaleDec(double scalingFactor)
                 { _cv_scaleDec = scalingFactor; }
  static void SetScaleRA(double scalingFactor)
                 { _cv_scaleRA = scalingFactor; }
  static void SetNullValue(double nullValue)
                 { _cv_nullValue = nullValue; }
  static void SetSuppressMetaDataDictionary(bool flag)
                 { _cv_suppressMetaDataDictionary = flag; }
  static void SetVerbose(bool flag)
                 { _cv_verbose = flag; }
					
  // Static getter methods:
  static bool   GetSuppressWCS()
                  { return _cv_suppressWCS; }
  static int    GetDebugLevel()
                  { return _cv_debugLevel; }
  static bool   GetRIPOrientation()
                  { return _cv_RIPOrientationFlag; }
  static double GetRotateSky()
                  { return _cv_rotateSky; }
  static double GetScaleVelocity()
                  { return _cv_scaleVelocity; }
  static bool   GetAutoScaleVelocityAxis()
                  { return _cv_autoScaleVelocityAxis; }
  static double GetScaleVoxelValues()
                  { return _cv_scaleVoxelValues; }
  static double GetScaleAllAxes()
                  { return _cv_scaleAllAxes; }
  static double GetScaleDec()
                  { return _cv_scaleDec; }
  static double GetScaleRA()
                  { return _cv_scaleRA; }
  static double GetNullValue()
                  { return _cv_nullValue; }
  static bool   GetSuppressMetaDataDictionary()
                  { return _cv_suppressMetaDataDictionary; }
  static bool   GetVerbose()
                  { return _cv_verbose; }

  // Virtual methods implementing pure virtual methods of ImageIOBase:
  virtual bool CanReadFile(const char* filename);
  virtual void ReadImageInformation();
  virtual void Read(void* buffer);
  virtual bool CanWriteFile(const char*);
  virtual void WriteImageInformation();
  virtual void Write(const void* buffer);

protected:

  // Constructors, etc:
  FITSImageIO() {};

  // Virtual methods overriding partially implemented methods of ImageIOBase:
  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  
private:

  // Instance variables:
  fitsfile*                                      m_fitsFile;
  FITSWCSTransform<double, c_dims>::Pointer      m_WCSTransform;

				// Note: the memory for 'm_fitsFile' is managed
				// by CFITSIO.  I.e., its memory is freed when
				// when m_fitsFile is closed.

  // Private class variables:   // TODO
  static int    _cv_debugLevel;
  static bool   _cv_suppressWCS;
  static double _cv_nullValue;
  static bool   _cv_RIPOrientationFlag;
  static double _cv_rotateSky;
  static double _cv_rollRA;
  static double _cv_rollDec;
  static double _cv_scaleVoxelValues;
  static double _cv_scaleRA;
  static double _cv_scaleDec;
  static double _cv_scaleVelocity;
  static bool   _cv_autoScaleVelocityAxis;
  static double _cv_scaleAllAxes;
  static bool   _cv_suppressMetaDataDictionary;
  static bool	_cv_verbose;

  // Deactivate object copying:
  FITSImageIO(const Self&);      // Intentionally not implemented.
  void operator=(const Self&);   // Intentionally not implemented.

  // Private methods:
  std::string getFitsHeader();
};

} // END namespace itk

#endif // __itkFITSImageIO_h
