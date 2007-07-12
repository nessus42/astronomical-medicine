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

// TODO: Delete the following pragma.  I'm pretty sure that it's not
// needed, as it's included in itkWin32Header.h:

// #ifdef _MSC_VER
// #pragma warning ( disable : 4786 )
// #endif

#include <fstream>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <fitsio.h>
#include <itkImageIOBase.h>

#include <itkFITSWCSTransform.h>


namespace itk
{

//   namespace fitsio {
//     void _calcCoordinateFrame(
// 	 const std::string& fitsHeader,
// 	 const long lengthOfAxisInPixels[],
// 	 double origin[],
// 	 double spacing[],
// 	 std::vector<double> directionCosines[],
// 	 FITSWCSTransform<double, 3>::Pointer& transform);
//   }


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

  // Standard class typedefs:
  typedef FITSImageIO         Self;
  typedef ImageIOBase         Superclass;
  typedef SmartPointer<Self>  Pointer;
  
  //! Method for creation through the object factory.
  itkNewMacro(Self);

  //! Run-time type information (and related methods).
  itkTypeMacro(FITSImageIO, ImageIOBase);

  // Class setter procedures:
  static void SetSuppressWCS(bool flag) { _cv_suppressWCS = flag; }
  static void SetDebugLevel(int debugLevel)
                 { _cv_debugLevel = debugLevel; }
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
  static void SetScaleRA(double scalingFactor)
                 { _cv_scaleRA = scalingFactor; }
  static void SetNullValue(double nullValue)
                 { _cv_nullValue = nullValue; }
  static void SetSuppressMetaDataDictionary(bool flag)
                 { _cv_suppressMetaDataDictionary = flag; }
					
  // Static getter procedures:
  static bool   GetSuppressWCS() { return _cv_suppressWCS; }
  static int    GetDebugLevel()
                 { return _cv_debugLevel; }
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
  static double GetScaleRA()
                 { return _cv_scaleRA; }
  static double GetNullValue()
                 { return _cv_nullValue; }
  static bool   GetSuppressMetaDataDictionary()
                 { return _cv_suppressMetaDataDictionary; }

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
  FITSWCSTransform<double, 3>::Pointer           m_transform;

				// Note: the memory for 'm_fitsFile' is managed
				// by CFITSIO.  I.e., its memory is freed when
				// when m_fitsFile is closed.

  // Private class variables:   // TODO
  static int    _cv_debugLevel;
  static bool   _cv_suppressWCS;
  static double _cv_nullValue;
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

  // Deactivate object copying:
  FITSImageIO(const Self&);      // Intentionally not implemented.
  void operator=(const Self&);   // Intentionally not implemented.

  // Private methods:
  std::string getFitsHeader();
};

} // end namespace itk

#endif // __itkFITSImageIO_h
