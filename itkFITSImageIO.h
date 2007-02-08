/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFITSImageIO.h,v $
  Language:  C++
  Date:      $Date: 2005/09/28 15:41:54 $
  Version:   $1.0$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFITSImageIO_h
#define __itkFITSImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <fstream>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <fitsio.h>
#include <itkImageIOBase.h>

namespace itk
{
  
//! \class FITSImageIO
//!
//! \brief Read FITS file format. 
//!
//! TODO: Add here a link to documentation...
//!
//! \ingroup IOFilters

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

  // Class procedures:
  static void SetDebugLevel(int debugLevel)
                 { _cv_debugLevel = debugLevel; }
  static void SetRotateSky(double degrees)
                 { _cv_rotateSky = degrees; }
  static void SetScaleVelocityAxis(double scalingFactor)
                 { _cv_scaleVelocityAxis = scalingFactor; }
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
    

  // Virtual methods implementing pure virtual methods of ImageIOBase:
  virtual bool CanReadFile(const char* filename);
  virtual void ReadImageInformation();
  virtual void Read(void* buffer);
  virtual bool CanWriteFile(const char*);
  virtual void WriteImageInformation();
  virtual void Write(const void* buffer);

protected:

  FITSImageIO() {}; // default ctor

  // Virtual methods overriding partially implemented methods of ImageIOBase:
  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  
private:

  // Instance variables:
  fitsfile* m_fitsFile;

  // Private class variables:   // TODO
  static int    _cv_debugLevel;
  static double _cv_nullValue;
  static double _cv_rotateSky;
  static double _cv_rotateDecIntoVelocityAxis;
  static double _cv_rotateRAIntoVelocityAxis;
  static double _cv_scaleVoxelValues;
  static double _cv_scaleRA;
  static double _cv_scaleDec;
  static double _cv_scaleVelocityAxis;
  static bool   _cv_autoScaleVelocityAxis;
  static double _cv_scaleAllAxes;

  // Deactivate object copying:
  FITSImageIO(const Self&) { assert(false); }
  void operator=(const Self&) { assert(false); }

  // Private methods:
  std::string getFitsHeader();
};

} // end namespace itk

#endif // __itkFITSImageIO_h
