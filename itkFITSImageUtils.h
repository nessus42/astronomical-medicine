// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Module: itkFITSImageUtils
//
// Description: Utilities for manipulating itk::FITSImage objects.
//
// Author:
//      Douglas Alan <douglas_alan AT harvard.edu>
//                   <doug AT alum.mit.edu>
//      Initiative in Innovative Computing at Harvard University
//
// Copyright (c) 2006-2007 Harvard University
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details.
//=============================================================================
#ifndef _itkFITSImageUtils_h
#define _itkFITSImageUtils_h

#include <itkMatrix.h>
#include <itkFITSWCSTransform.h>
#include <itkFITSImageIO.h>

// BEGIN
namespace itk {
namespace fits {
namespace _internal {

// Internal typedefs:
typedef Matrix<double, 3, 3> Matrix;


// Global functions:
void setNullValue(double nullValue);


//*****************************************************************************
//*****                                                                   *****
//*****             FITSImage: concrete type                              *****
//*****                                                                   *****
//*****************************************************************************

template <class ImageType>
class FITSImage
{

  typedef FITSWCSTransform<double, ImageType::ImageDimension> WCS;

public:
  typedef          FITSImage<ImageType>     Self;
  typedef typename WCS::InputPointType      IjkPoint;
  typedef typename IjkPoint::VectorType     IjkVector;
  typedef typename WCS::OutputPointType     WcsPoint;
  typedef typename WcsPoint::VectorType     WcsVector;
  typedef typename WCS::ConstPointer	    WcsTransformConstPtr;
  typedef typename WCS::Pointer		    WcsTransformPtr;


  struct Params {
    typename ImageType::Pointer           itkImage;
    double                                angularUnitsInMicroDegrees;

    Params()
      : itkImage(0),
	angularUnitsInMicroDegrees(1)
    {}
  };
    

private:

  // Instance variables:
  const Params		  _params;
  ImageType&	          _itkImage;
  WcsTransformConstPtr    _wcsTransform;
  IjkPoint                _ijkCenter;
  WcsPoint                _wcsCenter;
  WcsVector               _unitIInWcs;
  WcsVector               _unitJInWcs;
  WcsVector               _unitIInApproximateAngularSpace;
  WcsVector               _unitJInApproximateAngularSpace;
  IjkVector            	  _ijkNorthVector;
  double               	  _rotationOfJFromIjkNorthVectorInDegrees;
  double               	  _raAngularScalingFactor;

  // Private static methods:
  WcsTransformConstPtr    makeWcsTransform();

  // Private methods:
  void	 initializeInstanceVars();
  Matrix ijkToNorthOrientedEquiangularMatrix();
  Matrix raDecVToLpsMatrix();

  // Deactivate copy ctor and and assignment:
  FITSImage(const FITSImage&);
  void operator=(const FITSImage&);

public:

  // Constructors:
  explicit FITSImage(const Params& params);

  // YOU ARE HERE: You need to modify the constructor to accept the params.

  // Accessor methods::
  typename ImageType::Pointer
     getITKImage() { return _params.itkImage; }

  typename ImageType::ConstPointer 
     getITKImage() const
        { return typename ImageType::ConstPointer(_params.itkImage); }

  IjkPoint                ijkCenter() const    { return _ijkCenter; }
  WcsPoint                wcsCenter() const    { return _wcsCenter; }
  WcsVector               unitIInWcs() const   { return _unitIInWcs; }
  WcsVector               unitJInWcs() const   { return _unitJInWcs; }
  WcsTransformConstPtr    wcsTransform() const { return _wcsTransform; }

  WcsVector     unitIInApproximateAngularSpace() const
                   { return _unitIInApproximateAngularSpace; }
  WcsVector     unitJInApproximateAngularSpace() const
                   { return _unitJInApproximateAngularSpace; }
  IjkVector     ijkNorthVector() const
                   { return _ijkNorthVector; }
  double        rotationOfJFromIjkNorthVectorInDegrees() const
                   { return _rotationOfJFromIjkNorthVectorInDegrees; }
  double        raAngularScalingFactor() const
                   { return _raAngularScalingFactor; }

};

// Internal functions:
void     fillMatrix(Matrix& m, const double vals[3][3]);
Matrix   rotationMatrix(double degrees);
Matrix	 scalingMatrix(double xScale, double yScale, double zScale);

// Inline internal functions:

inline double
degreesToRadians(double degrees)
{
  return degrees * M_PI/180;
}

inline double
radiansToDegrees(double radians)
{
  return radians * 180/M_PI;
}

inline double
square(double x)
{
  return x * x;
}

inline double
cartesianLength(double x, double y)
{
  return sqrt(square(x) + square(y));
}

inline bool
isOdd(size_t num)
{
  return num & 1;
}

} // END namespace _internal

//-----------------------------------------------------------------------------
// Export symbols from _internal into itk::fits
//-----------------------------------------------------------------------------

using _internal::setNullValue;
using _internal::FITSImage;

} } // END namespaces

#ifndef ITK_MANUAL_INSTANTIATION
#include <itkFITSImageUtils.txx>
#endif

#endif // _itkFITSImageUtils_h
