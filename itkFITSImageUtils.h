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
#include <itkImageRegionIterator.h>
#include <itkFITSWCSTransform.h>
#include <itkFITSImageIO.h>

// BEGIN
namespace itk {
namespace fits {

// Typedefs:
typedef Matrix<double, 4, 4> HMatrix;

// Constants:
enum { e_i = FITSImageIO::c_i,
       e_j = FITSImageIO::c_j,
       e_k = FITSImageIO::c_k };

enum 




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
    bool				  wcsP;  //?
    bool				  equiangularP; //?
    bool				  northUpP; //?
    bool				  eastLeftP; //?
    bool				  autoscaleZAxisP; //?
    bool				  lpsP; //?
    double 				  xAxisScale; //?
    double			          yAxisScale; //?
    double				  zAxisScale; //?
    double				  rotateSky; //?

    Params()
      : itkImage(0)
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
  HMatrix ijkToEquiangularMatrix();
  HMatrix xyzToLpsMatrix();

  // Deactivate copy ctor and and assignment:
  FITSImage(const FITSImage&);
  void operator=(const FITSImage&);

public:

  // Constructors:
  explicit FITSImage(const Params& params);

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


//-----------------------------------------------------------------------------
// Functions defined in itk::fits
//-----------------------------------------------------------------------------


void     setNullValue(double nullValue);
void     fillMatrix(HMatrix& m, const double vals[4][4]);
HMatrix  rotationMatrix(double degrees);
HMatrix	 scalingMatrix(double xScale, double yScale, double zScale);

template <class ImageT> void
   writeImageInfo(const FITSImage<ImageT>& image, ostream& out);

template <class PixelT> void
   scalePixelValues(Image<PixelT, c_dims>& image, double multiplier);

template <class PixelT> void
   reflectPixels(Image<PixelT, c_dims>& image,
                 bool flipRAFlag, bool flipDecFlag, bool flipVFlag);

template <class ImageT> void
   rightConcatenateTransformation(ImageT& image, const HMatrix& m);

template <class ImageT> void
   leftConcatenateTransformation(ImageT& image, const HMatrix& m);

template <class ImageT> void
   setCoordinateFrameTransformation(ImageT& image, const HMatrix& m);


} } // END namespaces

#ifndef ITK_MANUAL_INSTANTIATION
#include <itkFITSImageUtils.txx>
#endif

#endif // _itkFITSImageUtils_h
