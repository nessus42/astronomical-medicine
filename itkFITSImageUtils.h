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

//-----------------------------------------------------------------------------
// Typedefs
//-----------------------------------------------------------------------------

typedef Matrix<double, 4, 4> HMatrix;
typedef Vector<double, 4> HVector;

//-----------------------------------------------------------------------------
// Constants
//-----------------------------------------------------------------------------

enum {e_i, e_j, e_k };
enum {e_x, e_y, e_z}; 
enum {e_ra, e_dec, e_vel};
enum {e_left, e_posterior, e_superior};


//-----------------------------------------------------------------------------
// Internal definitions
//-----------------------------------------------------------------------------

// BEGIN
namespace _internal {

enum { c_dims = FITSImageIO::c_dims };

using std::cout;
using std::endl;
using std::ostream;
using std::string;


//*****************************************************************************
//*****                                                                   *****
//*****             FITSImage: concrete type                              *****
//*****                                                                   *****
//*****************************************************************************

template <class ImageT>
class FITSImage
{

  typedef FITSWCSTransform<double, ImageT::ImageDimension> WCS;

public:
  typedef          FITSImage<ImageT>        Self;
  typedef typename WCS::InputPointType      IjkPoint;
  typedef typename IjkPoint::VectorType     IjkVector;
  typedef typename WCS::OutputPointType     WcsPoint;
  typedef typename WcsPoint::VectorType     WcsVector;
  typedef typename WCS::ConstPointer	    WcsTransformConstPtr;
  typedef typename WCS::Pointer		    WcsTransformPtr;


  struct Params {
    typename ImageT::Pointer           itkImage;
    bool				  wcsP;
    bool				  equiangularP;
    bool				  northUpP;
    bool				  autoscaleZAxisP;
    double 				  xAxisScale;
    double			          yAxisScale;
    double				  zAxisScale;
    double				  rotateSky;

    Params()
      : itkImage(0)
    {}
  };
    

private:

  // Instance variables:
  const Params		  _params;
  ImageT&	          _itkImage;
  WcsTransformConstPtr    _wcsTransform;
  IjkPoint                _ijkCenter;
  IjkPoint       	  _ijkSize;
  WcsPoint		  _wcsImageOrigin;
  WcsPoint                _wcsImageCenter;
  WcsVector               _unitIInWcs;
  WcsVector               _unitJInWcs;
  WcsVector               _unitIInApproximateAngularSpace;
  WcsVector               _unitJInApproximateAngularSpace;
  double		  _unitKInVelocity;         // In km/sec
  double		  _velocityAtIjkOrigin;     // In km/sec
  IjkVector            	  _ijkNorthVector;
  IjkVector		  _ijkEastVector;
  double               	  _rotationOfJFromIjkNorthVectorInDegrees;
  double		  _rotationOfIFromIjkEastVectorInDegrees;
  double               	  _raAngularScalingFactor;

  // Private static methods:
  WcsTransformConstPtr    makeWcsTransform();

  // Private methods:
  void	  initializeInstanceVars();
  void    setVelocityInstanceVars();
  double  zAxisAutoscale(bool, const HMatrix&);

  // Deactivate copy ctor and and assignment:
  FITSImage(const FITSImage&);
  void operator=(const FITSImage&);

public:

  // Constructors:
  explicit FITSImage(const Params& params);

  // Accessor methods::
  typename ImageT::Pointer
     getITKImage() { return _params.itkImage; }

  typename ImageT::ConstPointer 
     getITKImage() const
        { return typename ImageT::ConstPointer(_params.itkImage); }

  IjkPoint		  ijkSize() const        { return _ijkSize; }
  IjkPoint                ijkCenter() const      { return _ijkCenter; }
  WcsPoint	          wcsImageOrigin() const { return _wcsImageOrigin; }
  WcsPoint                wcsImageCenter() const { return _wcsImageCenter; }
  WcsVector               unitIInWcs() const     { return _unitIInWcs; }
  WcsVector               unitJInWcs() const     { return _unitJInWcs; }
  WcsTransformConstPtr    wcsTransform() const   { return _wcsTransform; }

  WcsVector     unitIInApproximateAngularSpace() const
                   { return _unitIInApproximateAngularSpace; }
  WcsVector     unitJInApproximateAngularSpace() const
                   { return _unitJInApproximateAngularSpace; }
  double	unitKInVelocity() const          { return _unitKInVelocity; }
  double	velocityAtIjkOrigin() const      { return _velocityAtIjkOrigin;}

  IjkVector     ijkNorthVector() const
                   { return _ijkNorthVector; }
  IjkVector     ijkEastVector() const
                   { return _ijkEastVector; }
  double        rotationOfJFromIjkNorthVectorInDegrees() const
                   { return _rotationOfJFromIjkNorthVectorInDegrees; }
  double        rotationOfIFromIjkEastVectorInDegrees() const
                   { return _rotationOfIFromIjkEastVectorInDegrees; }
  double        raAngularScalingFactor() const
                   { return _raAngularScalingFactor; }


  // Nonvirtual methods:
  HMatrix ijkToEquiangularMatrix() const;
  HMatrix ijkToWcsMatrix() const;
  HMatrix ijkToNorthUpMatrix() const;
  HMatrix autoscaleZMatrix() const;
};

//-----------------------------------------------------------------------------
// Functions defined in itk::fits::_internal
//-----------------------------------------------------------------------------

template <class ImageT> void
   writeImageInfo(const FITSImage<ImageT>& image, std::ostream& out);

template <class ImageT> void
   writeFitsHeader(const ImageT& image, ostream& out);

template <class PixelT> void
   scalePixelValues(Image<PixelT, _internal::c_dims>& image,
                    double multiplier);

template <class PixelT> void
   reflectPixels(Image<PixelT, _internal::c_dims>& image,
                 bool flipRAFlag, bool flipDecFlag, bool flipVFlag);

template <class ImageT> void
   rightConcatenateTransformation(ImageT& image, const HMatrix& m);

template <class ImageT> void
   leftConcatenateTransformation(ImageT& image, const HMatrix& m);

template <class ImageT> void
   setCoordinateFrameTransformation(ImageT& image, const HMatrix& m);

template <class ImageT> HMatrix
   getCoordinateFrameTransformation(const ImageT& image);


} // END namespace _internal


//-----------------------------------------------------------------------------
// Export definitions from itk::fits::_internal into itk::fits
//-----------------------------------------------------------------------------

using _internal::FITSImage;
using _internal::writeImageInfo;
using _internal::writeFitsHeader;
using _internal::scalePixelValues;
using _internal::reflectPixels;
using _internal::rightConcatenateTransformation;
using _internal::leftConcatenateTransformation;
using _internal::setCoordinateFrameTransformation;
using _internal::getCoordinateFrameTransformation;


//-----------------------------------------------------------------------------
// Functions defined in itk::fits
//-----------------------------------------------------------------------------

void           setNullValue(double nullValue);
void	       setFITSImageIODebugLevel(int debugLevel);
void           fillMatrix(HMatrix& m, const double vals[4][4]);
HMatrix        rotationMatrix(double degrees);
HMatrix	       scalingMatrix(double xScale, double yScale, double zScale);
const HMatrix& xyzToLpsMatrix();


//-----------------------------------------------------------------------------
// Inline functions
//-----------------------------------------------------------------------------

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


} } // END namespace itk::fits

#ifndef ITK_MANUAL_INSTANTIATION
#include <itkFITSImageUtils.txx>
#endif

#endif // _itkFITSImageUtils_h
