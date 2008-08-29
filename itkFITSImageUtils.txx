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
#ifndef _itkFITSImageUtils_txx
#define _itkFITSImageUtils_txx

#include <assert.h>
#include <iostream>

#include <wcs.h>

#include <itkImage.h>
#include <itkMetaDataObject.h>

#include <itkFlipImageFilter.h>
#include <itkBinomialBlurImageFilter.h>
// #include <itkDerivativeImageFilter.h>
// #include <itkMeanImageFilter.h>
// #include <itkBinaryMedianImageFilter.h>
// #include <itkGradientAnisotropicDiffusionImageFilter.h>

#include <itkFITSImageIO.h>
#include <itkFITSImageUtils.h>

// BEGIN
namespace itk {
namespace fits { 

// BEGIN
namespace _internal {

//-----------------------------------------------------------------------------
// Internal definitions
//-----------------------------------------------------------------------------

enum { c_dims = itk::FITSImageIO::c_dims };

using itk::Image;

using std::cout;
using std::endl;
using std::ostream;

using std::string;

//-----------------------------------------------------------------------------
// Constructor of FITSImage
//-----------------------------------------------------------------------------

/*ctor*/
template <class ImageT>
FITSImage<ImageT>::FITSImage(const typename Self::Params& params)
  : _params(params), _itkImage(*params.itkImage)
{
  assert(&_itkImage);
  initializeInstanceVars();

  // Set the ITK Image's coordinate transformation parameters:
  {
    const double angularUnitsScaling = 
      1e6 / _params.angularUnitsInMicroDegrees;
    setCoordinateFrameTransformation(
       _itkImage,
       this->raDecVToLpsMatrix()
       * scalingMatrix(angularUnitsScaling * params.raScale,
                       angularUnitsScaling,
                       1e4)
       * this->ijkToNorthOrientedEquiangularMatrix());
  }

  // YOU ARE HERE: You want to make the ijkToNorthOrientedEquiangularMatrix be
  // null if the suppress WCS flag is set to be true.  Also, you probably want to
  // suppress reading the WCS information to begin with if this is the case, as
  // Jens has some FITS images that have no WCS data.

  // YOU ARE HERE: You need to add back in autoscaling for Vel, rather than
  // just multiplying it by 1e4.
}


//-----------------------------------------------------------------------------
// initializeInstanceVars(): private method of FITSImage template class
//-----------------------------------------------------------------------------

template <class ImageT>
void 
FITSImage<ImageT>::initializeInstanceVars()
{
  typename ImageT::RegionType
    allOfImage = _itkImage.GetLargestPossibleRegion();
  typename ImageT::SizeType imageSize = allOfImage.GetSize();
  typename ImageT::IndexType imageOrigin = allOfImage.GetIndex();

  enum { c_i = FITSImageIO::c_i,
 	 c_j = FITSImageIO::c_j,
 	 c_k = FITSImageIO::c_k };

  _ijkCenter[c_i] = imageOrigin[c_i] + imageSize[c_i]/2.0 - 0.5;
  _ijkCenter[c_j] = imageOrigin[c_j] + imageSize[c_j]/2.0 - 0.5;
  _ijkCenter[c_k] = imageOrigin[c_k] + imageSize[c_k]/2.0 - 0.5;

  WcsTransformConstPtr wcs = makeWcsTransform(); 
  _wcsTransform = wcs;
  _wcsCenter = wcs->TransformPoint(_ijkCenter);
  
  // TODO: It's confusing that in some situations V is 1 and dec is 2, and
  // in others, dec is 1 and V is 2.  We need a better way to denote this.

  enum {ra, dec, v};

  // Calculate lengths of i and j vectors in RA/Dec space:
  { 
    IjkPoint ijkLeftHalfAPixel  = _ijkCenter;
    IjkPoint ijkRightHalfAPixel = _ijkCenter;
    IjkPoint ijkDownHalfAPixel  = _ijkCenter;
    IjkPoint ijkUpHalfAPixel    = _ijkCenter;
    ijkLeftHalfAPixel[c_i]  -= .5;
    ijkRightHalfAPixel[c_i] += .5;
    ijkDownHalfAPixel[c_j]  -= .5;
    ijkUpHalfAPixel[c_j]    += .5;
      
    WcsPoint wcsLeftHalfAPixel  = wcs->TransformPoint(ijkLeftHalfAPixel);
    WcsPoint wcsRightHalfAPixel = wcs->TransformPoint(ijkRightHalfAPixel);
    WcsPoint wcsDownHalfAPixel  = wcs->TransformPoint(ijkDownHalfAPixel);
    WcsPoint wcsUpHalfAPixel    = wcs->TransformPoint(ijkUpHalfAPixel);
      
    _unitIInWcs = wcsRightHalfAPixel - wcsLeftHalfAPixel;
    _unitJInWcs = wcsUpHalfAPixel - wcsDownHalfAPixel;
  }

  _raAngularScalingFactor = cos(degreesToRadians(_wcsCenter[dec]));
  _unitIInApproximateAngularSpace = _unitIInWcs;
  _unitJInApproximateAngularSpace = _unitJInWcs;
  _unitIInApproximateAngularSpace[ra] *= _raAngularScalingFactor;
  _unitJInApproximateAngularSpace[ra] *= _raAngularScalingFactor;

  // Calculate the image's rotation by finding a northward-oriented vector in
  // WCS space, and then transforming it into IJK space.  We can then use trig
  // to determine the amount of rotation of the image from north in IJK space:
  {
    WcsTransformPtr inverseWcs = WCS::New();
    wcs->GetInverse(inverseWcs);

    // We calculate wcsNorthVector, just for the purpose of getting a Dec
    // increment that is about the size of a pixel or so.  We could, of course,
    // accomplish much the same thing just by selecting an arbitrary small
    // increment northward, but then it might be too small or it might be too
    // big for the image:
    WcsVector wcsNorthVector =
      fabs(_unitJInWcs[dec]) > fabs(_unitIInWcs[dec])
      ? _unitJInWcs
      : _unitIInWcs;
    wcsNorthVector[ra] = 0;
    wcsNorthVector[dec] = fabs(wcsNorthVector[dec]);

    // We now add the northward increment vector to our center point expressed
    // in RA/Dec coordinates.  This gives us a slightly northward point in WCS
    // coordinates, and we then use an inverse WCS transform to get the IJK
    // pixel position (with floating point IJK index values) of this
    // fastidiously calculated northward point:
    WcsPoint wcsPointNorthOfCenter = _wcsCenter + wcsNorthVector;
    IjkPoint ijkPointNorthOfCenter =
      inverseWcs->TransformPoint(wcsPointNorthOfCenter);

    // We then subtract the center point in IJK coordinates from the northward
    // point to get a northward pointing vector in IJK coordinates:
    _ijkNorthVector = ijkPointNorthOfCenter - _ijkCenter;
    _ijkNorthVector[v] = 0;

    // And finally, with the northward pointing vector in IJK coordinates, we
    // can determine how much the image is rotated from north:
    _rotationOfJFromIjkNorthVectorInDegrees =
      radiansToDegrees(atan2(-1 * _ijkNorthVector[c_i], _ijkNorthVector[c_j]));
  }
}


//-----------------------------------------------------------------------------
// makeWcsTransform(): private static method of FITSImage template class
//-----------------------------------------------------------------------------

/*proc*/
template <class ImageT>
typename FITSImage<ImageT>::WcsTransformConstPtr
FITSImage<ImageT>::makeWcsTransform()
{
  using itk::MetaDataDictionary;
  using itk::MetaDataObject;
  const MetaDataDictionary& mdd = _itkImage.GetMetaDataDictionary();

  // We cast away const here conly because ExposeMetaData() is unfortunately
  // not const correct:
  string fitsHeader;
  ExposeMetaData(const_cast<MetaDataDictionary&>(mdd),
		 "FITS Header",
		 fitsHeader);
  assert(!fitsHeader.empty());

  // TODO: Add more error checking here in case 'wcs' ends up in some sort of
  // erroneous state.
  const ConstRcMallocPointer<WorldCoor> wcs = wcsinit(fitsHeader.c_str());
  WcsTransformPtr retval = WCS::New();
  retval->SetWCS(wcs);
  WcsTransformConstPtr constRetval (retval);
  return constRetval;
}


//-----------------------------------------------------------------------------
// writeImageInfo(): template function
//-----------------------------------------------------------------------------

/*proc*/
template <class ImageT>
void
writeImageInfo(const FITSImage<ImageT>& image, ostream& out)
{
  const ImageT& itkImage = *image.getITKImage();
  out << "Image center, in IJK space with (0,0,0) index origin: "
      << image.ijkCenter() << "\n"
      << "Image center, in RA/Dec: " << image.wcsCenter() << "\n"
      << "I vector, in RA/Dec: "
      << image.unitIInWcs() << "\n"
      << "J vector, in RA/Dec: "
      << image.unitJInWcs() << "\n"
      << "I vector, in approximate angular space: "
      << image.unitIInApproximateAngularSpace() << "\n"
      << "J vector, in approximate angular space: "
      << image.unitJInApproximateAngularSpace() << "\n"
      << "|I|, in approximate angular space: "
      << image.unitIInApproximateAngularSpace().GetNorm() << "\n"
      << "|J|, in approximate angular space: "
      << image.unitJInApproximateAngularSpace().GetNorm() << "\n"
      << "North vector in IJK space: " << image.ijkNorthVector() << "\n"
      << "Rotation of J from North:  "
      << image.rotationOfJFromIjkNorthVectorInDegrees() << "\n"
      << "Direction cosines:\n"
      << itkImage.GetDirection()
      << "Image spacing: " << itkImage.GetSpacing() << "\n"
      << "Image origin: " << itkImage.GetOrigin() << "\n"
    ;
}


//-----------------------------------------------------------------------------
// applyFlipImageFilter(): template function
//-----------------------------------------------------------------------------

/*proc*/
template <class PixelT>
typename Image<PixelT, c_dims>::Pointer
applyFlipImageFilter(const typename Image<PixelT, c_dims>::Pointer& image)
{
  typedef Image<PixelT, c_dims> ImageT;
  typedef itk::FlipImageFilter<ImageT> FilterType;
  typedef typename FilterType::FlipAxesArrayType FlipAxesArrayType;
  static typename FilterType::Pointer filter = FilterType::New();
  FlipAxesArrayType flipArray;
  flipArray[0] = 0;
  flipArray[1] = 1;
  flipArray[2] = 0;
  filter->SetFlipAxes(flipArray);
  filter->SetInput(image);
  // filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// applyBinomialBlur(): template function
//-----------------------------------------------------------------------------

/*proc*/ 
template <class PixelT>
typename Image<PixelT, c_dims>::Pointer
applyBinomialBlurFilter(
      const typename Image<PixelT, c_dims>::Pointer& image)
{
  typedef Image<PixelT, c_dims> ImageT;
  typedef itk::BinomialBlurImageFilter<ImageT, ImageT> FilterType;
  static typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetRepetitions(1);
  // filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// applyMeanFilter(): template function
//-----------------------------------------------------------------------------

// template <class PixelT>
// local proc void
// doMeanFilter(Image<PixelT, c_dims>& image)
// {
// //     typedef itk::MeanImageFilter<ImageT, ImageT> FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     typename ImageT::SizeType indexRadius;
// //     indexRadius[0] = 5;
// //     indexRadius[1] = 5;
// //     indexRadius[2] = 5;
// //     filter->SetRadius(indexRadius);



//-----------------------------------------------------------------------------
// scalePixelValues(): template function
//-----------------------------------------------------------------------------

// This function scales the values of all the pixels of an ITK image by
// multiplying them all by specified value.

/*proc*/
template <class PixelT>
void
scalePixelValues(Image<PixelT, c_dims>& image,
                 double multiplier)
{
  if (multiplier == 1) return;

  std::cerr << "multiplier=" << multiplier << std::endl; //d

  // Make sure that we actually have the pixels loaded into the image:
  image.Update();

  typedef Image<PixelT, c_dims> ImageT;
  typedef ImageRegionIterator<ImageT> IterT;
  typedef typename ImageT::RegionType  RegionT;
  RegionT allOfImage = image.GetLargestPossibleRegion();
  IterT it(&image, allOfImage);
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    PixelT& pixel = it.Value();
    pixel *= multiplier;
  }
}


//-----------------------------------------------------------------------------
// reflectPixels(): template functiona
//-----------------------------------------------------------------------------

// This functions flips the image around the specified axes.  It does this by
// actually moving the pixels around, rather than by messing with the direction
// cosines.  It is, unfortunately, quite hairy.

// CAVEAT: This function assumes that RA is aligned with the i-axis, that Dec
// is aligned with the j-axis, and that V is aligned with the k-axis.  This is
// not always the case, however.  Sometimes a FITS image is rotated in RA/Dec
// space.  Furthermore, even then, this is only an approximation that works for
// small areas of the sky that are not near a pole.  As this function is only
// designed to be used with OsiriX, we can live with this caveat.

/*proc*/
template <class PixelT>
void
reflectPixels(Image<PixelT, c_dims>& image,
	      bool flipRAFlag, bool flipDecFlag, bool flipVFlag)
{
  // CAVEAT: This code will break if c_dims ever changes from 3.

  // Return immediately if a NOP is requested:
  if (!flipRAFlag and !flipDecFlag and !flipVFlag) return;

  // Make sure that we actually have the pixels loaded into the image:
  image.Update();

  typedef Image<PixelT, c_dims> ImageT;
  typedef typename ImageT::IndexType   IndexType;
  typedef typename ImageT::SizeType    SizeType;
  typedef typename ImageT::RegionType  RegionType;

  RegionType allOfImage = image.GetLargestPossibleRegion();
  SizeType imageSize = allOfImage.GetSize();
  IndexType imageOrigin = allOfImage.GetIndex();

  enum { c_i = FITSImageIO::c_i,
	 c_j = FITSImageIO::c_j,
	 c_k = FITSImageIO::c_k };

  const size_t minRaIndex  = imageOrigin[c_i];
  const size_t minDecIndex = imageOrigin[c_j];
  const size_t minVelIndex = imageOrigin[c_k];

  const size_t maxRaIndex  = minRaIndex  + imageSize[c_i] - 1;
  const size_t maxDecIndex = minDecIndex + imageSize[c_j] - 1;
  const size_t maxVelIndex = minVelIndex + imageSize[c_k] - 1;

  const size_t middleRaIndex  = minRaIndex  + ((imageSize[c_i] + 1) / 2) - 1;
  const size_t middleDecIndex = minDecIndex + ((imageSize[c_j] + 1) / 2) - 1;
  const size_t middleVelIndex = minVelIndex + ((imageSize[c_k] + 1) / 2) - 1;
  
  const bool lastFlipIsV = flipVFlag;
  const bool lastFlipIsDec = !lastFlipIsV and flipDecFlag;
  const bool lastFlipIsRA = !lastFlipIsV and !lastFlipIsDec and flipRAFlag;

  const size_t raStopIndex = lastFlipIsRA   ? middleRaIndex
                                            : maxRaIndex;
  const size_t decStopIndex = lastFlipIsDec ? middleDecIndex
                                            : maxDecIndex;
  const size_t velStopIndex = lastFlipIsV   ? middleVelIndex
                                            : maxVelIndex;

  const bool needToWorryAboutMiddleVel = flipVFlag and isOdd(imageSize[c_k]);
  const bool needToWorryAboutMiddleDec = flipDecFlag and isOdd(imageSize[c_j]);

  for (size_t vel_i = minVelIndex, velReverse_i = maxVelIndex;
       vel_i <= velStopIndex;
       ++vel_i, --velReverse_i)
    {
      for (size_t dec_i= minDecIndex, decReverse_i = maxDecIndex;
	   dec_i <= decStopIndex;
	   ++dec_i, --decReverse_i)
	{
	  for (size_t ra_i = minRaIndex, raReverse_i = maxRaIndex;
	       ra_i <= raStopIndex;
	       ++ra_i, --raReverse_i)
	    {
	      // Break out of the loop in situations in which we'd be swapping
	      // a pair of pixels back to their original locations due to
	      // swapping the pixels more than once:
	      if (needToWorryAboutMiddleVel) {
		if (vel_i == middleVelIndex) {
		  if (flipDecFlag) {
		    if (dec_i > middleDecIndex) break;
		    else if (needToWorryAboutMiddleDec and
			     flipRAFlag and
			     dec_i == middleDecIndex and
			     ra_i > middleRaIndex) {
		      break;
		    }
		  } else if (flipRAFlag and ra_i > middleRaIndex) break;
		}
	      }

	      IndexType thisPixelIndex;
	      IndexType oppositePixelIndex;
	      thisPixelIndex[c_i] = ra_i;
	      thisPixelIndex[c_j] = dec_i;
	      thisPixelIndex[c_k] = vel_i;
	      oppositePixelIndex[c_i] = flipRAFlag  ? raReverse_i  : ra_i;
	      oppositePixelIndex[c_j] = flipDecFlag ? decReverse_i : dec_i;
	      oppositePixelIndex[c_k] = flipVFlag   ? velReverse_i : vel_i;
	      PixelT tmp = image.GetPixel(thisPixelIndex);
	      image.SetPixel(thisPixelIndex,
			      image.GetPixel(oppositePixelIndex));
	      image.SetPixel(oppositePixelIndex, tmp);
	    }
	}
    }
}


//-----------------------------------------------------------------------------
// initializeChangeOfBasis(): private method of FITSImage template class
//-----------------------------------------------------------------------------

// /*private method*/
// template <typename ImageT>
// void
// FITSImage<ImageT>::initializeChangeOfBasis()
// {
//   const size_t dims = ImageT::ImageDimension;

//   // Set the origin to (0, 0, 0):
//   typename ImageT::PointType origin;
//   for (size_t d = 0; d < dims; ++d) origin[d] = 0;
//   _itkImage.SetOrigin(origin);

//   // Set the spacing vector to (1, 1, 1):
//   typename ImageT::SpacingType spacing;
//   for (size_t d = 0; d < dims; ++d) spacing[d] = 1;
//   _itkImage.SetSpacing(spacing);

//   // TODO: Replace these constants with something somewhere that is more
//   // globally accessible.
//   enum {ra, dec, v};
//   enum {l, p, s};

//   // Set the direction cosine matrix to properly orient RA, Dec, and V into LPS
//   // space:
//   typename ImageT::DirectionType direction;
//   direction(l, ra) = -1;
//   direction(p, v) = 1;
//   direction(s, dec) = 1;
//   _itkImage.SetDirection(direction);
// }


//-----------------------------------------------------------------------------
// raDecVToLpsMatrix(): private method of FITSImage template class
//-----------------------------------------------------------------------------

/*private method*/
template <typename ImageT>
HMatrix
FITSImage<ImageT>::raDecVToLpsMatrix()
{
  // TODO: Replace these constants with something somewhere that is more
  // globally accessible.
  enum {ra, dec, v};
  enum {l, p, s};

  HMatrix retval;
  retval(l, ra) = -1;
  retval(p, v) = 1;
  retval(s, dec) = 1;
  retval(3, 3) = 1;
  return retval;
}


//-----------------------------------------------------------------------------
// rightConcatenateTransformation(): template function
//-----------------------------------------------------------------------------

// This function right-concatenates `m` onto the coordinate frame
// transformation that is already stored in `image`.

/*proc*/ 
template <class ImageT>
void
rightConcatenateTransformation(ImageT& image, const HMatrix& m)
{
  setCoordinateFrameTransformation(
     image,
     getCoordinateFrameTransformation(image) * m
     );
}


//-----------------------------------------------------------------------------
// leftConcatenateTransformation(): template function
//-----------------------------------------------------------------------------

// This function left-concatenates `m` onto the coordinate frame
// transformation that is already stored in `image`.

/*proc*/ 
template <class ImageT>
void
leftConcatenateTransformation(ImageT& image, const HMatrix& m)
{
  setCoordinateFrameTransformation(
     image,
     getCoordinateFrameTransformation(image) * m
     );
}


//-----------------------------------------------------------------------------
// setCoordinateFrameTransformation(): template function
//-----------------------------------------------------------------------------

// This function sets the matrix stored within an ITK image that is used to
// convert coordinates from Index Space to LPS space.

// Note: There is alternative conception of the term "coordinate frame
// transformation", according to which the matrix, instead of transforming
// coordinates from Index Space to LPS Space would transform the basis vectors
// for Index Space to the basis vectors for LPS space.  The matrix for this
// alternative conception is merely the inverse of the matrix that we are using
// here.

/*proc*/ 
template <class ImageT>
void
setCoordinateFrameTransformation(ImageT& image, const HMatrix& m)
 			       
{
  typename ImageT::PointType origin;
  for (size_t d = 0; d < c_dims; ++d) origin[d] = m(d, 3);
  image.SetOrigin(origin);

  typename ImageT::SpacingType spacing;
  for (size_t col = 0; col < c_dims; ++col) {
    spacing[col] = sqrt(square(m(0, col)) +
			square(m(1, col)) +
			square(m(2, col)));
  }
  image.SetSpacing(spacing);

  Matrix<double, 3, 3> directionCosines;
  for (size_t row = 0; row < c_dims; ++row) {
    for (size_t col = 0; col < c_dims; ++col) {
      directionCosines(row, col) = m(row, col) / spacing[col]; 
    }
  }
  image.SetDirection(directionCosines);
}


//-----------------------------------------------------------------------------
// getCoordinateFrameTransformation(): template function
//-----------------------------------------------------------------------------

// This function gets the matrix stored within an ITK image that is used to
// convert coordinates from Index Space to LPS space.

// Note: There is alternative conception of the term "coordinate frame
// transformation", according to which the matrix, instead of transforming
// coordinates from Index Space to LPS Space would transform the basis vectors
// for Index Space to the basis vectors for LPS space.  The matrix for this
// alternative conception is merely the inverse of the matrix that we are using
// here.

/*proc*/ 
template <class ImageT>
HMatrix
getCoordinateFrameTransformation(const ImageT& image)
{
  typedef typename ImageT::PointType Point;
  typedef typename ImageT::SpacingType Spacing;
  typedef typename ImageT::DirectionType Direction;

  Point& origin = image.GetOrigin();
  Spacing& spacing = image.GetSpacing();
  Direction& direction = image.GetDirection();

  // Set the origin column of the retval matrix:
  HMatrix retval;
  for (size_t row = 0; row < c_dims; ++row) {
    retval(row, 3) = origin[row];
  }
  retval(3, 3) = 1;
  
  // Set the basis vector columns of the retval matrix:
  for (size_t row = 0; row < c_dims; ++row) {
    for (size_t col = 0; col < c_dims; ++col) {
      retval(row, col) = direction(row, col) * spacing[col];
    }
  }
  return retval;
}


//-----------------------------------------------------------------------------
// transformToNorthOrientedEquiangular():
//    private method of FITSImage template class
//-----------------------------------------------------------------------------

// /*method*/
// template <class ImageT>
// void
// FITSImage<ImageT>::transformToNorthOrientedEquiangular()
// {
//   enum {ra, dec};
//   Matrix m;

//   m(0, 0) = this->unitIInApproximateAngularSpace()[ra] * 1000 * 1000;
//   m(0, 1) = this->unitIInApproximateAngularSpace()[dec] * 1000 * 1000;
//   m(0, 2) = 0;

//   m(1, 0) = this->unitJInApproximateAngularSpace()[ra] * 1000 * 1000;
//   m(1, 1) = this->unitJInApproximateAngularSpace()[dec] * 1000 * 1000;
//   m(1, 2) = 0;

//   m(2, 0) = 0;
//   m(2, 1) = 0;
//   m(2, 2) = 1;

//   rightConcatinateTransformation(*this->getITKImage(), m);
// }


//-----------------------------------------------------------------------------
// ijkToNorthOrientedEquiangularMatrix():
//    private method of FITSImage template class
//-----------------------------------------------------------------------------

/*method*/
template <class ImageT>
HMatrix
FITSImage<ImageT>::ijkToNorthOrientedEquiangularMatrix()
{
  // YOU ARE HERE: I don't think that this is right at all.  Make it right.

  enum {ra, dec};
  HMatrix retval;

  retval(0, 0) = this->unitIInApproximateAngularSpace()[ra];
  retval(0, 1) = this->unitIInApproximateAngularSpace()[dec];
  retval(0, 2) = 0;
  retval(0, 3) = 0;

  retval(1, 0) = this->unitJInApproximateAngularSpace()[ra];
  retval(1, 1) = this->unitJInApproximateAngularSpace()[dec];
  retval(1, 2) = 0;
  retval(1, 3) = 0;

  retval(2, 0) = 0;
  retval(2, 1) = 0;
  retval(2, 2) = 1;
  retval(2, 3) = 0;

  retval(3, 0) = 0;
  retval(3, 1) = 0;
  retval(3, 2) = 0;
  retval(3, 3) = 1;
    
  return retval;
}


// //-----------------------------------------------------------------------------
// // transformToUnreorientedEquiangular(): template function
// //-----------------------------------------------------------------------------

// /*proc*/
// template <class ImageT>
// void
// transformToUnreorientedEquiangular(ImageT& image)
// {
  
//   const FITSImage<ImageT> info(image);
//   const double iAngleLen = info.unitIInApproximateAngularSpace().GetNorm();
//   const double jAngleLen = info.unitJInApproximateAngularSpace().GetNorm();

//   typename ImageT::SpacingType spacing;
//   spacing[0] = iAngleLen * 1000 * 1000;
//   spacing[1] = jAngleLen * 1000 * 1000;
//   spacing[2] = (iAngleLen + jAngleLen) / 2 * 1000 * 1000;

//   image.SetSpacing(spacing);
// }


// //-----------------------------------------------------------------------------
// // reorientNorth(): template function
// //-----------------------------------------------------------------------------

// /*proc*/
// template <class ImageT>
// void
// reorientNorth(ImageT& image)
// {
//   const FITSImage<ImageT> info(image);

//   // Multiply the direction cosine matrix by a rotation matrix to compensate
//   // for image rotation:
//   image.SetDirection(
//      image.GetDirection() *
//      rotationMatrix(-1 * info.rotationOfJFromIjkNorthVectorInDegrees())
//      );
// }


} // END namespace itk::fits::_internal


//-----------------------------------------------------------------------------
// Export symbols from _internal into itk::fits
//-----------------------------------------------------------------------------

using _internal::applyFlipImageFilter;
using _internal::applyBinomialBlurFilter;
using _internal::reflectPixels;
using _internal::scalePixelValues;
using _internal::writeImageInfo;

using _internal::setCoordinateFrameTransformation;
using _internal::getCoordinateFrameTransformation;
using _internal::rightConcatenateTransformation;
using _internal::leftConcatenateTransformation;

// using _internal::reorientNorth;
// using _internal::transformToUnreorientedEquiangular;

} } // END namespace itk::fits

#endif // _itkFITSImageUtils_txx
