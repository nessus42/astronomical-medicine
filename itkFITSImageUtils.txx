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
#include <itkMatrix.h>
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
template <class ImageType>
FITSImage<ImageType>::FITSImage(const typename Self::Params& params)
  : _params(params), _itkImage(*params.itkImage)
{
  assert(&_itkImage);

  typename ImageType::RegionType
    allOfImage = _itkImage.GetLargestPossibleRegion();
  typename ImageType::SizeType imageSize = allOfImage.GetSize();
  typename ImageType::IndexType imageOrigin = allOfImage.GetIndex();

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


  // Set the ITK Image's coordinate transformation parameters:
  {

    setCoordinateFrameTransformation(
       _itkImage,
       this->raDecVToLpsMatrix() *
       this->ijkToNorthOrientedEquiangularMatrix());
  }
}


//-----------------------------------------------------------------------------
// makeWcsTransform(): private static method of FITSImage template class
//-----------------------------------------------------------------------------

/*proc*/
template <class ImageType>
typename FITSImage<ImageType>::WcsTransformConstPtr
FITSImage<ImageType>::makeWcsTransform()
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
template <class ImageType>
void
writeImageInfo(const FITSImage<ImageType>& image, ostream& out)
{
  const ImageType& itkImage = *image.getITKImage();
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
template <class PixelType>
typename Image<PixelType, c_dims>::Pointer
applyFlipImageFilter(const typename
		     Image<PixelType, c_dims>::Pointer& image)
{
  typedef Image<PixelType, c_dims> ImageType;
  typedef itk::FlipImageFilter<ImageType> FilterType;
  typedef typename FilterType::FlipAxesArrayType FlipAxesArrayType;
  typename FilterType::Pointer filter = FilterType::New();
  FlipAxesArrayType flipArray;
  flipArray[0] = 0;
  flipArray[1] = 1;
  flipArray[2] = 0;
  filter->SetFlipAxes(flipArray);
  filter->SetInput(image);
  filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// applyBinomialBlur(): template function
//-----------------------------------------------------------------------------

/*proc*/ 
template <class PixelType>
typename Image<PixelType, c_dims>::Pointer
applyBinomialBlurFilter(const typename Image<PixelType, c_dims>::Pointer& image)
{
  typedef Image<PixelType, c_dims> ImageType;
  typedef itk::BinomialBlurImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetRepetitions(1);
  filter->Update();
  return filter->GetOutput();
}


//-----------------------------------------------------------------------------
// applyMeanFilter(): template function
//-----------------------------------------------------------------------------

// template <class PixelType>
// local proc void
// doMeanFilter(Image<PixelType, c_dims>& image)
// {
// //     typedef itk::MeanImageFilter<ImageType, ImageType> FilterType;
// //     typename FilterType::Pointer filter = FilterType::New();
// //     typename ImageType::SizeType indexRadius;
// //     indexRadius[0] = 5;
// //     indexRadius[1] = 5;
// //     indexRadius[2] = 5;
// //     filter->SetRadius(indexRadius);


//-----------------------------------------------------------------------------
// reflectPixels(): template function
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
template <class PixelType>
void
reflectPixels(Image<PixelType, c_dims>& image,
	      bool flipRAFlag, bool flipDecFlag, bool flipVFlag)
{
  // CAVEAT: This code will break if c_dims ever changes from 3.

  // Return immediately if a NOP is requested:
  if (!flipRAFlag and !flipDecFlag and !flipVFlag) return;

  // Make sure that we actually have the pixels loaded into the image:
  image.Update();

  typedef Image<PixelType, c_dims> ImageType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType SizeType;

  typename ImageType::RegionType allOfImage = image.GetLargestPossibleRegion();
  typename ImageType::SizeType imageSize = allOfImage.GetSize();
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
	      PixelType tmp = image.GetPixel(thisPixelIndex);
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
// template <typename ImageType>
// void
// FITSImage<ImageType>::initializeChangeOfBasis()
// {
//   const size_t dims = ImageType::ImageDimension;

//   // Set the origin to (0, 0, 0):
//   typename ImageType::PointType origin;
//   for (size_t d = 0; d < dims; ++d) origin[d] = 0;
//   _itkImage.SetOrigin(origin);

//   // Set the spacing vector to (1, 1, 1):
//   typename ImageType::SpacingType spacing;
//   for (size_t d = 0; d < dims; ++d) spacing[d] = 1;
//   _itkImage.SetSpacing(spacing);

//   // TODO: Replace these constants with something somewhere that is more
//   // globally accessible.
//   enum {ra, dec, v};
//   enum {l, p, s};

//   // Set the direction cosine matrix to properly orient RA, Dec, and V into LPS
//   // space:
//   typename ImageType::DirectionType direction;
//   direction(l, ra) = -1;
//   direction(p, v) = 1;
//   direction(s, dec) = 1;
//   _itkImage.SetDirection(direction);
// }


//-----------------------------------------------------------------------------
// raDecVToLpsMatrix(): private method of FITSImage template class
//-----------------------------------------------------------------------------

/*private method*/
template <typename ImageType>
Matrix
FITSImage<ImageType>::raDecVToLpsMatrix()
{
  // TODO: Replace these constants with something somewhere that is more
  // globally accessible.
  enum {ra, dec, v};
  enum {l, p, s};

  Matrix retval;
  retval(l, ra) = -1;
  retval(p, v) = 1;
  retval(s, dec) = 1;
  return retval;
}


//-----------------------------------------------------------------------------
// rightConcatinateTransformation(): internal template function
//-----------------------------------------------------------------------------

// /*internal proc*/ 
// template <class ImageType>
// void
// rightConcatinateTransformation(ImageType& image, const Matrix& m)
 			       
// {
//   typename ImageType::SpacingType spacingMultiplier;

//   for (size_t col = 0; col < 3; ++col) {
//     spacingMultiplier[col] = sqrt(square(m(0, col)) +
// 				  square(m(1, col)) +
// 				  square(m(2, col)));
//   }

//   // Do an element-by-element multiplication:
//   image.SetSpacing(image.GetSpacing() * spacingMultiplier);

//   Matrix directionMultiplier;
//   for (size_t row = 0; row < 3; ++row) {
//     for (size_t col = 0; col < 3; ++col) {
//       directionMultiplier(row, col) = m(row, col) / spacingMultiplier[col]; 
//     }
//   }
//   image.SetDirection(image.GetDirection() * directionMultiplier);
// }


//-----------------------------------------------------------------------------
// setCoordinateFrameTransformation(): internal template function
//-----------------------------------------------------------------------------

// This function sets the matrix that is used to convert coordinates from Index
// Space to LPS space.

// Note: There is alternative conception of the term "coordinate frame
// transformation", according to which the matrix, instead of transforming
// coordinates from Index Space to LPS Space would transform the basis vectors
// for Index Space to the basis vectors for LPS space.  The matrix for this
// alternative conception is merely the inverse of the matrix that we are using
// here.

/*internal proc*/ 
template <class ImageType>
void
setCoordinateFrameTransformation(ImageType& image, const Matrix& m)
 			       
{
  const size_t dims = ImageType::ImageDimension;
  const size_t mDims = 3;

  // At the moment, we don't preserve the origin at all, so we'll just set it
  // to (0, 0, 0):
  typename ImageType::PointType origin;
  for (size_t d = 0; d < dims; ++d) origin[d] = 0;
  image.SetOrigin(origin);

  typename ImageType::SpacingType spacing;
  for (size_t col = 0; col < dims; ++col) {
    spacing[col] = sqrt(square(m(0, col)) +
			square(m(1, col)) +
			square(m(2, col)));
  }
  image.SetSpacing(spacing);

  Matrix directionCosines;
  for (size_t row = 0; row < mDims; ++row) {
    for (size_t col = 0; col < mDims; ++col) {
      directionCosines(row, col) = m(row, col) / spacing[col]; 
    }
  }
  image.SetDirection(directionCosines);
}


//-----------------------------------------------------------------------------
// transformToNorthOrientedEquiangular():
//    private method of FITSImage template class
//-----------------------------------------------------------------------------

// /*method*/
// template <class ImageType>
// void
// FITSImage<ImageType>::transformToNorthOrientedEquiangular()
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
template <class ImageType>
Matrix
FITSImage<ImageType>::ijkToNorthOrientedEquiangularMatrix()
{
  enum {ra, dec};
  Matrix retval;

  retval(0, 0) = this->unitIInApproximateAngularSpace()[ra] * 1000 * 1000;
  retval(0, 1) = this->unitIInApproximateAngularSpace()[dec] * 1000 * 1000;
  retval(0, 2) = 0;

  retval(1, 0) = this->unitJInApproximateAngularSpace()[ra] * 1000 * 1000;
  retval(1, 1) = this->unitJInApproximateAngularSpace()[dec] * 1000 * 1000;
  retval(1, 2) = 0;

  retval(2, 0) = 0;
  retval(2, 1) = 0;
  // retval(2, 2) = 1;  //d This has been replaced just as a kludge until we
                        //d put back in some sort of auto-scaling for V.
  retval(2, 2) = 10000; //d

  // YOU ARE HERE

  return retval;
}


// //-----------------------------------------------------------------------------
// // transformToUnreorientedEquiangular(): template function
// //-----------------------------------------------------------------------------

// /*proc*/
// template <class ImageType>
// void
// transformToUnreorientedEquiangular(ImageType& image)
// {
  
//   const FITSImage<ImageType> info(image);
//   const double iAngleLen = info.unitIInApproximateAngularSpace().GetNorm();
//   const double jAngleLen = info.unitJInApproximateAngularSpace().GetNorm();

//   typename ImageType::SpacingType spacing;
//   spacing[0] = iAngleLen * 1000 * 1000;
//   spacing[1] = jAngleLen * 1000 * 1000;
//   spacing[2] = (iAngleLen + jAngleLen) / 2 * 1000 * 1000;

//   image.SetSpacing(spacing);
// }


// //-----------------------------------------------------------------------------
// // reorientNorth(): template function
// //-----------------------------------------------------------------------------

// /*proc*/
// template <class ImageType>
// void
// reorientNorth(ImageType& image)
// {
//   const FITSImage<ImageType> info(image);

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
using _internal::writeImageInfo;

// using _internal::reorientNorth;
// using _internal::transformToUnreorientedEquiangular;

} } // END namespace itk::fits

#endif // _itkFITSImageUtils_txx
