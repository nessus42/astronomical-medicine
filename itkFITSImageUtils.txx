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
#include <algorithm>

#include <wcs.h>

#include <itkImage.h>
#include <itkMetaDataObject.h>
#include <itkFITSImageIO.h>
#include <itkFITSImageUtils.h>

#include <da_usual.h>
using douglasAlan::runTimeError;

namespace itk {
  namespace fits {
    namespace _internal {

//-----------------------------------------------------------------------------
// Constructor of FITSImage
//-----------------------------------------------------------------------------

/*ctor*/
template <class ImageT>
FITSImage<ImageT>::FITSImage(const typename Self::Params& params)
  : _params(params),
    _itkImage(*params.itkImage),
    _wcsTransform(0),
    _unitKInVelocity(1),
    _velocityAtIjkOrigin(0),
    _rotationOfJFromIjkNorthVectorInDegrees(0),
    _rotationOfIFromIjkEastVectorInDegrees(0)
{
  assert(&_itkImage);
  initializeInstanceVars();

  assert(!(params.northUpP and (params.wcsP or params.equiangularP)) and
         (params.wcsP or !params.autoscaleZAxisP));

  // TODO: Clean up the comments below.

  // These are the coordinate transformations that we need to worry about
  // below:

  //    - IJK to North-up rotation
  //    - IJK-to-WCS transformation
  //    -    including velocity axis WCS
  //    - Equiangular adjustment to WCS
  //    - Autoscaling of velocity axis
  //    - Specified rotation
  //    - Scaling of X, Y, and Z axis
  //    - Conversion to LPS
  //    - RAS adjustment to LPS  // TODO: Do I really need this?

  // The matrix below is a conversion of the coordinates, not the frame, and
  // therefore, will affect the original coordinates starting from the last
  // matrix in the product chain first.  Also keep in mind that these
  // transformations are the inverses of the frame transformations.

  // TODO: Find out from someone what the term is for the transformation that
  // transforms the coordinates rather than the basis vectors and change any
  // code and coments (and my little primer) to reflect this terminology.

  // Set the ITK Image's coordinate transformation parameters:
  {
    HMatrix ijkToXyzMatrix;
    if (params.equiangularP) {
      ijkToXyzMatrix = this->ijkToEquiangularMatrix();
    } else if (params.wcsP) {
      ijkToXyzMatrix = this->ijkToWcsMatrix();
    } else if (params.northUpP) {
      ijkToXyzMatrix = this->ijkToNorthUpMatrix();
    } else {
      ijkToXyzMatrix.SetIdentity();
    }

    HMatrix cft =
      xyzToLpsMatrix()
      * scalingMatrix(params.xAxisScale,
                      params.yAxisScale,
                      params.zAxisScale)
      * rotationMatrix(params.rotateSky)
      * scalingMatrix(1, 1,
                      this->zAxisAutoscale(
                         params.autoscaleZAxisP,
                         ijkToXyzMatrix))
      * ijkToXyzMatrix;

    setCoordinateFrameTransformation(_itkImage, cft);
  }
}


//-----------------------------------------------------------------------------
// initializeInstanceVars(): private method of FITSImage template class
//-----------------------------------------------------------------------------

template <class ImageT>
void 
FITSImage<ImageT>::initializeInstanceVars()
{
  this->setVelocityInstanceVars();

  typename ImageT::RegionType
    allOfImage = _itkImage.GetLargestPossibleRegion();
  typename ImageT::SizeType imageSize = allOfImage.GetSize();
  typename ImageT::IndexType imageOrigin = allOfImage.GetIndex();

  _ijkSize[e_i] = imageSize[e_i];
  _ijkSize[e_j] = imageSize[e_j];
  _ijkSize[e_k] = imageSize[e_k];
  _ijkCenter[e_i] = imageOrigin[e_i] + imageSize[e_i]/2.0 - 0.5;
  _ijkCenter[e_j] = imageOrigin[e_j] + imageSize[e_j]/2.0 - 0.5;
  _ijkCenter[e_k] = imageOrigin[e_k] + imageSize[e_k]/2.0 - 0.5;

  WcsTransformConstPtr wcs = makeWcsTransform(); 
  if (wcs) {
    _wcsTransform = wcs;
    const IjkPoint zeroZeroZero;
    _wcsImageOrigin = wcs->TransformPoint(zeroZeroZero);
    _wcsImageOrigin[e_vel] = _velocityAtIjkOrigin;
    _wcsImageCenter = wcs->TransformPoint(_ijkCenter);
    _wcsImageCenter[e_vel] = 
      _velocityAtIjkOrigin + _ijkCenter[e_k] * _unitKInVelocity;

    // TODO: It's confusing that in some situations V is 1 and dec is 2, and
    // in others, dec is 1 and V is 2.  We need a better way to denote this.

    // Calculate lengths of i and j vectors in RA/Dec space:
    { 
      IjkPoint ijkLeftHalfAPixel  = _ijkCenter;
      IjkPoint ijkRightHalfAPixel = _ijkCenter;
      IjkPoint ijkDownHalfAPixel  = _ijkCenter;
      IjkPoint ijkUpHalfAPixel    = _ijkCenter;
      ijkLeftHalfAPixel[e_i]  -= .5;
      ijkRightHalfAPixel[e_i] += .5;
      ijkDownHalfAPixel[e_j]  -= .5;
      ijkUpHalfAPixel[e_j]    += .5;

      WcsPoint wcsLeftHalfAPixel  = wcs->TransformPoint(ijkLeftHalfAPixel);
      WcsPoint wcsRightHalfAPixel = wcs->TransformPoint(ijkRightHalfAPixel);
      WcsPoint wcsDownHalfAPixel  = wcs->TransformPoint(ijkDownHalfAPixel);
      WcsPoint wcsUpHalfAPixel    = wcs->TransformPoint(ijkUpHalfAPixel);

      _unitIInWcs = wcsRightHalfAPixel - wcsLeftHalfAPixel;
      _unitJInWcs = wcsUpHalfAPixel - wcsDownHalfAPixel;
    }

    _raAngularScalingFactor = cos(degreesToRadians(_wcsImageCenter[e_dec]));
    _unitIInApproximateAngularSpace = _unitIInWcs;
    _unitJInApproximateAngularSpace = _unitJInWcs;
    _unitIInApproximateAngularSpace[e_ra] *= _raAngularScalingFactor;
    _unitJInApproximateAngularSpace[e_ra] *= _raAngularScalingFactor;

    // Calculate the image's rotation by finding a northward-oriented vector in
    // WCS space, and then transforming it into IJK space.  We can then use trig
    // to determine the amount of rotation of the image from north in IJK space:
    {
      WcsTransformPtr inverseWcs = WCS::New();
      wcs->GetInverse(inverseWcs);

      // We calculate wcsNorthVector, just for the purpose of getting a Dec
      // increment that is about the size of a pixel or so.  We could, of
      // course, accomplish much the same thing just by selecting an arbitrary
      // small increment northward, but then it might be too small or it might
      // be too big for the image:
      WcsVector wcsNorthVector =
        fabs(_unitJInWcs[e_dec]) > fabs(_unitIInWcs[e_dec])
        ? _unitJInWcs
        : _unitIInWcs;
      wcsNorthVector[e_ra] = 0;
      wcsNorthVector[e_dec] = fabs(wcsNorthVector[e_dec]);

      // We now add the northward increment vector to our center point
      // expressed in RA/Dec coordinates.  This gives us a slightly northward
      // point in WCS coordinates, and we then use an inverse WCS transform to
      // get the IJK pixel position (with floating point IJK index values) of
      // this fastidiously calculated northward point:
      WcsPoint wcsPointNorthOfCenter = _wcsImageCenter + wcsNorthVector;
      IjkPoint ijkPointNorthOfCenter =
        inverseWcs->TransformPoint(wcsPointNorthOfCenter);

      // We then subtract the center point in IJK coordinates from the
      // northward point to get a northward pointing vector in IJK coordinates:
      _ijkNorthVector = ijkPointNorthOfCenter - _ijkCenter;
      _ijkNorthVector[e_vel] = 0;
      _ijkNorthVector.Normalize();

      // And finally, with the northward pointing vector in IJK coordinates, we
      // can determine how much the image is rotated from north
      // (counterclockwise):
      _rotationOfJFromIjkNorthVectorInDegrees =
        radiansToDegrees(atan2(-1 * _ijkNorthVector[e_i],
                               _ijkNorthVector[e_j]));

      // We're now going to do all the same stuff as above, only for east,
      // rather than north:
      WcsVector wcsEastVector =
        fabs(_unitJInWcs[e_ra]) > fabs(_unitIInWcs[e_ra])
        ? _unitJInWcs
        : _unitIInWcs;
      wcsEastVector[e_ra] = fabs(wcsEastVector[e_ra]);
      wcsEastVector[e_dec] = 0;
      WcsPoint wcsPointEastOfCenter = _wcsImageCenter + wcsEastVector;
      IjkPoint ijkPointEastOfCenter =
        inverseWcs->TransformPoint(wcsPointEastOfCenter);
      _ijkEastVector = ijkPointEastOfCenter - _ijkCenter;
      _ijkEastVector[e_vel] = 0;
      _ijkEastVector.Normalize();
      _rotationOfIFromIjkEastVectorInDegrees =
        radiansToDegrees(atan2(-1 * _ijkEastVector[e_j],
                               -1 * _ijkEastVector[e_i]));

//       // Calculate the rotation of the east unit vector from the north unit
//       // vector by calculating the dot product:
//       _rotationOfIjkEastVectorFromIjkNorthVector =
//         acos(_ijkEastVector[c_i] * _ijkNorthVector[c_i] +
//              _ijkEastVector[c_j] * _ijkEastVector[c_j]);
    }
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

  // Fetch the FITS header out of the metadata dictionary.  We cast away const
  // here conly because ExposeMetaData() is unfortunately not const correct:
  string fitsHeader;
  ExposeMetaData(const_cast<MetaDataDictionary&>(mdd),
		 "FITS Header",
		 fitsHeader);
  assert(!fitsHeader.empty());

  // TODO: Add more error checking here in case 'wcs' ends up in some sort of
  // erroneous state.
  const ConstRcMallocPointer<WorldCoor> wcs = wcsinit(fitsHeader.c_str());
  if (wcs) {
    WcsTransformPtr retval = WCS::New();
    retval->SetWCS(wcs);
    WcsTransformConstPtr constRetval (retval);
    return constRetval;
  } else {
    return 0;
  }
}


//-----------------------------------------------------------------------------
// setVelocityInstanceVars(): private method of FITSImage template class
//-----------------------------------------------------------------------------

/*proc*/
template <class ImageT>
void
FITSImage<ImageT>::setVelocityInstanceVars()
{
  using itk::MetaDataDictionary;
  using itk::MetaDataObject;
  const MetaDataDictionary& mdd = _itkImage.GetMetaDataDictionary();

  // Fetch the velocity information out of the metadata dictionary.  We cast
  // away const here conly because ExposeMetaData() is unfortunately not const
  // correct:
  double velocityAtIndexOrigin;
  double velocityDelta = 0;
  ExposeMetaData(const_cast<MetaDataDictionary&>(mdd),
		 "fits2itk.velocityAtIndexOrigin",
                 velocityAtIndexOrigin);
  ExposeMetaData(const_cast<MetaDataDictionary&>(mdd),
                 "fits2itk.velocityDelta",
                 velocityDelta);
  if (velocityDelta != 0) {
    _unitKInVelocity = velocityDelta;
    _velocityAtIjkOrigin = velocityAtIndexOrigin;
  }
}


//-----------------------------------------------------------------------------
// ijkToEquiangularMatrix():
//    nonvirtual method of FITSImage template class
//-----------------------------------------------------------------------------

/*method*/
template <class ImageT>
HMatrix
FITSImage<ImageT>::ijkToEquiangularMatrix() const
{
  HMatrix retval;
  retval.SetIdentity();
  if (!_wcsTransform) return retval;
  retval(e_ra,  e_i) = _unitIInApproximateAngularSpace[e_ra];
  retval(e_dec, e_i) = _unitIInApproximateAngularSpace[e_dec];
  retval(e_ra,  e_j) = _unitJInApproximateAngularSpace[e_ra];
  retval(e_dec, e_j) = _unitJInApproximateAngularSpace[e_dec];
  return retval;
}


//-----------------------------------------------------------------------------
// ijkToWcsMatrix():
//    nonvirtual method of FITSImage template class
//-----------------------------------------------------------------------------

/*method*/
template <class ImageT>
HMatrix
FITSImage<ImageT>::ijkToWcsMatrix() const
{
  HMatrix retval;
  retval.SetIdentity();
  if (!_wcsTransform) return retval;
  retval(e_ra,  e_i) = _unitIInWcs[e_ra];
  retval(e_dec, e_i) = _unitIInWcs[e_dec];
  retval(e_vel, e_i) = _unitIInWcs[e_vel];
  retval(e_ra,  e_j) = _unitJInWcs[e_ra];
  retval(e_dec, e_j) = _unitJInWcs[e_dec];
  retval(e_vel, e_j) = _unitJInWcs[e_vel];
  retval(e_ra, c_dims) = _wcsImageCenter[e_ra];
  retval(e_dec, c_dims) = _wcsImageCenter[e_dec];
  retval(e_vel, c_dims) = _wcsImageCenter[e_vel];
  return retval;
}


//-----------------------------------------------------------------------------
// ijkToNorthUpMatrix(): 
//    nonvirtual method of FITSImage tremplace class
//-----------------------------------------------------------------------------

template <class ImageT>
HMatrix
FITSImage<ImageT>::ijkToNorthUpMatrix() const
{
  return rotationMatrix(_rotationOfJFromIjkNorthVectorInDegrees);
}


//-----------------------------------------------------------------------------
// zAxisAutoscale():
//    private method of FITSImage template class
//-----------------------------------------------------------------------------

/*proc*/ 
template <class ImageT>
double
FITSImage<ImageT>::zAxisAutoscale(bool autoscaleZAxisP,
                                  const HMatrix& ijkToXyzMatrix)
{
  if (!autoscaleZAxisP) return 1;
  
  HVector ijkSizeVector;
  ijkSizeVector[e_i] = _ijkSize[e_i];
  ijkSizeVector[e_j] = _ijkSize[e_j];
  ijkSizeVector[e_k] = _ijkSize[e_k];

  HVector xyzSizeVector = ijkToXyzMatrix * ijkSizeVector;
  double xyMax = std::max(fabs(xyzSizeVector[e_x]), fabs(xyzSizeVector[e_y]));
  return xyMax / xyzSizeVector[e_z];
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

  out << "Pixel coordinates:\n";
  out << "   Image size ("
      << image.ijkSize()[e_i] << ", "
      << image.ijkSize()[e_j] << ", "
      << image.ijkSize()[e_k] << ")\n"
      << "   Image origin: (0, 0, 0)\n"
      << "   Center pixel: ("
      << image.ijkCenter()[e_ra] << ", "
      << image.ijkCenter()[e_dec] << ", "
      << image.ijkCenter()[e_vel] << ")\n";
  
  if (image.wcsTransform()) {
    out 
      << "RA/Dec/Velocity:\n"
      << "   Image origin: ("
      << image.wcsImageOrigin()[e_ra] << ", "
      << image.wcsImageOrigin()[e_dec] << ", "
      << image.wcsImageOrigin()[e_vel] << ")\n"
      << "   Center pixel: ("
      << image.wcsImageCenter()[e_ra] << ", "
      << image.wcsImageCenter()[e_dec] << ", "
      << image.wcsImageCenter()[e_vel] << ")\n"
      << "   Pixel width: ("
      << image.unitIInWcs()[e_ra] << ", "
      << image.unitIInWcs()[e_dec] << ")\n"
      << "   Pixel height: ("
      << image.unitJInWcs()[e_ra] << ", "
      << image.unitJInWcs()[e_dec] << ")\n"
      << "   Pixel depth: " << image.unitKInVelocity() << " km/s\n"
      << "Approximate Angular Space:\n"
      << "   Pixel width: ("
      << image.unitIInApproximateAngularSpace()[e_ra] << ", "
      << image.unitIInApproximateAngularSpace()[e_dec] << ")\n"
      << "   Pixel height: ("
      << image.unitJInApproximateAngularSpace()[e_ra] << ", "
      << image.unitJInApproximateAngularSpace()[e_dec] << "\n"
      << "Approximate Angular Distance:\n"
      << "   Pixel width: "
      << image.unitIInApproximateAngularSpace().GetNorm() << "\n"
      << "   Pixel height: "
      << image.unitJInApproximateAngularSpace().GetNorm() << "\n"
      << "Orientation of image:\n"
      << "   Unit north vector in pixel coordinates: ("
      << image.ijkNorthVector()[e_i] << ", "
      << image.ijkNorthVector()[e_j] << ")\n"
      << "   Rotation of unit north vector clockwise:  "
      << image.rotationOfJFromIjkNorthVectorInDegrees() << " degrees\n"
      << "   Rotation of unit east vector clockwise:  "
      << image.rotationOfIFromIjkEastVectorInDegrees() << " degrees\n"
      << "   Angle between unit north vector and unit east vector: " 
      << (image.rotationOfJFromIjkNorthVectorInDegrees()
          - image.rotationOfIFromIjkEastVectorInDegrees()
          + 90)
      << " degrees\n"
      << "Pixel to WCS coordinate transformation matrix:";

    HMatrix t = image.ijkToWcsMatrix();
    for (int row = 0; row <= c_dims; ++row) {
      out << "\n  ";
      for (int col = 0; col <= c_dims; ++col) {
        out << " " << t(row, col);
      }
    }
    out << endl;

  } else {
    out << "There is no WCS information in the input image." << endl;
  }

//       << "Direction cosines:\n"
//       << itkImage.GetDirection()
//       << "Image spacing: " << itkImage.GetSpacing() << "\n"
//       << "Image origin: " << itkImage.GetOrigin() << "\n"

}


//-----------------------------------------------------------------------------
// writeFitsHeader(): template function
//-----------------------------------------------------------------------------

/*proc*/
template <class ImageT>
void
writeFitsHeader(const ImageT& image, ostream& out)
{
  using itk::MetaDataDictionary;
  using itk::MetaDataObject;
  const MetaDataDictionary& mdd = image.GetMetaDataDictionary();
  string fitsHeader;
  ExposeMetaData(const_cast<MetaDataDictionary&>(mdd),
		 "FITS Header",
		 fitsHeader);
  out << fitsHeader;
}


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

  const size_t minRaIndex  = imageOrigin[e_i];
  const size_t minDecIndex = imageOrigin[e_j];
  const size_t minVelIndex = imageOrigin[e_k];

  const size_t maxRaIndex  = minRaIndex  + imageSize[e_i] - 1;
  const size_t maxDecIndex = minDecIndex + imageSize[e_j] - 1;
  const size_t maxVelIndex = minVelIndex + imageSize[e_k] - 1;

  const size_t middleRaIndex  = minRaIndex  + ((imageSize[e_i] + 1) / 2) - 1;
  const size_t middleDecIndex = minDecIndex + ((imageSize[e_j] + 1) / 2) - 1;
  const size_t middleVelIndex = minVelIndex + ((imageSize[e_k] + 1) / 2) - 1;
  
  const bool lastFlipIsV = flipVFlag;
  const bool lastFlipIsDec = !lastFlipIsV and flipDecFlag;
  const bool lastFlipIsRA = !lastFlipIsV and !lastFlipIsDec and flipRAFlag;

  const size_t raStopIndex = lastFlipIsRA   ? middleRaIndex
                                            : maxRaIndex;
  const size_t decStopIndex = lastFlipIsDec ? middleDecIndex
                                            : maxDecIndex;
  const size_t velStopIndex = lastFlipIsV   ? middleVelIndex
                                            : maxVelIndex;

  const bool needToWorryAboutMiddleVel = flipVFlag and isOdd(imageSize[e_k]);
  const bool needToWorryAboutMiddleDec = flipDecFlag and isOdd(imageSize[e_j]);

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
	      thisPixelIndex[e_i] = ra_i;
	      thisPixelIndex[e_j] = dec_i;
	      thisPixelIndex[e_k] = vel_i;
	      oppositePixelIndex[e_i] = flipRAFlag  ? raReverse_i  : ra_i;
	      oppositePixelIndex[e_j] = flipDecFlag ? decReverse_i : dec_i;
	      oppositePixelIndex[e_k] = flipVFlag   ? velReverse_i : vel_i;
	      PixelT tmp = image.GetPixel(thisPixelIndex);
	      image.SetPixel(thisPixelIndex,
			      image.GetPixel(oppositePixelIndex));
	      image.SetPixel(oppositePixelIndex, tmp);
	    }
	}
    }
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

} } } // END namespace ::itk::fits::_internal

#endif // _itkFITSImageUtils_txx
