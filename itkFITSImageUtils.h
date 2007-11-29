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

#include <itkFITSImageIO.h>

// BEGIN
namespace itk {
namespace fits {

  // Global functions:

  /*proc*/ void
  setNullValue(double nullValue);

  /*proc*/ FITSImageIO::WCSTransform::ConstPointer
  deprecated_getWCSTransform();

  /*proc*/ void
  setNullValue(double);

namespace _internal {

  // Internal functions:

  /*internal proc*/ void*
  loadFITSImageIO();

  /*internal proc*/ Matrix<double, 3, 3>
  rotationMatrix(double degrees);

  /*internal proc*/ void*
  deprecated_getWCSTransform();

} // END namespace _internal


} } // END namespace itk::fits


#ifndef ITK_MANUAL_INSTANTIATION
#include <itkFITSImageUtils.txx>
#endif

#endif // _itkFITSImageUtils_h
