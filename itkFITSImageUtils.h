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
#ifndef __itkFITSImageUtils_h
#define __itkFITSImageUtils_h


// BEGIN
namespace itk {
namespace fits {
namespace __internal {


/*internal proc*/ Matrix<double, 3, 3>
rotationMatrix(double degrees);

} } } // END namespace itk::fits::__internal


#ifndef ITK_MANUAL_INSTANTIATION
#include <itkFITSImageUtils.txx>
#endif

#endif // __itkFITSImageUtils_h
