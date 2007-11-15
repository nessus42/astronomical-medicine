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

#include <cmath>
#include <itkMatrix.h>
#include <da_sugar.h>

namespace itk {
namespace fits {
namespace __internal {


//-----------------------------------------------------------------------------
// rotationMatrix(): internal proc
//-----------------------------------------------------------------------------

Matrix<double, 3, 3>
rotationMatrix(double degrees)
{
  const double s = sin(degrees/180 * M_PI);
  const double c = cos(degrees/180 * M_PI);
  Matrix<double, 3, 3> retval;
  retval(0, 0) = c;
  retval(0, 1) = -s;
  retval(1, 0) = s;
  retval(1, 1) = c;
  retval(2, 2) = 1;
  return retval;
}


} } } // END namespace itk::fits::__internal
