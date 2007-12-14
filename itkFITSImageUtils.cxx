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
#include <dlfcn.h>                         // For dlopen()
#include <pathToExecutable.h>
#include <itkFITSImageIO.h>
#include <itkFITSImageUtils.h>
#include <da_sugar.h>

namespace itk {
namespace fits {
namespace _internal {


//-----------------------------------------------------------------------------
// loadFITSImageIO(): internal function
//-----------------------------------------------------------------------------

local proc void*
loadFITSImageIO()
{
  static void* library = 0;
  if (!library) {
    const char* const libFilePath = pteJoinPath(pathToExecutableDir(),
						"libitkFITSImageIO.so");
    library = dlopen(libFilePath, RTLD_LAZY);
    if (!library) {
      runTimeError("Could not load libitkFITSImageIO.so");
    }
  }
  return library;
}


//-----------------------------------------------------------------------------
// fillMatrix(): internal function
//-----------------------------------------------------------------------------

/*internal proc*/
void
fillMatrix(Matrix& m, const double vals[3][3])
{
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      m(row, col) = vals[row][col];
    }
  }
}

//-----------------------------------------------------------------------------
// rotationMatrix(): internal function
//-----------------------------------------------------------------------------

/*internal proc*/
Matrix
rotationMatrix(double degrees)
{
  const double radians = degreesToRadians(degrees);
  const double s = sin(radians);
  const double c = cos(radians);
  const double matrixVals[3][3] =
    { 
        c, 0, s,
	0, 1, 0,
	s, 0, c 
    };
  Matrix retval;
  fillMatrix(retval, matrixVals);
  return retval;
}


//-----------------------------------------------------------------------------
// setNullValue(): function
//-----------------------------------------------------------------------------

proc void
setNullValue(double nullValue)
{
  FITSImageIO::NullValueSetter setNullValue =
    reinterpret_cast<FITSImageIO::NullValueSetter>(
      dlsym(_internal::loadFITSImageIO(), "itkFITSImageIO_setNullValue")
      );
  if (!setNullValue) {
    runTimeError("Could not find function itkFITSImageIO_setNullValue()");
  }
  (*setNullValue)(nullValue);
}


} } } // END namespaces
