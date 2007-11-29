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
#include <itkFITSImageIO.h>
#include <itkMatrix.h>
#include <pathToExecutable.h>
#include <da_sugar.h>

namespace itk {
namespace fits {
namespace _internal {


//-----------------------------------------------------------------------------
// loadFITSImageIO(): local function
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
// deprecated_getWCSTransform(): internal function
//-----------------------------------------------------------------------------

proc void*
deprecated_getWCSTransform()
{
  FITSImageIO::WCSTransformGetter getTransform =
    reinterpret_cast<FITSImageIO::WCSTransformGetter>(
      dlsym(loadFITSImageIO(),
	    "itkFITSImageIO_deprecatedGetWCSTransform")
      );
  if (!getTransform) {
    runTimeError("Could not find function "
		 "itkFITSImageIO_deprecatedGetWCSTransform()");
  }
  return (*getTransform)();
}


//-----------------------------------------------------------------------------
// rotationMatrix(): internal function
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

} // END namespace _internal


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

} } // END namespace itk::fits
