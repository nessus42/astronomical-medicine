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

//-----------------------------------------------------------------------------
// loadFITSImageIO(): local function
//-----------------------------------------------------------------------------

local proc void*
loadFITSImageIO()
{
  static void* library = 0;
  if (!library) {
    const char* const libFilePath = pteJoinPath(pathToExecutableDir(),
                                                "libitkFITSImageIO.dylib");
    library = dlopen(libFilePath, RTLD_LAZY);
    if (!library) {
      runTimeError("Could not load libitkFITSImageIO.dylib");
    }
  }
  return library;
}


//-----------------------------------------------------------------------------
// fillMatrix(): function
//-----------------------------------------------------------------------------

proc void
fillMatrix(HMatrix& m, const double vals[4][4])
{
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      m(row, col) = vals[row][col];
    }
  }
}


//-----------------------------------------------------------------------------
// rotationMatrix(): function
//-----------------------------------------------------------------------------

proc HMatrix
rotationMatrix(double degrees)
{
  const double radians = degreesToRadians(degrees);
  const double s = sin(radians);
  const double c = cos(radians);
  HMatrix retval;
  retval.SetIdentity();
  retval[e_i][e_i] = c;
  retval[e_i][e_j] = -s;
  retval[e_j][e_i] = s;
  retval[e_j][e_j] = c;
  return retval;
}


//-----------------------------------------------------------------------------
// scalingMatrix(): function
//-----------------------------------------------------------------------------

proc HMatrix
scalingMatrix(double xScale, double yScale, double zScale)
{
  const double matrixVals[4][4] =
    { 
      xScale,   0,        0,       0,
      0,        yScale,   0,       0,
      0,        0,        zScale,  0,
      0,        0,        0,       1,
    };
  HMatrix retval;
  fillMatrix(retval, matrixVals);
  return retval;
}


//-----------------------------------------------------------------------------
// xyzToLpsMatrix(): function
//-----------------------------------------------------------------------------

proc const HMatrix&
xyzToLpsMatrix()
{
  static HMatrix retval;
  static bool firstTime = true;
  if (firstTime) {
    firstTime = false;
    retval(e_left,       e_x) = 1;
    retval(e_posterior,  e_z) = 1;
    retval(e_superior,   e_y) = 1;
    retval(3, 3) = 1;
  }
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
      dlsym(loadFITSImageIO(), "itkFITSImageIO_setNullValue")
      );
  if (!setNullValue) {
    runTimeError("Could not find function itkFITSImageIO_setNullValue()");
  }
  (*setNullValue)(nullValue);
}


//-----------------------------------------------------------------------------
// setFITSImageIODebugLevel(): function
//-----------------------------------------------------------------------------

proc void
setFITSImageIODebugLevel(int debugLevel)
{
  FITSImageIO::DebugLevelSetter setDebugLevel =
    reinterpret_cast<FITSImageIO::DebugLevelSetter>(
      dlsym(loadFITSImageIO(), "itkFITSImageIO_setDebugLevel")
      );
  if (!setDebugLevel) {
    runTimeError("Could not find function itkFITSImageIO_setDebugLevel()");
  }
  (*setDebugLevel)(debugLevel);
}

} } // END namespaces
