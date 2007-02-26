// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
//   Program:   FITS Reader for ITK
//   Module:    itkFITSImageIOFactory.cxx
//   Package: 	FITS IO
//   Author:    Douglas Alan <doug AT alum.mit.edu>
//              Initiative in Innovative Computing at Harvard University
//
//   Copyright (c) 2006 Douglas Alan
//
//   This software is freely distributable under the open source MIT X11
//   License.
//
//   See
//
//      http://www.opensource.org/licenses/mit-license
//
//   for details.
//
//=============================================================================

#include "itkFITSImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkFITSImageIO.h"
#include "itkVersion.h"
  
namespace itk
{

FITSImageIOFactory::FITSImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkFITSImageIO",
                         "FITS Image IO",
                         1,
                         CreateObjectFunction<FITSImageIO>::New());
}


FITSImageIOFactory::~FITSImageIOFactory()
{
  // Intentionally left blank.
}


const char* 
FITSImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}


const char* 
FITSImageIOFactory::GetDescription() const
{
  return "FITS ImageIO Factory, allows the loading of FITS images into Insight";
}

} // end namespace itk
