// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  FITS Reader for ITK
// Module:   itkFITSImageIOFactory.cxx
// Package:  FITS IO
// Author:   Douglas Alan <douglas_alan AT harvard.edu>
//                        <doug AT alum.mit.edu>
//           Initiative in Innovative Computing at Harvard University
//
// Copyright (c) 2006-2007 Harvard University
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details.
//=============================================================================

#include "itkFITSImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkFITSImageIO.h"
#include "itkVersion.h"

#include <iostream> //d
using std::cerr;    //d 
using std::endl;    //d  

namespace itk
{

/*ctor*/
FITSImageIOFactory::FITSImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkFITSImageIO",
                         "FITS Image IO",
                         1,
                         CreateObjectFunction<FITSImageIO>::New());
}


/*dtor*/
FITSImageIOFactory::~FITSImageIOFactory()
{
  // Intentionally left blank.
}


/*method*/ const char* 
FITSImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}


/*method*/ const char* 
FITSImageIOFactory::GetDescription() const
{
  return "FITS ImageIO Factory, allows the loading of FITS images into ITK";
}

} // end namespace itk


// Function that is called when the shared library is loaded by
// itk::ObjectFactoryBase::LoadDynamicFactories():

/*proc*/ extern "C" itk::ObjectFactoryBase*
itkLoad()
{
  cerr << "Hello, I'm your friendly neighborhood itkLoad()" << endl; //d
  
  static itk::FITSImageIOFactory::Pointer f = itk::FITSImageIOFactory::New();
  return f;
}
