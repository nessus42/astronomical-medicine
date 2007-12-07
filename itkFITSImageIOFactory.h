// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  FITS Reader for ITK
// Module:   itkFITSImageIOFactory.h
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

#ifndef _itkFITSImageIOFactory_h
#define _itkFITSImageIOFactory_h

#include <itkObjectFactoryBase.h>
#include <itkImageIOBase.h>

namespace itk
{

//! Create instances of FITSImageIO objects using an object factory.

class ITK_EXPORT FITSImageIOFactory : public ObjectFactoryBase
{
public:  

  // Standard class typedefs:
  typedef FITSImageIOFactory   Self;
  typedef ObjectFactoryBase  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  // Class methods used to interface with the registered factories:
  virtual const char* GetITKSourceVersion(void) const;
  virtual const char* GetDescription(void) const;
  
  // Method for class instantiation:
  itkFactorylessNewMacro(Self);

  // Run-time type information (and related methods):
  itkTypeMacro(FITSImageIOFactory, ObjectFactoryBase);

  // Register one factory of this type:
  // static void RegisterOneFactory();

  // Note: We no longer need the RegisterOneFactory() class method, as we are
  // now having ITK autoload FITSImageIO as a plugin, and the plugin mechanism
  // looks for the function itkLoad() in libitkFITSImageIO.so instead.
  // itkLoad() is defined in itkFITSImageIOFactory.cxx.

protected:
  FITSImageIOFactory();
  ~FITSImageIOFactory();

private:
  FITSImageIOFactory(const Self&); // Purposely not implemented.
  void operator=(const Self&);     // Purposely not implemented.

};
  
  
} // end namespace itk

#endif
