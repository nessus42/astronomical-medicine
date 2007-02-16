//=============================================================================
//
//   Program:   FITS Reader for ITK
//   Module:    itkFITSImageIOFactory.h
//   Language:  C++
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

#ifndef __itkFITSImageIOFactory_h
#define __itkFITSImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

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
  static void RegisterOneFactory(void)
  {
    FITSImageIOFactory::Pointer FITSFactory = FITSImageIOFactory::New();
    ObjectFactoryBase::RegisterFactory(FITSFactory);
  }

protected:
  FITSImageIOFactory();
  ~FITSImageIOFactory();

private:
  FITSImageIOFactory(const Self&); // Purposely not implemented.
  void operator=(const Self&);     // Purposely not implemented.

};
  
  
} // end namespace itk

#endif
