// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  FITS Reader for ITK
// Module:   itkFITSWCSTransform.h
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

#ifndef __itkFITSWCSTransform_h
#define __itkFITSWCSTransform_h

#include <iostream>
#include <itkTransform.h>
#include <RcPointer.h>

struct WorldCoor;

namespace itk
{

  const int c_FITSWCSTransformNumOfDimensions = 3;

//*****************************************************************************
//*****                                                                   *****
//*****         FITSWCSTransform<T, N>:                                   *****
//*****              leaf template subclass of Transform<T, N, N>         *****
//*****                                                                   *****
//*****************************************************************************

//! @class FITSWCSTransform

//! FITS WCS (World Coordinate System) Transform

//! Transforms IJK image coordinates to astronomical physical coordinates
//! (i.e., right ascension [RA], declination [Dec], and velocity [V]), and
//! vice versa (via the back transform).

//! NOTE: At the moment, only FITSWCSTransform<double, c_dims> is implemented.

template<class TScalarType, unsigned int NDimensions>
class FITSWCSTransform {};

template<>
class ITK_EXPORT FITSWCSTransform<double, c_FITSWCSTransformNumOfDimensions>: 
    public Transform<double,
		     c_FITSWCSTransformNumOfDimensions,
		     c_FITSWCSTransformNumOfDimensions>
{
public:
  
  // The number of dimensions in an image:
  enum { c_dims = 3};

  // ----- Standard ITK class typedefs -----
  typedef FITSWCSTransform                      Self;
  typedef Transform<double, c_dims, c_dims>     Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  // ----- Import typedefs from superclass -----
  typedef Superclass::InputPointType            InputPointType;
  typedef Superclass::OutputPointType           OutputPointType;
      
  // ----- Standard ITK class declaration macro invocations -----
  itkNewMacro(Self);
  itkTypeMacro(FITSWCSTransform, Transform);

  // ----- Constructors, etc -----
  Self& 			operator=(const Self&);
 
  // ----- Setter methods -----
  void		  		SetWCS(const ConstRcMallocPointer<WorldCoor>&);

  // ----- Accessor methods -----
  ConstRcMallocPointer<WorldCoor>
				GetWCS() const { assert(m_wcs);
				                 return m_wcs; }

  // ----- Virtual methods inherited from Transform -----
  OutputPointType     	      	TransformPoint(const InputPointType  &point )
                                  const;
  bool        	      		IsLinear() const { return false; }

  // ----- Non-virtual methods -----
  bool                	      	GetInverse(Self* inverse) const;
  void				Invert();


protected:

  // ----- Protected constructors, etc. -----
  FITSWCSTransform(): Superclass(c_dims, 0), m_wcs(0), m_isInverted(false) {};

  // ----- Protected virtual methods inherited from LightObject -----
  void 				PrintSelf(std::ostream &os, Indent indent)
                                   const;

private:

  // Deactivate copy constructor:
  FITSWCSTransform(const Self&); // Intentionally not implemented.

  // ----- Instance variables -----
  ConstRcMallocPointer<WorldCoor>    m_wcs;
  bool				     m_isInverted;

}; // class FITSWCSTransform<double, c_dims>

}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_FITSWCSTransform(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT FITSWCSTransform< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef FITSWCSTransform< ITK_TEMPLATE_2 x > FITSWCSTransform##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkFITSWCSTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkFITSWCSTransform.txx"
#endif

#endif /* __itkFITSWCSTransform_h */
