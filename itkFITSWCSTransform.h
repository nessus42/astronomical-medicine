//=============================================================================
//
//   Program:   FITS Reader for ITK
//   Module:    itkFITSWCSTransform.h
//   Language:  C++
//   Author:    Douglas Alan <doug AT alum.mit.edu>
//              Initiative in Innovative Computing at Harvard University
//
//   Copyright (c) 2007 Douglas Alan
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

#ifndef __itkFITSWCSTransform_h
#define __itkFITSWCSTransform_h

#include <iostream>
#include "itkTransform.h"

struct WorldCoor;

namespace itk
{

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

//! NOTE: At the moment, only FITSWCSTransform<double, 3> is implemented.

template<class TScalarType, unsigned int NDimensions>
class FITSWCSTransform {};

template<>
class ITK_EXPORT FITSWCSTransform<double, 3>: 
    public Transform<double, 3, 3>
{
public:

  // Standard ITK class typedefs:
  typedef FITSWCSTransform                      Self;
  typedef Transform<double, 3, 3>               Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  // Import typedefs from superclass:
  typedef Superclass::InputPointType            InputPointType;
  typedef Superclass::OutputPointType           OutputPointType;
      
  // Standard ITK class declaration macro invocations:
  itkNewMacro(Self);
  itkTypeMacro(FITSWCSTransform, Transform);

  // Virtual methods inherited from Transform:
  OutputPointType     TransformPoint(const InputPointType  &point ) const;
  virtual bool        IsLinear() const { return false; }

  // Non-virtual methods inherited from Transform:
  bool                GetInverse(Self* inverse) const;

  //! Call when done setting or changing attributes.
  void		      Update() { assert(m_wcs); }

  // Setter methods:
  void		      SetWCS(WorldCoor*);

  // Accessor methods:
  const WorldCoor&    GetWCS() const { assert(m_wcs); return *m_wcs; }

protected:

  // Constructors, etc:
  FITSWCSTransform(): Superclass(3, 0), m_wcs(0) {};
  ~FITSWCSTransform();

  // Protected virtual methods inherited from LightObject:
  void PrintSelf(std::ostream &os, Indent indent) const;

private:

  // Deactivate copying:
  FITSWCSTransform(const Self&); // Intentionally not implemented.
  void operator=(const Self&);   // Intentionally not implemented.

  // Instance variables:
  WorldCoor*          m_wcs;


}; // class FITSWCSTransform<double, 3>

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
