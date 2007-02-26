// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
//   Program:   FITS Reader for ITK
//   Module:    itkFITSWCSTransform.cxx
//   Package: 	FITS IO
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

#include <wcs.h>
#include <itkFITSWCSTransform.h>
#include <da_sugar.h>

namespace itk
{

//*****************************************************************************
//*****                                                                   *****
//*****      FITSWCSTransform<double, 3>:                                 *****
//*****         leaf subclass of Transform<double, 3, 3>                  *****
//*****                                                                   *****
//*****************************************************************************

//----------------------------------------------------------------------------
// Constructors, etc.
//----------------------------------------------------------------------------

method FITSWCSTransform<double, 3>&
FITSWCSTransform<double, 3>::
operator=(const Self& orig)
{
  if (this != &orig) {
    m_wcs = orig.m_wcs;
    m_isInverted = orig.m_isInverted;
  }
  return *this;
}


//----------------------------------------------------------------------------
// Setter methods
//----------------------------------------------------------------------------

method void
FITSWCSTransform<double, 3>::
SetWCS(const ConstRcMallocPointer<WorldCoor>& wcs) {
  assert(!m_wcs);
  m_wcs = wcs;
}


//----------------------------------------------------------------------------
// PrintSelf(): inherited virtual method 
//----------------------------------------------------------------------------

method void
FITSWCSTransform<double, 3>::
PrintSelf(std::ostream &os, Indent indent) const
{
  // TODO: Make this real.
  Superclass::PrintSelf(os, indent);
  os << indent << "Temporarily left blank." << endl;
}


//----------------------------------------------------------------------------
// TransformPoint(): inherited virtual method
//----------------------------------------------------------------------------

method FITSWCSTransform<double, 3>::OutputPointType
FITSWCSTransform<double, 3>::TransformPoint(const InputPointType &point) const
{

  OutputPointType retval;
   
  // Don't worry -- we are not going to modify m_wcs.  We just need to
  // cast away const so that pix2wcs will accept it: 
  WorldCoor* wcs = const_cast<WorldCoor*>(m_wcs.rawPointer());

  if (m_isInverted) {
    int offscl;
    wcs2pix(wcs, point[0], point[1], &retval[0], &retval[1], &offscl);
    retval[2] = 0;
  } else {
    pix2wcs(wcs, point[0], point[1], &retval[0], &retval[1]);
    retval[2] = 0;
  }
  return retval;
}


//----------------------------------------------------------------------------
// GetInverse(): inherited non-virtual method
//----------------------------------------------------------------------------

method bool
FITSWCSTransform<double, 3>::GetInverse(Self* inverse) const
{
  assert(inverse);
  *inverse = *this;
  inverse->Invert();
  return true;
}


//----------------------------------------------------------------------------
// Invert(): non-virtual method
//----------------------------------------------------------------------------

method void
FITSWCSTransform<double, 3>::Invert()
{
  m_isInverted = !m_isInverted;
}

} // namespace itk
