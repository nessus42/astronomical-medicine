// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  FITS Reader for ITK
// Module:   itkFITSWCSTransform.cxx
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

#include <wcs.h>
#include <itkFITSWCSTransform.h>
#include <itkFITSImageIO.h>
#include <da_sugar.h>

namespace itk
{

const int c_dims = FITSImageIO::c_dims;


//*****************************************************************************
//*****                                                                   *****
//*****      FITSWCSTransform<double, c_dims>:                            *****
//*****         leaf subclass of Transform<double, c_dims, c_dims>        *****
//*****                                                                   *****
//*****************************************************************************

//----------------------------------------------------------------------------
// Constructors, etc.
//----------------------------------------------------------------------------

method FITSWCSTransform<double, c_dims>&
FITSWCSTransform<double, c_dims>::
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
FITSWCSTransform<double, c_dims>::
SetWCS(const ConstRcMallocPointer<WorldCoor>& wcs) {
  assert(!m_wcs);
  m_wcs = wcs;
}


//----------------------------------------------------------------------------
// PrintSelf(): inherited virtual method 
//----------------------------------------------------------------------------

method void
FITSWCSTransform<double, c_dims>::
PrintSelf(std::ostream &os, Indent indent) const
{
  // TODO: Make this real.
  Superclass::PrintSelf(os, indent);
  os << indent << "Temporarily left blank." << endl;
}


//----------------------------------------------------------------------------
// TransformPoint(): inherited virtual method
//----------------------------------------------------------------------------

method FITSWCSTransform<double, c_dims>::OutputPointType
FITSWCSTransform<double, c_dims>::
TransformPoint(const InputPointType &point) const
{
  assert(m_wcs);

  OutputPointType retval;
   
  // Don't worry -- we are not going to modify m_wcs.  We just need to
  // cast away const so that pix2wcs will accept it: 
  WorldCoor* wcs = const_cast<WorldCoor*>(m_wcs.rawPointer());

  if (m_isInverted) {
    int offscl;
    wcs2pix(wcs, point[0], point[1], &retval[0], &retval[1], &offscl);

    // Adjust for the fact that the index origin for FITS is (1, 1, 1) and the
    // index origin for ITK is (0, 0, 0):
    retval[0] -= 1;
    retval[1] -= 1;

    // We don't do WCS for the velocity axis, so for now, we just set it to 0
    // (TODO: Fix this):
    retval[2] = 0;

  } else {
    // The +1's below are to compensate for the fact that the index origin for
    // FITS is (1, 1, 1) and the index origin for ITK is (0, 0, 0):
    pix2wcs(wcs, point[0] + 1, point[1] + 1, &retval[0], &retval[1]);

    // We don't do WCS for the velocity axis, so for now, we just set it to 0
    // (TODO: Fix this):
    retval[2] = 0;
  }
  return retval;
}


//----------------------------------------------------------------------------
// GetInverse(): inherited non-virtual method
//----------------------------------------------------------------------------

method bool
FITSWCSTransform<double, c_dims>::GetInverse(Self* inverse) const
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
FITSWCSTransform<double, c_dims>::Invert()
{
  m_isInverted = !m_isInverted;
}

} // namespace itk
