//=============================================================================
//
//   Program:   FITS Reader for ITK
//   Module:    itkFITSWCSTransform.cxx
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

#include <wcs.h>
#include "itkFITSWCSTransform.h"
#include "da_sugar.h"

using std::cout;
using std::cerr;
using std::endl;

namespace itk
{

//*****************************************************************************
//*****                                                                   *****
//*****      FITSWCSTransform<double, 3>:                                 *****
//*****         leaf subclass of Transform<double, 3, 3>                  *****
//*****                                                                   *****
//*****************************************************************************

//----------------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------------

dtor
FITSWCSTransform<double, 3>::~FITSWCSTransform()
{
  free(m_wcs);
}


//----------------------------------------------------------------------------
// Setter methods
//----------------------------------------------------------------------------

//! A FITSWCSTransform object takes possession of the WorldCoor object
//! that is passed to it.

method void
FITSWCSTransform<double, 3>::SetWCS(WorldCoor* wcs) {
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

  // YOU ARE HERE: Use libwcs to transform the coordinates.

  OutputPointType retval;
  pix2wcs(m_wcs, point[0], point[1], &retval[0], &retval[1]);
  retval[2] = 0;
  return retval;
}


//----------------------------------------------------------------------------
// GetInverse(): inherited non-virtual method
//----------------------------------------------------------------------------

method bool
FITSWCSTransform<double, 3>::GetInverse(Self* inverse) const
{
  if (!inverse) {
    return false;
  }
  itkExceptionMacro("FITSWCSTransform<...>::GetInverse() is not yet "
		    "implemented");
}

} // namespace itk
