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

#include "itkFITSWCSTransform.h"
#include "da_sugar.h"

using std::cout;
using std::cerr;
using std::endl;

namespace itk
{

//*****************************************************************************
//*****                                                                   *****
//*****         FITSWCSTransform<double, 3>:                              *****
//*****              subclass of Transform<double, 3, 3 >                 *****
//*****                                                                   *****
//*****************************************************************************

//----------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------

ctor
FITSWCSTransform<double, 3>::FITSWCSTransform()
: Superclass(0, 0)
{}


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
  retval[0] = 5;
  retval[1] = 6;
  retval[2] = 55;
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
