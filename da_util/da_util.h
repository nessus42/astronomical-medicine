// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
//   Module:    da_util.h
//   Author:    Douglas Alan <douglas_alan AT harvard.edu>
//                           <doug AT alum.mit.edu>
//
//   Copyright (c) 2007 Douglas Alan
//
//   This software is freely distributable under the open source MIT X11
//   License.
//
//   See
//
//      http://www.opensource.org/licenses/mite-license
//
//   for details.
//
// ============================================================================

#ifndef __da_util_h
#define __da_util_h

#include <iostream>
#include <fstream>
#include <memory>
#include <da_usual.h>


namespace da {

//-----------------------------------------------------------------------------
// Global functions
//-----------------------------------------------------------------------------

enum ErrorAction { dieOnError, raiseExceptionOnError, returnCodeOnError};

SuccessFlag
copyStream(std::istream& srcStream, std::ostream& destStream,
	   ErrorAction errorAction=raiseExceptionOnError,
	   size_t bufferSize=8192)
  throw(std::ios::failure);


std::auto_ptr<std::ifstream>
openFileForReading(const char* filepath,
		   ErrorAction errorAction=raiseExceptionOnError)
  throw(std::ios::failure);


std::auto_ptr<std::ofstream>
openFileForWriting(const char* filepath,
		   ErrorAction errorAction=raiseExceptionOnError)
  throw(std::ios::failure);


//-----------------------------------------------------------------------------
// Classes
//-----------------------------------------------------------------------------

// IstreamExceptionActivator and OstreamExceptionActivator turn on exception
// throwing for, respectively, istreams and ostreams.  When an "exception
// activator" destructor is called, the stream is put back into whatever state
// is was in (with respect to exception throwing) prior to the construction of
// the exception activator.

class IstreamExceptionActivator
{
  // Instance variables:
  std::istream&         _stream;
  std::ios::iostate     _originalState;

public:

  IstreamExceptionActivator(std::istream&);
  ~IstreamExceptionActivator();

private:

  // Deactivate assignment and copying:
  void operator=(const IstreamExceptionActivator&);
  IstreamExceptionActivator(const IstreamExceptionActivator&);
};


class OstreamExceptionActivator
{
  // Instance variables:
  std::ostream&         _stream;
  std::ios::iostate     _originalState;

public:

  OstreamExceptionActivator(std::ostream&);
  ~OstreamExceptionActivator();

private:

  // Deactivate assignment and copying:
  void operator=(const OstreamExceptionActivator&);
  OstreamExceptionActivator(const OstreamExceptionActivator&);
};


} // end namespace da

#endif // __da_util_h
