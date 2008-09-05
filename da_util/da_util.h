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


namespace douglasAlan {

using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::auto_ptr;
using std::string;
  
//-----------------------------------------------------------------------------
// Functions
//-----------------------------------------------------------------------------

enum ErrorAction { dieOnError, raiseExceptionOnError, returnCodeOnError};

SuccessFlag
copyStream(istream& srcStream, ostream& destStream,
           ErrorAction errorAction=raiseExceptionOnError,
           size_t bufferSize=8192)
  throw(ios::failure);


auto_ptr<ifstream>
openFileForReading(const char* filepath,
		   ErrorAction errorAction=raiseExceptionOnError)
  throw(ios::failure);


auto_ptr<ofstream>
openFileForWriting(const char* filepath,
		   ErrorAction errorAction=raiseExceptionOnError)
  throw(ios::failure);

void toLower(string& s);
bool endMatchesP(const string& filepath, const string& extension);


inline double square(double x)
{
  return x * x;
}


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
  istream&         _stream;
  ios::iostate     _originalState;

public:

  IstreamExceptionActivator(istream&);
  ~IstreamExceptionActivator();

private:

  // Deactivate assignment and copying:
  void operator=(const IstreamExceptionActivator&);
  IstreamExceptionActivator(const IstreamExceptionActivator&);
};


class OstreamExceptionActivator
{
  // Instance variables:
  ostream&         _stream;
  ios::iostate     _originalState;

public:

  OstreamExceptionActivator(ostream&);
  ~OstreamExceptionActivator();

private:

  // Deactivate assignment and copying:
  void operator=(const OstreamExceptionActivator&);
  OstreamExceptionActivator(const OstreamExceptionActivator&);
};


} // end namespace douglasAlan

#endif // __da_util_h
