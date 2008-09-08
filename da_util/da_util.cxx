// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
//   Module:    da_util.cxx
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


#include <iostream>
using std::istream;
using std::ostream;
using std::ios;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <memory>
using std::auto_ptr;

#include <string>
using std::string;

#include <stdexcept>
#include <da_util.h>
#include <da_sugar.h>


// BEGIN
namespace douglasAlan {


//-----------------------------------------------------------------------------
// Constructor of IstreamExceptionActivator
//-----------------------------------------------------------------------------

ctor
IstreamExceptionActivator::
IstreamExceptionActivator(istream& stream)
  : _stream(stream), _originalState(stream.exceptions())
{
  _stream.exceptions(ios::badbit);
}


//-----------------------------------------------------------------------------
// Destructor of IstreamExceptionActivator
//-----------------------------------------------------------------------------

dtor
IstreamExceptionActivator::
~IstreamExceptionActivator()
{
  _stream.exceptions(_originalState);
}


//-----------------------------------------------------------------------------
// Constructor of OstreamExceptionActivator
//-----------------------------------------------------------------------------

ctor
OstreamExceptionActivator::
OstreamExceptionActivator(ostream& stream)
  : _stream(stream), _originalState(stream.exceptions())
{
  _stream.exceptions(ios::badbit | ios::eofbit);
}


//-----------------------------------------------------------------------------
// Destructor of OstreamExceptionActivator
//-----------------------------------------------------------------------------

dtor
OstreamExceptionActivator::
~OstreamExceptionActivator()
{
  _stream.exceptions(_originalState);
}


//-----------------------------------------------------------------------------
// copyStream(): function
//-----------------------------------------------------------------------------

proc SuccessFlag
copyStream(istream& srcStream, ostream& destStream,
	   const ErrorAction errorAction, const size_t bufferSize)
  throw(ios::failure)
{
  IstreamExceptionActivator srcActivator (srcStream); 
  OstreamExceptionActivator destActivator (destStream);
  char buffer[bufferSize];
  while (true) {
    try {
      srcStream.read(buffer, bufferSize);
      destStream.write(buffer, srcStream.gcount());
    }
    catch (ios::failure& e) {
      switch (errorAction) {
      case dieOnError:
	{ 
	  string errorMessage =
	    "Exception while copying help file to terminal: ";
	  errorMessage += e.what();
	  syscallError(errorMessage);
	}
      case raiseExceptionOnError: throw;
      case returnCodeOnError: break;
      default: error("Bug!");
      }
    } // catch
    if (srcStream.eof()) return succeeded;
  } // while
  error("Bug!");
}


//-----------------------------------------------------------------------------
// openFileForReading(): function
//-----------------------------------------------------------------------------

proc auto_ptr<ifstream>
openFileForReading(const char* const filepath,
		   const ErrorAction errorAction)
  throw(ios::failure)
{
  auto_ptr<ifstream> file (new ifstream(filepath));
  if (errorAction == returnCodeOnError) return file;
  else if (!*file) {
    string errorMessage = "Could not open file \"";
    errorMessage += filepath;
    errorMessage += "\" for reading";
    switch (errorAction) {
    case dieOnError:
      syscallError(errorMessage);
    case raiseExceptionOnError:
      const char* syserror = syscallErrorMessage();
      if (syserror) {
	errorMessage += ": ";
	errorMessage += syserror;
	errorMessage += ".";
      }
      throw ios::failure(errorMessage);
    default:
      error("Bug!");
    }
  } else return file;
}


//-----------------------------------------------------------------------------
// openFileForWriting(): function
//-----------------------------------------------------------------------------

proc auto_ptr<ofstream>
openFileForWriting(const char* const filepath,
		   const ErrorAction errorAction)
  throw(ios::failure)
{
  auto_ptr<ofstream> file (new ofstream(filepath));
  if (errorAction == returnCodeOnError) return file;
  else if (!*file) {
    string errorMessage = "Could not open file \"";
    errorMessage += filepath;
    errorMessage += "\" for writing";
    switch (errorAction) {
    case dieOnError:
      syscallError(errorMessage);
    case raiseExceptionOnError:
      const char* syserror = syscallErrorMessage();
      if (syserror) {
	errorMessage += ": ";
	errorMessage += syserror;
	errorMessage += ".";
      }
      throw ios::failure(errorMessage);
    default:
      error("Bug!");
    }
  } else return file;
}


//-----------------------------------------------------------------------------
// toLower(): function
//-----------------------------------------------------------------------------

// Converts \a s to lowercase.

proc void
toLower(string& s)
{
  for (string::iterator p = s.begin(); p != s.end( ); ++p) {
    *p = tolower(*p);
  }
}


//-----------------------------------------------------------------------------
// endMatchesP(): function
//-----------------------------------------------------------------------------

// Returns true iff \a extension is on the end of \a filepath.

proc bool
endMatchesP(const string& filepath, const string& extension)
{
  typedef string::size_type StrLen;
  const StrLen extensionLength = extension.length();
  const StrLen filepathLength = filepath.length();
  if (extensionLength >= filepathLength) return false;
  string filepathEnd = filepath.substr(filepathLength - extensionLength);
  if (filepathEnd == extension) return true;
  else return false;
}


} // END namespace douglasAlan
