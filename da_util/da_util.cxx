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

// For showFactoryClasses() //d
#include <itkFITSImageIOFactory.h>
int itk::FITSImageIOFactory::_cv_testValue = 0; //d
void //d
itk::FITSImageIOFactory::SetTestValue(int value)
{ _cv_testValue = value; } //d

int //d
itk::FITSImageIOFactory::GetTestValue()
{ return _cv_testValue; } //d





namespace douglasAlan {


ctor
IstreamExceptionActivator::
IstreamExceptionActivator(istream& stream)
  : _stream(stream), _originalState(stream.exceptions())
{
  _stream.exceptions(ios::badbit);
}

dtor
IstreamExceptionActivator::
~IstreamExceptionActivator()
{
  _stream.exceptions(_originalState);
}


ctor
OstreamExceptionActivator::
OstreamExceptionActivator(ostream& stream)
  : _stream(stream), _originalState(stream.exceptions())
{
  _stream.exceptions(ios::badbit | ios::eofbit);
}

dtor
OstreamExceptionActivator::
~OstreamExceptionActivator()
{
  _stream.exceptions(_originalState);
}


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


} // end namespace da
