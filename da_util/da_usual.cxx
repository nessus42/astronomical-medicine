// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Descpription : Definitions that are used by neary everything written by
//                Douglas Alan.
// Author:        Douglas Alan <doug AT alum.mit.edu>
//
// Copyright (c) 1993-2007 Douglas Alan
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details
//=============================================================================

#include <errno.h>
#include <libgen.h>		// For basename()
#include <stdlib.h>             // For abort()
#include <string.h>
#include <iostream>
#include <string>

#include <pathToExecutable.h>

#include <da_usual.h>
#include <da_sugar.h>

using std::cerr;
using std::endl;
using std::string;

namespace da {

// File-scope variables:
local string _programName;

//-----------------------------------------------------------------------------
// runTimeError(): procedure
//-----------------------------------------------------------------------------

proc void proc
runTimeError(const string& message)
{
  cerr << endl;
  if (!programName().empty()) cerr << programName() << ": ";
  cerr << "FATAL ERROR: " << message << endl;
  exit(1);
}


//-----------------------------------------------------------------------------
// _internalError(): procedure used by the error() macro
//-----------------------------------------------------------------------------

proc void 
_internalError(const string& message, const char* filename,
	       const int line_number)
{
  cerr << endl;
  if (!programName().empty()) cerr << programName() << ": ";
  cerr << "FATAL INTERNAL ERROR: " << message << endl
       << "  File: " << filename << endl
       << "  Line: " << line_number << endl;
  abort();
}


//-----------------------------------------------------------------------------
// heapError(): procedure used by the checkHeap() procedure
//-----------------------------------------------------------------------------

proc void 
heapError()
{
  cerr << endl;
  if (!programName().empty()) cerr << programName() << ": ";
  cerr << "FATAL ERROR: Heap exhausted!" << endl;
  abort();
}


//-----------------------------------------------------------------------------
// _assertm(): procedure used by the assertm() macro
//-----------------------------------------------------------------------------

// TODO: get rid of the * const's here and in runTimeError()

proc void 
_assertm(const string& expr, const string& message,
	 const string& filename, const unsigned line_number)
{
  cerr << endl;
  if (!programName().empty()) cerr << programName() << ": ";
  cerr << "FATAL INTERNAL ERROR: Assertion Failed!" << endl
       << "  Error message: " << message << endl
       << "  Assertion: " << expr << endl
       << "  File: " << filename << endl
       << "  Line: " << line_number << endl;
  abort();
}


//-----------------------------------------------------------------------------
// badInvariant(): procedure
//-----------------------------------------------------------------------------

//# daBadInvariant() prints out an error message saying that the rep invariant
//# is bad and terminates the program, dumping core.

proc void 
badInvariant()
{
  error("Bad rep invariant!");
}


//-----------------------------------------------------------------------------
// breakpoint(): procedure
//-----------------------------------------------------------------------------

//# breakpoint() is a NOP for compiling in a breakpoint.

proc void 
breakpoint()
{}


//-----------------------------------------------------------------------------
// setProgramName(): procedure
//-----------------------------------------------------------------------------

proc void
setProgramName(const string& progname)
{
  assertm(_programName.empty(),
          "setProgramName() must not be invoked more than once");
  _programName = progname;
}


//-----------------------------------------------------------------------------
// programName(): procedure
//-----------------------------------------------------------------------------

proc string
programName()
{
  if (_programName.empty()) {
    _programName = ::basename(::pathToExecutable());
  }
  return _programName;
}


//-----------------------------------------------------------------------------
// warning(): procedure
//-----------------------------------------------------------------------------

proc void
warning(const string& warningMessage)
{
   if (!programName().empty()) cerr << programName() << ": ";
  cerr << "WARNING: " << warningMessage << endl;
}


//-----------------------------------------------------------------------------
// syscallWarning(): procedure
//-----------------------------------------------------------------------------

proc void
syscallWarning(const string& warningMessage)
{
  if (!programName().empty()) cerr << programName() << ": ";
  cerr << "WARNING: " << warningMessage << ": ";
  perror("");
}


//-----------------------------------------------------------------------------
// syscallError(): procedure
//-----------------------------------------------------------------------------

proc void
syscallError(const string& errorMessage)
{
  if (!programName().empty()) cerr << programName() << ": ";
  cerr << "FATAL ERROR: " << errorMessage;
  const char* const strerrorMessage = strerror(errno);
  if (strerrorMessage and errno != EINVAL) cerr << ": " << strerrorMessage;
  else cerr << ".";
  cerr << endl;
  exit(1);
}


//-----------------------------------------------------------------------------
// syscallErrorMessage(): procedure
//-----------------------------------------------------------------------------

proc const char*
syscallErrorMessage()
{
  const char* const strerrorMessage = strerror(errno);
  if (strerrorMessage and errno != EINVAL) return strerrorMessage;
  else return 0;
}


//-----------------------------------------------------------------------------
// raiseHeapExhaustedException(): procedure
//-----------------------------------------------------------------------------

proc void
raiseHeapExhaustedException()
{
  error("Heap exhausted!");
}

} // end namespace da
