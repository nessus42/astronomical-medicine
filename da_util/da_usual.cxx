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

namespace douglasAlan {

//-----------------------------------------------------------------------------
// Internal variables
//-----------------------------------------------------------------------------

// File-scope variables:
local string f_programName;


namespace internal {
   int debugLevel = 0;
   int verbosityLevel = 0;
}


//-----------------------------------------------------------------------------
// runTimeError(): function
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
// internalError(): internal function used by the error() macro
//-----------------------------------------------------------------------------

namespace internal {

  proc void
  internalError(const string& message, const char* filename,
		const int line_number)
  {
    cerr << endl;
    if (!programName().empty()) cerr << programName() << ": ";
    cerr << "FATAL INTERNAL ERROR: " << message << endl
	 << "  File: " << filename << endl
	 << "  Line: " << line_number << endl;
    abort();
  }
}


//-----------------------------------------------------------------------------
// heapError(): function used by the checkHeap() function
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
// assertm(): internal function used by the assertm() macro
//-----------------------------------------------------------------------------

namespace internal {
proc void 
assertm_(const string& expr, const string& message,
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
} }


//-----------------------------------------------------------------------------
// badInvariant(): function
//-----------------------------------------------------------------------------

//# daBadInvariant() prints out an error message saying that the rep invariant
//# is bad and terminates the program, dumping core.

proc void 
badInvariant()
{
  error("Bad rep invariant!");
}


//-----------------------------------------------------------------------------
// breakpoint(): function
//-----------------------------------------------------------------------------

//# breakpoint() is a NOP for compiling in a breakpoint.

proc void 
breakpoint()
{}


//-----------------------------------------------------------------------------
// setProgramName(): function
//-----------------------------------------------------------------------------

proc void
setProgramName(const string& progName)
{
  assertm(f_programName.empty(),
          "setProgramName() must not be invoked more than once");
  f_programName = progName;
}


//-----------------------------------------------------------------------------
// programName(): function
//-----------------------------------------------------------------------------

proc string
programName()
{
  if (f_programName.empty()) {
    f_programName = basename(pathToExecutable());
  }
  return f_programName;
}

//-----------------------------------------------------------------------------
// setDebugLevel(): function
//-----------------------------------------------------------------------------

proc void
setDebugLevel(int debugLevel)
{
  assertm(!internal::debugLevel,
	  "setDebugLevel() must not be be invoked more than once");
  internal::debugLevel = debugLevel;
}


//-----------------------------------------------------------------------------
// setVerbosityLevel(): function
//-----------------------------------------------------------------------------

proc void
setVerbosityLevel(int verbosityLevel)
{
  assertm(!internal::verbosityLevel,
	  "setVerbosityLevel() must not be be invoked more than once");
  internal::verbosityLevel = verbosityLevel;
}


//-----------------------------------------------------------------------------
// warning(): function
//-----------------------------------------------------------------------------

proc void
warning(const string& warningMessage)
{
   if (!programName().empty()) cerr << programName() << ": ";
  cerr << "WARNING: " << warningMessage << endl;
}


//-----------------------------------------------------------------------------
// syscallWarning(): function
//-----------------------------------------------------------------------------

proc void
syscallWarning(const string& warningMessage)
{
  if (programName().empty()) cerr << programName() << ": ";
  cerr << "WARNING: " << warningMessage << ": ";
  perror("");
}


//-----------------------------------------------------------------------------
// syscallError(): function
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
// syscallErrorMessage(): function
//-----------------------------------------------------------------------------

proc const char*
syscallErrorMessage()
{
  const char* const strerrorMessage = strerror(errno);
  if (strerrorMessage and errno != EINVAL) return strerrorMessage;
  else return 0;
}


//-----------------------------------------------------------------------------
// raiseHeapExhaustedException(): function
//-----------------------------------------------------------------------------

proc void
raiseHeapExhaustedException()
{
  error("Heap exhausted!");
}

} // end namespace douglas_alan
