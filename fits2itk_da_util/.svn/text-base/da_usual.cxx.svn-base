// -*- Mode: C++; fill-column: 79; fill-prefix: "//# " -*-

static const char rcsid[] = 
"$Id: da_usual.C,v 1.2 2004/12/10 22:36:05 nessus Exp $";

// Description : Definitions that are used by everything written by Douglas
//               Alan.
// Author      : Douglas Alan <doug@alum.mit.edu>

#include <stdlib.h>             // For abort()
#include <stdio.h>
#include <iostream>
#include <string>

#include <da_usual.h>
#include <da_sugar.h>

using std::cerr;
using std::endl;
using std::string;

// File-scope variables:
local string programName;

//-----------------------------------------------------------------------------
// daRunTimeError(): procedure
//-----------------------------------------------------------------------------

proc void proc
daRunTimeError(const string& message)
{
  cerr << endl;
  if (!::programName.empty()) cerr << ::programName << ": ";
  cerr << "FATAL ERROR: " << message << endl;
  exit(1);
}


//-----------------------------------------------------------------------------
// _daInternalError(): procedure used by the error() macro
//-----------------------------------------------------------------------------

proc void 
_daInternalError(const string& message, const char* filename,
                 const int line_number)
{
  cerr << endl;
  if (!::programName.empty()) cerr << ::programName << ": ";
  cerr << "FATAL INTERNAL ERROR: " << message << endl
       << "  File: " << filename << endl
       << "  Line: " << line_number << endl;
  abort();
}


//-----------------------------------------------------------------------------
// daHeapError(): procedure used by the daCheckHeap() procedure
//-----------------------------------------------------------------------------

proc void 
daHeapError()
{
  cerr << endl;
  if (::programName.size()) cerr << ::programName << ": ";
  cerr << "FATAL ERROR: Heap exhausted!" << endl;
  abort();
}


//-----------------------------------------------------------------------------
// _daAssertm(): procedure used by the assertm() macro
//-----------------------------------------------------------------------------

// TODO: get rid of the * const's here and in daRunTimeError()

proc void 
_daAssertm(const string& expr, const string& message,
           const string& filename, const unsigned line_number)
{
  cerr << endl;
  if (::programName.size()) cerr << ::programName << ": ";
  cerr << "FATAL INTERNAL ERROR: Assertion Failed!" << endl
       << "  Error message: " << message << endl
       << "  Assertion: " << expr << endl
       << "  File: " << filename << endl
       << "  Line: " << line_number << endl;
  abort();
}


//-----------------------------------------------------------------------------
// daBadInvariant(): procedure
//-----------------------------------------------------------------------------

//# daBadInvariant() prints out an error message saying that the rep invariant
//# is bad and terminates the program, dumping core.

proc void 
daBadInvariant()
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
// daSetProgramName(): procedure
//-----------------------------------------------------------------------------

proc void
daSetProgramName(const string& progname)
{
  assertm(::programName.empty(),
          "daSetProgramName() must not be invoked more than once");
  ::programName = progname;
}


//-----------------------------------------------------------------------------
// daProgramName(): procedure
//-----------------------------------------------------------------------------

proc string
daProgramName()
{
  return ::programName;
}


//-----------------------------------------------------------------------------
// daWarning(): procedure
//-----------------------------------------------------------------------------

proc void
daWarning(const string& warningMessage)
{
   if (::programName.size()) cerr << ::programName << ": ";
  cerr << "WARNING: " << warningMessage << endl;
}


//-----------------------------------------------------------------------------
// daSyscallWarning(): procedure
//-----------------------------------------------------------------------------

proc void
daSyscallWarning(const string& warningMessage)
{
  if (::programName.size()) cerr << ::programName << ": ";
  cerr << "WARNING: " << warningMessage << ": ";
  perror("");
}


//-----------------------------------------------------------------------------
// daSyscallError(): procedure
//-----------------------------------------------------------------------------

proc void
daSyscallError(const string& errorMessage)
{
  if (::programName.size()) cerr << ::programName << ": ";
  cerr << "FATAL ERROR: " << errorMessage << ": ";
  perror("");
  exit(1);
}


//-----------------------------------------------------------------------------
// daRaiseHeapExhaustedException(): procedure
//-----------------------------------------------------------------------------

proc void
daRaiseHeapExhaustedException()
{
  error("Heap exhausted!");
}
