// -*- Mode: C++; fill-column: 79; fill-prefix: "//# " -*-
// $Id: da_usual.h,v 1.2 2004/12/10 22:36:05 nessus Exp $

// Description : This file contains declarations that are needed in just about
//               every source file written by Douglas Alan.
// Author      : Douglas Alan <nessus@mit.edu>
// Copyright    (C) 1993-2005  Douglas Alan


#ifndef DA_USUAL_H
#define DA_USUAL_H

#include <assert.h>
#include <string>

// Global procedures:

void            daRunTimeError(const std::string& message);
void            _daInternalError(const std::string& message,
                                 const char* filename,
                                 int line_number);
void            _daAssertm(const std::string&, const std::string&,
                           const std::string&, unsigned);
void            daBreakpoint();
void            daHeapError();
void            daRaiseHeapExhaustedException();
void            daSetProgramName(const std::string& progname);
void            daSyscallWarning(const std::string& warningMessage);
void            daSyscallError(const std::string& errorMessage);
std::string     daProgramName();
void            daWarning(const std::string& warning);


//-----------------------------------------------------------------------------
// daNop(): global inline procedure
//-----------------------------------------------------------------------------

//# daNop() does absolutely nothing.

inline void
daNop()
{}

//-----------------------------------------------------------------------------
// DA_NDEBUG: macro
//-----------------------------------------------------------------------------

//# If NDEBUG is defined, then we define DA_NDEBUG to be true, otherwise
//# we define DA_NDEBUG to be false.

#ifdef NDEBUG
#   define DA_NDEBUG 1
#else
#   define DA_NDEBUG 0
#endif // NDEBUG


//-----------------------------------------------------------------------------
// daIfDebug(): macro
//-----------------------------------------------------------------------------

//# This macro only evaluates the argument if NDEBUG is not defined.

#ifdef NDEBUG
#   define daIfDebug(ignore)
#else
#   define daIfDebug(expr) \
       expr
#endif // NDEBUG


//-----------------------------------------------------------------------------
// daSuperDebug(): macro
//-----------------------------------------------------------------------------

//# This macro only evaluates the argument if SUPER_DEBUG is defined.

#ifdef SUPER_DEBUG
#   define daSuperDebug(expr) \
       expr
#else
#   define daSuperDebug(ignore)
#endif // SUPER_DEBUG


//-----------------------------------------------------------------------------
// daCheckHeap(): global inline template procedure
//-----------------------------------------------------------------------------

//# This template procedure is kind of like the 'assert' macro, except that the
//# argument is always evaluated, the result of evaluating the argument is
//# returned, and if the result of evaluating the argument is false, then
//# the error message that is output indicates that the heap is exhausted.

template <class T> inline T
daCheckHeap(T arg)
{
  if (DA_NDEBUG && !arg) daHeapError();
  return arg;
}


//-----------------------------------------------------------------------------
// daAssertm(): macro
//-----------------------------------------------------------------------------

//# daAssertm() is just like the 'assert' macro, except that it takes an extra
//# argument, which is an error message to output if the assertion fails.

#ifdef NDEBUG
#   define daAssertm(ignore1, ignore2)
#else
#   define daAssertm(expr, message) \
        if (expr); else _daAssertm(#expr, message, __FILE__, __LINE__)
#endif // NDEBUG


//-----------------------------------------------------------------------------
// daError(): macro
//-----------------------------------------------------------------------------

//# daError() generates a fatal error, writing the string 'message' to 'cerr'
//# to let the user have some idea of what went wrong.

#define daError(message) \
        _daInternalError(message, __FILE__, __LINE__)


//-----------------------------------------------------------------------------
// daCheckSyscall(): macro
//-----------------------------------------------------------------------------

#define daCheckSyscall(expr, warningMessage) \
   if ((expr) != succeeded) daSyscallWarning(warningMessage); else daNop()

//-----------------------------------------------------------------------------
// daCheckSyscallE(): macro
//-----------------------------------------------------------------------------

#define daCheckSyscallE(expr, errorMessage) \
   if ((expr) != succeeded) daSyscallError(errorMessage); else daNop()

#endif // DA_USUAL_H
