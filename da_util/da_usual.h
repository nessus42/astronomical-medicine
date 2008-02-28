// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Description : This file contains declarations that are needed in just about
//               every source file written by Douglas Alan.
// Author      : Douglas Alan <nessus@mit.edu>
//
// Copyright    (c) 1993-2007  Douglas Alan
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details
//=============================================================================

#ifndef DA_USUAL_H
#define DA_USUAL_H

#include <assert.h>
#include <iostream>
#include <string>

// BEGIN
namespace douglasAlan {

//=============================================================================
// Global constants
//=============================================================================

enum SuccessFlag { succeeded, failed };


//=============================================================================
// Global functions
//=============================================================================

// Global functions:
void            breakpoint();
void            runTimeError(const std::string& message);
void            heapError();
std::string     programName();
void            raiseHeapExhaustedException();
void		setDebugLevel(int debugLevel);
void		setVerbosityLevel(int verbosityLevel);
void            syscallWarning(const std::string& warningMessage);
void            syscallError(const std::string& errorMessage);
const char*     syscallErrorMessage();
void            warning(const std::string& warning);


//-----------------------------------------------------------------------------
// Internal functions
//-----------------------------------------------------------------------------

namespace internal {

  void          internalError(const std::string& message,
			      const char* filename,
			      int line_number);
  void          assertm_(const std::string&, const std::string&,
			 const std::string&, unsigned);
}


//-----------------------------------------------------------------------------
// Internal variables
//-----------------------------------------------------------------------------

namespace internal {

  extern int debugLevel;
  extern int verbosityLevel;

}


//=============================================================================
// Freer: concrete data type
//=============================================================================

//! Frees a malloced object via RAII.

class Freer {

  void* _malloced_object;

  // Disable copying:
  Freer(const Freer&);           
  void operator=(const Freer&);

public:
  Freer(void* malloced_object) : _malloced_object(malloced_object) {}
  ~Freer() { free((char *)_malloced_object); }
};


//=============================================================================
// nop(): inline function
//=============================================================================

//# nop() does absolutely nothing.

/*proc*/ inline void
nop()
{}


//=============================================================================
// getDebugLevel(): inline function
//=============================================================================

/*proc*/ inline int
getDebugLevel()
{
  return internal::debugLevel;
}


//=============================================================================
// getDebugLevel(): inline function
//=============================================================================

/*proc*/ inline int
getVerbosityLevel()
{
  return internal::verbosityLevel;
}


//=============================================================================
// debugPrint(): macro
//=============================================================================

#define daDebugPrint(message) \
   { if (douglasAlan::internal::debugLevel) \
        std::cerr << message << std::endl; }

//=============================================================================
// printIfVerbose(): local macro
//=============================================================================


#define daPrintIfVerbose(message) \
   { if (douglasAlan::internal::verbosityLevel) { \
        std::cout << message << '\n'; \
     } \
   }


//=============================================================================
// DA_NDEBUG: macro
//=============================================================================

//# If NDEBUG is defined, then we define DA_NDEBUG to be true, otherwise
//# we define DA_NDEBUG to be false.

#ifdef NDEBUG
#   define DA_NDEBUG 1
#else
#   define DA_NDEBUG 0
#endif // NDEBUG


//=============================================================================
// daIfDebug(): macro
//=============================================================================

//# This macro only evaluates the argument if NDEBUG is not defined.

#ifdef NDEBUG
#   define daIfDebug(ignore)
#else
#   define daIfDebug(expr) \
       expr
#endif // NDEBUG


//=============================================================================
// daSuperDebug(): macro
//=============================================================================

//# This macro only evaluates the argument if SUPER_DEBUG is defined.

#ifdef SUPER_DEBUG
#   define daSuperDebug(expr) \
       expr
#else
#   define daSuperDebug(ignore)
#endif // SUPER_DEBUG


//=============================================================================
// daCheckHeap(): inline template function
//=============================================================================

//# This template function is kind of like the 'assert' macro, except that the
//# argument is always evaluated, the result of evaluating the argument is
//# returned, and if the result of evaluating the argument is false, then
//# the error message that is output indicates that the heap is exhausted.

template <class T> inline T
checkHeap(T arg)
{
  if (DA_NDEBUG && !arg) heapError();
  return arg;
}


//=============================================================================
// daAssertm(): macro
//=============================================================================

//# daAssertm() is just like the 'assert' macro, except that it takes an extra
//# argument, which is an error message to output if the assertion fails.

#ifdef NDEBUG
#   define daAssertm(ignore1, ignore2)
#else
#   define daAssertm(expr, message) \
        if (expr); else douglasAlan::internal::assertm_(#expr, message, \
					                __FILE__, __LINE__)
#endif // NDEBUG


//=============================================================================
// daError(): macro
//=============================================================================

//# daError() generates a fatal error, writing the string 'message' to 'cerr'
//# to let the user have some idea of what went wrong.

#define daError(message) \
        douglasAlan::internal::internalError(message, __FILE__, __LINE__)


//=============================================================================
// daCheckSyscall(): macro
//=============================================================================

#define daCheckSyscall(expr, warningMessage) \
   if ((expr) != succeeded) { \
       douglasAlan::syscallWarning(warningMessage) \
       } else douglasAlan::nop()


//=============================================================================
// daCheckSyscallE(): macro
//=============================================================================

#define daCheckSyscallE(expr, errorMessage) \
   if ((expr) != succeeded) douglasAlan::syscallError(errorMessage); \
   else douglasAlan::nop()


} // END namespace douglasAlan

#endif // DA_USUAL_H
