// -*- Mode: C++; fill-column: 79 -*-
#ifndef da_pure_sugar_h
#define da_pure_sugar_h
//=============================================================================
// Author      : Douglas Alan <nessus@mit.edu>
//
// Copyright    (c) 1993-2007  Douglas Alan
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details
//=============================================================================

//# 'proc', 'method', 'ctor', and 'dtor' are defined to be the null
//# string.  They are used for declaring procedures and methods to make them
//# easier to find with grep, etc:

#include <iostream>

#define proc
#define method
#define private_method
#define ctor
#define dtor
#define local static
#define null  NULL

//# Aliases I use in my .cxx files to make my code prettier:

#define assertm         daAssertm
#define warning         daWarning
#define error           daError
#define ifDebug         daIfDebug
#define superDebug      daSuperDebug

using da::checkHeap;
using da::runTimeError;

using std::cout;
using std::cerr;
using std::endl;

// Define operator!=() for every case in which "==" is defined:

// template <class T1, class T2> inline bool
// operator!=(const T1& operand1, const T2& operand2)
// {
//   return bool(!(operand1 == operand2));
// }

#endif // DA_PURE_SUGAR_H
