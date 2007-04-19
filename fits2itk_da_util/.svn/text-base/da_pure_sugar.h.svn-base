// -*- Mode: C++; fill-column: 79; fill-prefix: "//# " -*-
#ifndef da_pure_sugar_h
#define da_pure_sugar_h

#define da_pure_sugar_rcsId_H \
"$Id: da_pure_sugar.h,v 1.4 2002/10/10 17:28:15 nessus Exp $"

// Description : [TBD]
// Author      : Douglas Alan <doug@alum.mit.edu>

//# 'proc', 'method', 'ctor', and 'dtor' are defined to be the null
//# string.  They are used for declaring procedures and methods to make them
//# easier to find with grep, etc:

#define proc
#define method
#define private_method
#define ctor
#define dtor
#define local static

//# Aliases I use in my .cxx files to make my code prettier:

#define assertm         daAssertm
#define badInvariant    daBadInvariant
#define checkHeap       daCheckHeap
#define error           daError
#define ifDebug         daIfDebug
#define rtError         daRunTimeError
#define superDebug      daSuperDebug
#define warning         daWarning

// Define operator!=() for every case in which "==" is defined:

// template <class T1, class T2> inline bool
// operator!=(const T1& operand1, const T2& operand2)
// {
//   return bool(!(operand1 == operand2));
// }

#endif // DA_PURE_SUGAR_H
