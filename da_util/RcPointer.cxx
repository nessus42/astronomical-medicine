// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
//
// Description: An implementation of reference-counted pointers
// Module:      RcPointer.h
// Language:    C++
// Author :     Douglas Alan <doug AT alum.mit.edu>.  Plagiarized from an
//              article published in the March-April issue of *C++ Report*,
//              "Memory Management and Smart Pointers".
//
// Copyright (c) 1993 Douglas Alan
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details
//
//=============================================================================

#include <da_sugar.h>

void
rciObjectDeletionError()
{
  error("Must not delete an RciObject!");
}
