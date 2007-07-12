// -*- Mode: C++; fill-column: 79 -*-
//=============================================================================
// Program:  fits2itk
// Author:   Douglas Alan <douglas_alan AT harvard.edu>
//                        <doug AT alum.mit.edu>
//           Initiative in Innovative Computing at Harvard University
//
// Copyright (c) 2006-2007 Harvard University
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details.
//=============================================================================


extern const char fits2itkVersion[] = "Version 0.4dev.2pending.1";


//=============================================================================
// Change Log
//=============================================================================

// See http://wiki.iic.harvard.edu/AstroMed/Software_Version_Numbers for
// documentation on how these version numbers work.


//---------------------------
// Version 0.1dev.0
//---------------------------

// *** Tue Dec  5, 2006 ***

// Added version number to usage output and a change log comments section.

// Changed autoscaling of velocity axis to be average of raPerI and decPerI,
// rather than just raPerI.

// RA/DEC/V were getting mapped to L/P/S rather than to L/P/S.  Fixed this.

// Added option to turn on debugging messages and turned them off by default.
// They were all getting output all the time prior to this change.

// Note this version not checked into subversion due to repository being down.

//---------------------------
// Version 0.1dev.1
//---------------------------

// *** Thu Dec  7, 2006 ***

// Added option to scale all axes and removed hack that automatically scaled
// all axes by 1000.

// Implemented a better method of autoscaling the velocity axis by
// fitting it maximally into cube.

// Completely reworked the change-of-basis matrix stuff to be clearer
// and to actually work correctly.

// *** Fri Dec  8, 2006 ***

// Added option to scale RA axis.  The most important purpose for this option
// is the scale RA by -1 so that it will display properly in Slicer.

//---------------------------
// Version 0.1
//---------------------------

// *** Tue Dec 12, 2006 ***

// Changed version to "0.1".  Could not check it into the NA-MIC sandbox due to
// the corruption of their repository.  So started my own Mercurial repository.

//---------------------------
// Version 0.2dev.0
//---------------------------

// *** Wed Dec 31, 2007 ***

// itkFITSImageIO.cxx: Fixed bug that caused core dump when filename extention
// test was longer than actual filename.

// cfitsio/CMakeLists.txt: Modified to set the correct "-D" options to gcc.
// (At least for OS X.)

// libwcs/fitshead.h: Added a #define for "dec2str" and "str2dec" because
// functions were conflicting with functions of the same name somewhere in
// Carbon.

// Recompiled with ITK 3.0.1 so as to get turn off preemptive denial of reading
// non-existant files.  (This is needed to support cfitsio's extended file
// syntax, which allows you to specify "filenames" that indicate more than just
// the filename.

// *** Thu Feb  1, 2007 ***

// itkFITSImageIO.cxx: Modified checkExtension() to use a complicated regular
// expression, rather than just doing a string compare on the end of the
// string.

//---------------------------
// Version 0.2
//---------------------------

// *** Mon Feb  5, 2007 ***

// Added a command line option to set the FITS null value for reading FITS
// files.

// *** Thu Feb  8, 2007 ***

// Added support for debug level -1, which is for testing out my
// itk::Transform.  The transform has been copied from
// itk::TranslationTransform, but I haven't done anything to it yet.  I've just
// made sure that it compiles.

// *** Tue Feb 13, 2007 ***

// Made itkFITSWCSTransform.h and itkFITSWCSTransform.txx and got them
// to compile.  They don't do anything interesting yet.

// Prettied up all the FITS Reader source files.

// *** Thu Feb 15, 2007 ***

// Added note to usage string telling how to slice off a useless forth
// dimension.

//---------------------------
// Version 0.3dev.0
//---------------------------

// *** Tue Feb 20, 2007 ***

// Added WorldCoor to FITSWCSTransform.

// *** Fri Feb 23, 2007

// Got both forward and reverse transforms to work for FITSWCSTransform, and
// made the debugging output prove this fact.

//---------------------------
// Version 0.3dev.1pending
//---------------------------

// *** Mon Feb 26, 2007 ***

// Modified to optionally feed the image through an ITK FlipImageFilter and see
// if the metadata dictionary is preserved.  Alas, it wasn't.  Also, the
// resulting nrrd output crashes Slicer, and if I try to output an Analyze file
// instead, I get an error about the image not being in RPI, PIR, or RIP
// orientation.

//---------------------------
// Version 0.3dev.1
//---------------------------

// *** Mon Mar  5, 2007 ***

// Added an option to apply the ITK BinomialBlurImageFilter.

// Reworked the code that applies filters so that you can apply multiple
// filters.

// Moved command line parsing into a class.

//---------------------------
// Version 0.3
//---------------------------

// *** Mon Mar 19, 2007 ***

// Deleted some old commented-out code, and cleaned up a few things.

// Added support for long options.

// Added a much longer usage message when "--help" is specified.

//---------------------------
// Version 0.3.1
//---------------------------

// *** Fri Apr 13, 2007 ***

// Revamped CMake configuration.

// Put bundled CFITSIO library back under the control of autoconf.

// Changed the names of the library subdirs.

// Updated the file headers and copyright and licensing info.

//---------------------------
// Version 0.4dev.0
//---------------------------

// *** Fri Jun 15, 2007 ***

// Coded up pathToExecutable library, and made sure that it works.

// Made a pathToExecutable subdir for this library and configured CMake to
// build it.

//---------------------------
// Version 0.4dev.1
//---------------------------

// *** Fri Jun 22, 2007 ***

// Moved usage text out of source code and into text files that are copied to
// stdout or stderr at runtime.  This way don't have to maintain them as
// cumbersome strings within the program.

// Modified the CMake config so as to automatically copy the usage message
// files from the source code directory into the build directory, so that
// during development fits2itk can locate them at runtime.

// Moved some personal library routines into the "da" namespace.

// Added a few utility routines to da_util.

// Added a facility for turing on IO exceptions when using streams that will
// automatically reset the stream (with regard to exception generation) back to
// the state it was in previously.

// Restructured the source code directory structure a bit.

//---------------------------
// Version 0.4dev.2pending
//---------------------------

// *** Fri Jul  6, 2007 ***

// Added --rotate-sky option to allow the user to specify a v-axis roll.

// *** Wed Jul 11, 2007 ***

// Added --no-wcs flag in order to allow pixel coordinates to be read off in
// Slicer, rather than RA and Dec.  Also causes images to not be rotated or
// squooshed at all due to the WCS information.

//---------------------------
// Version 0.4dev.2pending.1
//---------------------------

// *** Wed Jul 11, 2007 ***

// Refactored code a bit as support for the "--no-wcs" flag was rather ugly.

// Now the "--no-wcs" flag can be used in conjunction with other scaling
// options, which wouldn't work before.

//----------------------------------------------------------------------
// *** Changes described above this line are checked in to Mercurial ***
//----------------------------------------------------------------------

// TODO: Adding "--scale-dec" option to allow the user to scale the declination
// axis.
