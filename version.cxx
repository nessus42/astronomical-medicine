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

extern const char fits2itkVersion[] = "0.4dev.13";

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

//---------------------------
// Version 0.4dev.2pending.2
//---------------------------

// Added "--scale-dec" option to allow the user to scale the declination axis.

//----------------------------------------------------
// Version 0.4dev.2pending.2make-analyze-output-work.0
//----------------------------------------------------

// *** Thu Jul 26, 2007 ***

// Added an --RIP option to fits2itk that forces the image to be in the "RIP"
// spatial orientation.  (See the Insight Journal paper, "A Definition of
// Spatial Orientation for the Insight Toolkit".)  The way that this option
// works is set the change of basis matrix accordingly.

// The purpose of this option is to allow Analyze file output without getting a
// heinous error message.  Unfortunately, it also causes the Velocity axis to
// be reversed, and there would seem to be no way around this without creating
// a new image with all the pixels moved around.  FlipImageFilter doesn't help
// any because apparently it just manipulates the direction cosines, or
// something.  Whatever the case, the Analyze writer knows somehow that the
// orientation has been flipped, and then the error message returns.

//---------------------------
// Version 0.4dev.3pending
//---------------------------

// *** Mon Aug  6, 2007 ***

// Working on the "--flipv" option, to flip the velocity axis without mucking
// with direction cosines.  I.e., it will work by actually moving pixels around
// in the image.  The purpose for this is so that Analyze files can be properly
// displayed in OsiriX, et. al.

//---------------------------
// Version 0.4dev.3flipv-any
//---------------------------

// *** Tue Aug  7, 2007 ***

// Modified the flipv() function to allow flipping of any axis.  Boy this was
// really hairy to do it all in place in only one pass!

//----------------------------
// Version 0.4dev.3flipv-any.1
//----------------------------

// *** Wed Aug  8, 2007 ***

// Reordered the loop index ordering so that the ram access will be more
// sequential than striding, as it was previously.  Hopefully this will speed
// things up a bit.

//---------------------------
// Version 0.4dev.4
//---------------------------

// *** Sat Aug 11, 2007 ***

// Changed isOdd() to be a bit more elegant.

//---------------------------
// Version 0.4dev.5
//---------------------------

// *** Mon Aug 13, 2007 ***

// Added support for the "-n" option, which causes fits2itk to not write any
// output.

//---------------------------
// Version 0.4dev.6
//---------------------------

// *** Mon Aug 13, 2007 ***

// Implemented the "--verbose" option, which current prints out the IJK center
// coordinates and the RA/Dec center coordinates.  In the future it will do
// more.

//---------------------------
// Version 0.4dev.7
//---------------------------

// *** Tue Aug 14, 2007 ***

// Added to the "--verbose" output, information on the "i" and "j" basis
// vectors tranformed into RA/Dec space.  Also it now outputs these basis
// vectors transformed into "approximate angular" space.  Or that's what I'm
// calling it in the meantime.

//---------------------------
// Version 0.4dev.8
//---------------------------

// *** Tue Aug 14, 2007 ***

// Cosmetic improvements to "--verbose" output.  Removed some extraneous "-D1"
// debugging output.

// *** Mon Aug 20, 2007 ***

// Fixed a bug where the WCS information being returned was being shifted down
// and to the left by one pixel due to the disparity between ITK's (0, 0)
// origin and FITS's (1, 1) origin.

//---------------------------
// Version 0.4dev.9pending.0
//---------------------------

// *** Mon Aug 20, 2007 ***

// Refactored code to move calculation code out of writeImageInfo() and into
// its own class, ImageInfo<ImageType>.

//---------------------------
// Version 0.4dev.9pending.1
//---------------------------

// *** Tue Aug 21, 2007 ***

// Added angular length of unit i and j vectors to "--verbose" output.

// Changed some code to use ITK Matrix objects and also implemented
// "--reorient-north".

//---------------------------
// Version 0.4dev.9pending.2
//---------------------------

// Figured out everything there is to know about coordinate frame
// transformations and used this knowledge to fix the fact that the rotation
// was happening in the wrong direction, and rotating into the velocity axis.
// Also found and fixed a bug in which the image was being distorted due to
// initializeCoordinateFrame() not calling image.SetDirection().

//---------------------------
// Version 0.4dev.9
//---------------------------

// *** Wed Aug 22, 2007 ***

// Added the "--equiangular" option.

// Also did some code refactoring.

//---------------------------
// Version 0.4dev.10
//---------------------------

// *** Wed Aug 29, 2007 ***

// Fixed the coordinate frame transformations associated with "--equiangular"
// and "--reorient-north" after spending a week coming to grok them more fully.

// *** Fri Aug 31, 2007 ***

// Changed the names of some variables and methods to make more sense.

//---------------------------
// Version 0.4dev.11
//---------------------------

// *** Fri Sep  7, 2007 ***

// Modified CMakeLists.txt so that itk::FITSImageIO is a dynamic library.

// Also, merged in the changes to CMakeLists.txt that I made during NAMIC
// Programming Week.  These changes make the usage message files depend on the
// fits2itk executable, rather than on version.o

//---------------------------
// Version 0.4dev.12pending
//---------------------------

// Added support for loading FITSImageIO as a "dynamically loaded library",
// rather than as a "dependent shared library" or as a statically linked
// library.  I tried out this support with Slicer and it worked (at least
// nominally) but I did not yet succeed in making fits2itk dynamically load
// FITSImageIO.

//---------------------------
// Version 0.4dev.12pending.1
//---------------------------

// *** Tue Nov 13, 2007 ***

// With this version I made FITSImageIO be a "dynamically loaded library".
// rather than a dependent library.  I.e., now FITSImageIO is loaded at
// runtime, under the control of fits2itk.

// (Well, actually not directly under the control of fits2itk, because, as it
// turns out, the only way to get this to work was to set the environment
// variable $ITK_AUTOLOAD_PATH so that ITK could find libitkFITSImageIO.so as a
// plugin.  Other approaches I tried would either fail to register
// FITSImageIOFactory in the real ITK universe, as it either would end up
// registered in a hidden parallel dynamically-loaded namespace, or it would
// cause core dumps due to data structures in ObjectFactoryBase not being
// initialized properly.)

// This version is really just a first stab at all this, as the code is
// currently a mess and has lots of stuff commented out (using "//d") in order
// to make it work.  The reason for this is that many of the fits2itk command
// line options worked by calling static methods of FITSImageIO, but now that
// FITSImageIO is dynamically loaded, these methods are not easily accessible.

//---------------------------
// Version 0.4dev.12pending.2
//---------------------------

// *** Tue Nov 13, 2007 ***

// Uncommented out stuff in fits2itk that depended only on the FITSImageIO
// header and not on any of FITSImageIO's code.

//---------------------------
// Version 0.4dev.12pending.3
//---------------------------

// *** Thu Nov 15, 2007 ***

// Moved all of the image manipulating stuff out of fits2itk.cxx
// into itkFITSImageUtil.{cxx,txx,h}.

//---------------------------
// Version 0.4dev.12pending.4
//---------------------------

// *** Tue Nov 27, 2007 ***

// This version is a check in of what I had done before starting work on the
// followng version.

//---------------------------
// Version 0.4dev.12pending.5
//---------------------------

// *** Thu Nov 29, 2007 ***

// FITSImageIO is now dynamically loaded in main() and then also later loaded
// automatically by ITK's ObjectFactoryBase.  As hoped, both dlopen()'s end up
// getting the same version of FITSImageIO, and they share the same static
// state.  So, to set state in FITSImageIO before it is used by the IO factory
// for reading files, main() calls C-linkage setters that are exported in
// itkFITSImageIO.cxx.

//---------------------------
// Version 0.4dev.12pending.6
//---------------------------

// *** Thu Dec  6, 2007 ***

// Commented out almost all of the deprecated code from FITSImageIO, and got
// everything to compile.  (I don't know whether it actually works yet, but it
// doesn't dump core.)  This involved moving constructing the WCS transform
// object out of FITSImageIO, which isn't where it belonged anyway.

// Cleaned up CMakeLists.txt.

//---------------------------
// Version 0.4dev.12pending.7
//---------------------------

// Went through all of the code annotated with "//d" as a first pass at making
// things sane again.

// Commented out most functions and variables with "deprecated" in their names.

// Made the default behavior be "--equiangular --reorient-north" so that you
// don't have to specify these options anymore.  With all this changing of
// things around, there were some bugs related to ITK's smart pointers causing
// core dumps, but I> fixed them.

//---------------------------
// Version 0.4dev.12pending.8
//---------------------------

// *** Thu Dec 13, 2007 ***

// Changed ImageInfo into FITSImage, so that I don't have to keep making the
// ImageInfo object over and over again.

// Also moved all of the coordinate transformation stuff into pure matrix
// arithmetic stuff, rather than having strange kludges.

// Status: Everything compiles and things sort of work, but there's still a
// bunch more straightening out to do, and reimplementation of lost features.

//---------------------------
// Version 0.4dev.12pending.9
//---------------------------

// Commented out unused code, etc.

//----------------------------
// Version 0.4dev.12pending.10
//----------------------------

// I worked on this version at NA-MIC 2008 Winter Project Week.

// Features added:

// The "--xml" option to support use of fits2itk as a Slicer3 CLI plugin.

//-----------------------------------------
// Version 0.4dev.12pending.10.gt.pending.0
//-----------------------------------------

// Made an initial stab at Mike's request for a program to generate grid
// transform approximation output.  I got it to compile and work, modulo
// velocities in the output file, which are all currently set to 0.  Also, I
// currently output a WCS pixel for each pixel in the original image.  First I
// need to change this to output only two slices, as the velocity information
// is completely linear (at least in the files that we are looking at), and
// then I need to modify it to output a WCS pixel for every nth pixel in the
// original image, where n is set on the command line.

// I examined the output data like so:

//    cat wcs-image.nrrd | unu save -f nrrd -e ascii | less

//-----------------------------------------
// Version 0.4dev.12pending.10.gt.pending.1
//-----------------------------------------

// In this version, we only output two slices into the WCS grid, and we
// decimate both the RA and Dec axes by factors of five (more or less).

// In the next version we should make the decimation factor be configurable on
// the command line.

//-----------------------------------------
// Version 0.4dev.12pending.10.gt.pending.2
//-----------------------------------------

// Fixed the problem wherein the velocity was always 0 in the WCS grid output.

//-----------------------------------------
// Version 0.4dev.12pending.10.gt
//-----------------------------------------

// Made the WCS grid image decimation factor be configurable on the command
// line.

//-----------------------------------------
// Version 0.4dev.12pending.10.gt.1
//-----------------------------------------

// Added the function mapRange(), which had the positive side-effect of fixing
// an insidious fencepost error.

//----------------------------------------------------
// Version 0.4dev.12pending.10.gt.2D-extrude.pending.0
//----------------------------------------------------

// Implemented the "fexstrude" script.

//----------------------------------------------------
// Version 0.4dev.12
//----------------------------------------------------

// The ultimate 0.4 release will embody a major effort in refactoring the code
// to be simpler, more elegant, and more maintainable.  The 0.4 release will
// also do its WCS linear interpolation based on just a few pixels at the
// center of the image.  This will allow it to be usable closer to the poles.

// This version is most of the way there.  It should just need a bit of
// polishing and a few bug fixes.

// I merged changes from fits2itk.fixes into this version.  I also
// reimplementing most of the features that were lost when I refactored the
// code.

// I revamped the command line options completely.

// I revamped the verbose output completely.

//----------------------------------------------------
// Version 0.4dev.13
//----------------------------------------------------

// Added the --tiff-output option.  I'm not sure that it works quite right,
// though.  Also, I think I need to come up with separate coercion options for
// input and output.


//----------------------------------------------------------------------
// *** Changes described above this line are checked in to Mercurial ***
//----------------------------------------------------------------------

// Things that I need to do soon:

// o Detailed notes about the current state of things is in NOTES.txt.

// o If we ever want the FITSImageIO DLL to actually work as a plugin with
//   something, then we'll have to pass it its settings (e.g., debug level and
//   null value) via an environment variable (or somesuch), rather than via
//   extern "C" calls.

// o Change CMake config to not dependent libraries being built by my CMake.
//   Instead do it the more normal way of telling CMake that my stuff depends
//   on other stuff outside of the CMake directory.  Then also write a build
//   script that will build everything.  (*Maybe* use scons to do this.)

// o Need to update the shell script wrapper that provides the GUI interface to
//   Slicer 3.

// o Need to rewrite the help output as it is now completely inaccurate.

// o Make grid image output work again if we ever need it.
