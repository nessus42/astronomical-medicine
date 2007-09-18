INTRODUCTION
------------

fits2itk is free software, whose source code is available under the
terms of the permissive "MIT License".  See LICENSE.txt for details.
It relies on libraries that are GPL'ed, however, so the executable is
licensed under the GPL.

fits2itk was written and is maintained by

      Douglas Alan <douglas_alan AT harvard.edu>
                   <doug AT alum.mit.edu>
      Initiative in Innovative Computing at Harvard University

Documentation and background on fits2itk, and the AstroMed project in
general, can be found at

      http://astromed.iic.harvard.edu


LICENSING ISSUES
----------------

The code that was written for fits2itk is being released under the
very permissive "MIT License", but it currently uses some GPL'ed
libraries.  As a result, fits2itk executables, as they are currently
being distributed (or as you are likely to build them, should you
chose to build the executables from source) acquire the requirements
of these copylefts.  At some point we may release a version of
fits2itk that does not carry copyleft requirements.

If you care about the gory details of the copyleft situation, read the
following: The copylefted libraries used by fits2itk are

   CFITSIO -- written by NASA and in the public domain, with the
      exception of a couple third-party subroutines that are GPL'ed or
      LGPL'ed.  CFITSIO can be compiled with these subroutines turned
      off, at the cost of some reduced functionality.  We have yet to
      attempt this exercise ourselves, however.  The copylefted files
      in question are:

         compress.c -- this is GPL'ed and provides CFITSIO's ability
            to read gzip'ed FITS files.

         wcsutil.c -- this is LGPL'ed and provides CFITSIO's ability
            to do some of its own World Coordinate Systems
            transformations.  I don't believe, however, that we are
            actually making use of any of this functionality, as we
            use libwcs instead to perform WCS transformations.
            Consequently, it may be safe to turn off this feature of
            CFITSIO.

   libwcs -- written by Doug Mink of the Smithsonian Astrophysical
      Observatory, this World Coordinate Systems library is LGPL'ed.
      As we currently do not dynamically link this library, however,
      fits2itk executables currently pick up all the requirements of
      the LGPL.

fits2itk also bundles the following libraries that are not distributed
under the MIT License, but which are distributed under licenses that
are for all intents and purposes equivalent to the MIT License:

   pathToExecutable -- Adapted for general use by Douglas Alan from
      Python's "getpath.c".  It carries the permissive Python Software
      Foundation License Version 2.


BUILDING FIT2ITK FROM SOURCE
----------------------------

At the moment, fits2itk can only be built on Unix-like OS'es.  A
Windows release is planned "sometime soon".

fits2itk uses CMake as its "make" system.  In order to build fits2itk
from source, you will first need to acquire CMake if it is not already
installed on your computer.  CMake is available under a BSD-like
permissive free license via most package repositories, such as Yum,
APT, pkgsrc, Fink, MacPorts, etc.  If you prefer, you can also
download the source release of CMake, or a prebuilt binary directly
from

   http://www.cmake.org

Once you have CMake on your computer, you will need to download ITK
(The National Library of Medicine Insight Segmentation and
Registration Toolkit).  It is available with a BSD-like permissive
free license from:

   http://itk.org

You will only need to download the single tarball

   InsightToolkit-3.2.0.tar.gz

from the ITK distribution.

(By the time you read this, there may or may not be a newer version of
ITK, which may or may not work--although it most probably will--with
fits2itk.)

You can install ITK to wherever in your filesystem you would like it
to live, but in the following instructions, I will assume you are
installing the sources into /usr/local/src/itk and that you will build
the toolkit into /usr/local/lib/itk.  If you end up choosing different
locations for the software to live, then you will have to modify the
following instructions accordingly.

Build ITK as root (or as another user that has r/w access to
/usr/local):

    % mkdir /usr/local/src/itk
    % cd /usr/local/src/itk
    % tar xzf /path/to/InsightToolkit-3.2.0.tar.gz
    % mkdir /usr/local/lib/itk
    % cd /usr/local/lib/itk
    % cmake -DCMAKE_BUILD_TYPE=Release /usr/local/src/itk/InsightToolkit-3.2.0
    % make

The above "make" will probably take a while to complete.  Once it has
finished building, you are ready to build fits2itk.  Again, you can
build and install fits2itk wherever you would like, but again I am
going to assume that you will be placing the sources into
/usr/local/src/fits2itk, that you will build fits2itk in
/usr/local/lib/fits2itk, and that you will install the executable into
/usr/local/bin:

    % cd /usr/local/src
    % tar xjf /path/to/fits2itk-src-0.3.1.tar.bz2
    % mkdir /usr/local/lib/fits2itk
    % cd /usr/local/lib/fits2itk
    % setenv ITK_DIR /usr/local/lib/itk
    % cmake -DCMAKE_BUILD_TYPE=Release /usr/local/src/fits2itk
    % make
    % strip fits2itk
    % mv fits2itk /usr/local/bin

For help on how to use fits2itk, type

    % fits2itk -h

and see additional documentation at

    http://astromed.iic.harvard.edu/FITS-reader

For a tutorial on how to use fits2itk with 3D Slicer, see

    http://astromed.iic.harvard.edu/UsingSlicer


SUBMITTING YOUR IMPROVEMENTS
----------------------------

We use the Mercurial distributed version control system.  If perchance
you become so motivated as to code up a significant improvement to
fits2itk on your own, you can use Mercurial to write out a patch file
containing your modifications.  You can then email us the patch file
for consideration into the mainline development tree.

Mercurial can be found in most package repositories, or can be
downloaded from

   http://www.selenic.com/mercurial

Only submissions that bear a permissive free license, unencumbered by
any restrictions will be considered for inclusion into the mainline
sources.  E.g., the 3-clause BSD License and the MIT License are OK,
but GPL or LGPL are not.  The "MIT License" that is used in our code
is preferred, as it places fewer restrictions on the use of the
software than any other common freeware license.

To use Mercurial to generate the appropriate patch file, you would
first "clone" the fits2itk source directory that you downloaded:

   % hg clone /usr/local/src/fits2itk /some/path/fits2itk-myclone
   % cd /some/path/fits2itk-myclone

Once you have cloned a source repository of fits2itk, Mercurial can
track any changes that you make to the cloned version as compared to
the version you originally downloaded.  So, to submit your
modifications, you would edit the sources files in the directory
fits2itk-myclone.  You would then test out your changes by compiling
and running your modified version of fits2itk.  To do this, you will
need to make a build directory for your changes like so:

   % mkdir /some/path/fits2itk-mybuild
   % cd /some/path/fits2itk-mybuild
   % cmake -DCMAKE_BUILD_TYPE=Debug /some/path/fits2itk-myclone
   % make

Once you are convinced that your enhancements are bug-free, you can
get a list of the files that you changed by doing

   % cd /some/path/fits2itk-myclone
   % hg status

If you made any new files (e.g., "newfile.cxx"), you can tell
Mercurial about the new file like so:

   % hg add newfile.cxx

Before you can create a patch file, you first need to commit your
changes into your clone.  (Doing so does not commit the changes into
/usr/local/src/fits2itk -- it only commits them into your clone
repository.  To commit your changes into /usr/local/src/fits2itk, you
would use "hg push".  But you don't want to do that, since you need to
keep /usr/local/src/fits2itk as an unsullied copy of *our* repository
in order for Mercurial to be able generate the proper diffs of your
mods.)  To commit your changes into the clone, do

   % hg commit

Once you have committed your changes, you can create the patch file
like so:

   % hg export tip > fits2itk-mods.patch

You would them email us fits2itk-mods.patch.
