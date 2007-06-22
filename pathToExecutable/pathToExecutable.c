// -*- Mode: C; fill-column: 79 -*-

//=============================================================================
// Description:
//
//      pathToExecutable.c is a module that allows a Unix program to find the
//      location of its executable.  This capability is extremely useful for
//      writing programs that don't have to recompiled in order to be relocated
//      within the filesystem.  Any auxiliary files (dynamically loaded
//      libraries, help files, configuration files, etc.) can just be placed in
//      the same directory as the executable, and the function
//      pathToExecutable() can be used by the program at runtime to locate its
//      executable file and from there the program can locate any auxiliary
//      files it needs in order to operate.
//
//      pathToExecutable() is smart enough to follow a symlink (or even a chain
//      of symlinks) in order to find the true location of the executable.  In
//      this manner, for instance, you might install all of the files used by a
//      program (let's say it's called "my-program"), including the executable,
//      into the directory /usr/local/lib/my-program, and then put a symlink
//      into /usr/local/bin that points to the executable
//      /usr/local/lib/my-program/my-program.  Initially pathToExecutable()
//      will identify /usr/local/bin/my-program as the executable, but it will
//      then notice that this "file" is really a symbolic link.
//      pathToExecutable() will then follow the symbolic link and return
//      "/usr/local/lib/my-program/my-pogram" instead.
//
//      Before a program can call pathToExecutable(), setArgv() must be called
//      (canonically in main()) so that pathToExecutable() can fetch the value
//      of argv[0] and use it to help figure out where the executable is
//      located.
//
// Copyright and licensing information:
//
//      This software is a heavily modified version of getpath.c from the
//      Python 2.5.1 release.  Both this software and the original software
//      from which it is derived are freely distributable under the terms of
//      the permissive freeware license, Python Software Foundation License
//      Version 2.  You can read more about this license here:
//
//           http://www.python.org/psf/license
//
//      The original software from which this software is derived carries the
//      following copyright:
//
//         Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007 Python
//         Software Foundation.
//
//      The modifications to the original software, which are contained herein,
//      are
//
//         Copyright (c) 2007 Douglas Alan <doug AT alum.mit.edu>
//
//=============================================================================

#include <libgen.h>       // For basename()
#include <stdio.h>
#include <stdlib.h>       // For getenv()
#include <string.h>
#include <unistd.h>

#include <sys/types.h>    // For stat()
#include <sys/stat.h>     //  "   "

#include <sys/param.h>    // For MAXPATHLEN

#include <pathToExecutable.h>

#define c_filenameSeperator  '/'
#define c_pathDelimiter ':'

//-----------------------------------------------------------------------------
// Local functions
//-----------------------------------------------------------------------------

static void
fatalError(const char* const message)
{
  if (getArgc()) {
    const char* progname = basename(getArgv()[0]);
    if (progname != NULL) fprintf(stderr, "%s: ", progname);
  }
  fprintf(stderr, "FATAL ERROR: %s\n", message);
  exit(1);
}


//=============================================================================
// Exported functions
//=============================================================================

static int fs_argc = 0;
static const char* const* fs_argv;

//-----------------------------------------------------------------------------
// Exported global functions:
//
//    setArgv()
//    pathToExecutable()
//
//    These are the two functions that you really care about.
//-----------------------------------------------------------------------------

void
setArgv(const int argc, const char* const argv[])
{
  static bool firstTime = true;
  if (!firstTime) fatalError("setArgv() called more than once.");
  firstTime = false;
  fs_argc = argc;
  fs_argv = argv;
}


// Before any searches are done, the location of the executable is determined.
// If argv[0] has one or more slashes in it, it is used unchanged.  Otherwise,
// it must have been invoked from the shell's path, so we search $PATH for the
// named executable and use that.  If the executable was not found on $PATH (or
// there was no $PATH environment variable), the original argv[0] string is
// used.

// Next, the executable location is examined to see if it is a symbolic link.
// If so, the link is chased (correctly interpreting a relative pathname if one
// is found) and the directory of the link target is used.

// Finally, the pathname to the executable is returned as pathToExecutable()'s
// return value.

const char*
pathToExecutable(void)
{
  if (!getArgc()) {
    fatalError("pathToExecutable() can't get argv[0].");
  }
  static bool firstTime = true;
  static char retval[MAXPATHLEN+1];
  if (firstTime) {
    firstTime = false;
    static const char delimiter[2] = {c_pathDelimiter, '\0'};
    static const char separator[2] = {c_filenameSeperator, '\0'};
    const char* const argv0 = getArgv()[0];
    char* path = getenv("PATH");
    size_t bufsz;


#ifdef __APPLE__
#if MAC_OS_X_VERSION_MAX_ALLOWED >= MAC_OS_X_VERSION_10_4
    uint32_t nsexeclength = MAXPATHLEN;
#else
    unsigned long nsexeclength = MAXPATHLEN;
#endif // MAC_OS_X_VERSION_10_4
#endif // __APPLE__

    // If there is a slash in argv0, then we know that location of the
    // executable was specified explicitly, rather than being run from
    // $PATH:
    if (strchr(argv0, c_filenameSeperator)) {
      strncpy(retval, argv0, MAXPATHLEN);
      retval[MAXPATHLEN] = '\0';
    }

#ifdef __APPLE__

    else if (0 == _NSGetExecutablePath(retval, &nsexeclength) &&
	     retval[0] == c_filenameSeperator)
      {
	// NOP: `retval` is set by _NSGetExecutablePath above, and that is
	// all we need to have done here.
      }

#endif // __APPLE__

    // If there is no slash in the argv0 path, then we have to assume the
    // executable is on the user's $PATH.  (If $PATH isn't exported, you
    // lose):
    else if (path) {

      // Loop over all the directories in $PATH looking for an executable file
      // whose name is the same as argv0.  If we find one, then we have, in
      // theory, located the executable:
      while (true) {
	char* delimPtr = strchr(path, c_pathDelimiter);
	if (delimPtr) {
	  size_t len = delimPtr - path;
	  if (len > MAXPATHLEN) len = MAXPATHLEN;
	  strncpy(retval, path, len);
	  *(retval + len) = '\0';
	} else {
	  strncpy(retval, path, MAXPATHLEN);
	  retval[MAXPATHLEN] = '\0';
	}
	pteJoinPathM(retval, argv0);
	if (pteIsAnExecutableFile(retval)) break;
	if (!delimPtr) {
	  retval[0] = '\0';
	  break;
	}
	path = delimPtr + 1;
      }
    }

    // If we've gotten here, then it appears that we are out of luck, and we'll
    // just have to set retval to the empty string:
    else retval[0] = '\0';
  
    // Hmmm, I see no reason to absolutize the path as is done in the
    // following commented-out line of code.  Doing so can be
    // detrimental when using some network filesystem automounters,
    // for which 'pwd' doesn't really work right.  (Actually, I think
    // actually do always work right, as long as you can be sure that
    // the filesystem won't be auto-unmounted while you might still
    // want to use the path.  I think that as long as the executable
    // is still running from the filesystem, then you can be pretty
    // darn sure that the filesystem won't be unmounted.)

    // if (retval[0] != c_filenameSeperator) absolutize(retval);

    // We've now (hopefully) located the executable, sort of, but we're not
    // completely done, as we may have really only located a symbolic link to
    // the executable.  If that's the case, we have to chase through any chain
    // of symlinks that we might find until we get to the real file:
    {
      char tmpbuffer[MAXPATHLEN+1];
      int linklen = readlink(retval, tmpbuffer, MAXPATHLEN);
      while (linklen != -1) {

	// The retval from readlink is not null terminated, so we need to
	// make it so:
	tmpbuffer[linklen] = '\0';
	if (tmpbuffer[0] == c_filenameSeperator) {
	  strncpy(retval, tmpbuffer, MAXPATHLEN);
	  retval[MAXPATHLEN] = '\0';
	} else {
	  // If we've gotten here, then the symlink is a relative link, so
	  // interpret it relative to retval:
	  pteDirname(retval);
	  pteJoinPathM(retval, tmpbuffer);
	}
	linklen = readlink(retval, tmpbuffer, MAXPATHLEN);
      }
    }

    // Hmmmm, I think we'll return the entire path to the executable, just in
    // case the caller might find that information useful, rather than just the
    // path to the containing directory.  Consequently, I've commented out the
    // following line:
    // pteDirname(retval);

  } // end if (firstTime)

  return retval;
}


//-----------------------------------------------------------------------------
// Other useful exported global functions
//-----------------------------------------------------------------------------

const char*
pathToExecutableDir(void)
{
  static bool firstTime = true;
  static char retval[MAXPATHLEN+1];
  if (firstTime) {
    firstTime = false;
    const char* pathToExe = pathToExecutable();
    strcpy(retval, pathToExecutable());
    pteDirname(retval);
  }
  return retval;
}

int
getArgc(void)
{
  return fs_argc;
}


const char* const*
getArgv(void)
{
  return fs_argv;
}


// Given a string such as "/path/to/a/file", modifies it to be "/path/to/a".

void
pteDirname(char *filepath)
{
    size_t i = strlen(filepath);
    while (i > 0 && filepath[i] != c_filenameSeperator)
        --i;
    filepath[i] = '\0';
}


// Add a filepath component, by appending `pathTail` to `pathHead`.
// `pathHead` must have at least MAXPATHLEN + 1 bytes allocated, and
// contain a null-terminated string with no more than MAXPATHLEN
// characters (not counting the trailing null).  It's a fatal error if
// it contains a string longer than that (callers must be careful!).
// If these requirements are met, it's guaranteed that `pathHead` will
// still be a NUL-terminated string with no more than MAXPATHLEN
// characters at exit.  If `pathTail` is too long, only as much of
// stuff as fits will be appended.

void
pteJoinPathM(char* const pathHead, const char* const pathTail)
{
    size_t n, k;
    if (pathTail[0] == c_filenameSeperator) {
      n = 0;
    } else {
        n = strlen(pathHead);
        if (n > 0 && pathHead[n-1] != c_filenameSeperator && n < MAXPATHLEN) {
	  pathHead[n++] = c_filenameSeperator;
	}
    }
    k = strlen(pathTail);
    if (n + k > MAXPATHLEN) {
      fatalError("Buffer overflow in pteJoinPathM().");
    }
    strncpy(pathHead + n, pathTail, k);
    pathHead[n + k] = '\0';
}


const char*
pteJoinPath(const char* const pathHead, const char* const pathTail)
{
  static bool firstTime = true;
  static char* buffer;
  if (firstTime) {
    firstTime = false;
    buffer = malloc(MAXPATHLEN+1);
  }
  strncpy(buffer, pathHead, MAXPATHLEN);
  buffer[MAXPATHLEN] = '\0';
  pteJoinPathM(buffer, pathTail);
  return buffer;
}


bool
pteIsFile(const char* const filepath)
{
    struct stat buf;
    if (stat(filepath, &buf) != 0) return false;
    else if (!S_ISREG(buf.st_mode)) return false;
    else return true;
}


bool
pteIsAnExecutableFile(const char* const filepath)
{
    struct stat buf;
    if (stat(filepath, &buf) != 0) return false;
    else if (!S_ISREG(buf.st_mode)) return false;
    else if ((buf.st_mode & 0111) == 0) return true;
    else return true;
}


bool
pteIsDir(const char* const filepath)
{
    struct stat buf;
    if (stat(filepath, &buf) != 0) return false;
    else if (!S_ISDIR(buf.st_mode)) return false;
    else return true;
}


// pteCopyAbsolute() requires that `pathHead` be allocated at least MAXPATHLEN
// + 1 bytes and that pathTail be no more than MAXPATHLEN bytes.

void
pteCopyAbsolute(char* const destStr, const char* filepath)
{
    if (filepath[0] == c_filenameSeperator)
        strcpy(destStr, filepath);
    else {
        getcwd(destStr, MAXPATHLEN);
        if (filepath[0] == '.' && filepath[1] == c_filenameSeperator)
            filepath += 2;
        pteJoinPathM(destStr, filepath);
    }
}


// pteAbsolutize() requires that `filepath` be allocated at least MAXPATHLEN+1
// bytes.

void
pteAbsolutize(char* const filepath)
{
    char buffer[MAXPATHLEN + 1];
    if (filepath[0] == c_filenameSeperator) return;
    pteCopyAbsolute(buffer, filepath);
    strcpy(filepath, buffer);
}
