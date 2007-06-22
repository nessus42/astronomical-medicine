# -*- Mode: CMake; fill-column: 79 -*-

# TODO: Break out stuff that is not configuration into an included macro
# definitions file, so that this file can contain almost only stuff that might
# change from day to day.

project (fits2itk)

# Specify the source files that constitute the fits2itk program:
set (fits2itk_source_files
  fits2itk.cxx version.cxx itkFITSImageIO.cxx itkFITSImageIOFactory.cxx
  itkFITSWCSTransform.cxx
  )

# Specify libraries used by fits2itk that are packaged with the fits2itk
# source files as subdirs:
set (subdirs
  da_util
  # fits2itk_cfitsio
  libwcs
  pathToExecutable
  )

# Put the fits2itk main source directory into the C/C++ include path.  I do
# this because my C++ code uses '#include <foo.h>', rather than '#include
# "foo.h"' for its own header files.  In the past I've had bad luck with
# inconsistent interpretation of the latter form of "#include" by compilers:
include_directories(${PROJECT_SOURCE_DIR})

# For each library packaged with the fits2itk sources...
foreach (subdir ${subdirs})

  # ... put its subdir into the include path, so our "#include" statements can
  # find the library's header files:
  include_directories (${PROJECT_SOURCE_DIR}/${subdir})

  # ... and build the library:
  add_subdirectory (${PROJECT_SOURCE_DIR}/${subdir})

  endforeach (subdir ${subdirs})

# Find the ITK build.  The location of the ITK build can be specified by
# setting the ITK_DIR environment variable or by specifying
# "-DITK_DIR=/path/to/itk-build" on the cmake command line.  We must find ITK
# like so before doing "add_executable" for a program that uses ITK:
find_package (ITK REQUIRED)
if (ITK_FOUND)
  include (${ITK_USE_FILE})
  endif (ITK_FOUND)


# The following link_directories command has to occur before the following
# add_executable command.  Strangely enough, the corresponding
# include_directories comand can occur anywhere in the file:
link_directories (${PROJECT_BINARY_DIR}/cfitsio)
include_directories (${PROJECT_BINARY_DIR}/cfitsio)

# Specify the executables to build:
add_executable (fits2itk
  ${fits2itk_source_files}
)

# Specify the libraries that are used by the fits2itk program...
target_link_libraries (fits2itk

  # ... including the ones that come bundled with fits2itk in its subdirs:
  ${subdirs}
  ITKIO 
  ITKCommon
  cfitsio
  )


# BEGIN building the CFITSIO library, using it's own autoconf-based build
# system.

set_source_files_properties (itkFITSImageIO.cxx
  PROPERTIES
  OBJECT_DEPENDS ${PROJECT_BINARY_DIR}/cfitsio/libcfitsio.a
  )

add_custom_command (

  OUTPUT ${PROJECT_BINARY_DIR}/cfitsio/libcfitsio.a
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${PROJECT_SOURCE_DIR}/cfitsio
                         ${PROJECT_BINARY_DIR}/cfitsio
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "Running \\\"./configure\\\" on CFITSIO library..."
  COMMAND ${CMAKE_COMMAND}
  ARGS -E chdir ${PROJECT_BINARY_DIR}/cfitsio /bin/sh -c "./configure"
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo
  "Ignore the above message about running \\\"make\\\"\; I will do that for"
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "you right now."
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "Building CFITSIO library with \\\"make\\\"..."
  COMMAND ${CMAKE_COMMAND}
  ARGS -E chdir ${PROJECT_BINARY_DIR}/cfitsio /bin/sh -c
             "env MAKEFLAGS='' make"
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "Done building CFITSIO."
  )

# END building CFITSIO.


# Here we specify the files that we need to copy into the build directory so
# that they can be found at runtime by fits2itk during code development.  We do
# this here by declaring the compilation of version.cxx to be dependent on
# these files.  As this dependency is not really true, this approach is a bit
# of a kludge.  Nonetheless, this mechanisms seems to be the most convenient
# way to specify this sort of thing in CMake.

# Please note that you can only call `set_source_files_properties` once for
# each soure file.  (If you call it again later, the next call will completely
# overwrite the previous settings.)  So, if you want to pile on more
# dependencies onto this version.cxx, you have to do it here.

# Also, you can only apply a single `add_custom_command` to a source file in
# this manner.  This restriction is counterintuitive, as `add_custom_command`
# seems to apply to the files that need to be copied or created, not to the
# source file that has been declared to depend on them.  But so it goes.
# Consequently, if any one of the files specified below need to be copied or
# created, all of them will be.  (This is one of the reasons that we have
# separated the following stuff out from the building of CFITSIO, above.  It
# didn't make sense to have to rebuild all of CFITSIO just because
# `shortUsageMessage.txt` needed to be copied in.

set_source_files_properties (version.cxx
  PROPERTIES
  OBJECT_DEPENDS ${PROJECT_BINARY_DIR}/shortUsageMessage.txt
  OBJECT_DEPENDS ${PROJECT_BINARY_DIR}/verboseUsageMessage.txt
  )

add_custom_command (

  # The following will make sure that files that are needed by fits2itk at
  # runtime will be copied into the build directory so that fits2itk can find
  # them:
  OUTPUT ${PROJECT_BINARY_DIR}/shortUsageMessage.txt
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_SOURCE_DIR}/shortUsageMessage.txt
               ${PROJECT_BINARY_DIR}
  OUTPUT ${PROJECT_BINARY_DIR}/verboseUsageMessage.txt
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_SOURCE_DIR}/verboseUsageMessage.txt
               ${PROJECT_BINARY_DIR}
  )
