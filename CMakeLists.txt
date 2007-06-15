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
  fits2itk_da_util
  # fits2itk_cfitsio
  fits2itk_libwcs
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


# This will build the bundled CFITSIO library, using it's own autoconf-based
# build system:

set_source_files_properties (${fits2itk_source_files}
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
