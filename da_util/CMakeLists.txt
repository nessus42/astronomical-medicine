# -*- Mode: CMake; fill-column: 79 -*-

# For showFactoryClasses():
find_package (ITK REQUIRED)  #d 
include (${ITK_USE_FILE}) #d
# include_directories (/Users/nessus/s/am/src/fits2itk.dev) #d

add_library (da_util 
  da_usual.cxx da_util.cxx RcPointer.cxx
  )

include_directories (${PROJECT_SOURCE_DIR}/pathToExecutable)
