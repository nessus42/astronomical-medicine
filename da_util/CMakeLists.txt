# -*- Mode: CMake; fill-column: 79 -*-

add_library (da_util 
  da_usual.cxx da_util.cxx RcPointer.cxx
  )

include_directories (${PROJECT_SOURCE_DIR}/pathToExecutable)
