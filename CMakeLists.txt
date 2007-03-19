PROJECT(FITSImageIO)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

SUBDIRS(
  da_util
  cfitsio
  libwcs
  )

INCLUDE_DIRECTORIES(
  .
  ${FITSImageIO_SOURCE_DIR}/da_util
  ${FITSImageIO_SOURCE_DIR}/cfitsio
  ${FITSImageIO_SOURCE_DIR}/libwcs
  )

IF(NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${FITSImageIO_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
IF(NOT EXECUTABLE_OUTPUT_PATH)
  SET (EXECUTABLE_OUTPUT_PATH ${FITSImageIO_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)
 
LINK_DIRECTORIES(
  ${FITSImageIO_BINARY_DIR}/bin
  )

ADD_EXECUTABLE(fits2itk
  fits2itk.cxx version.cxx itkFITSImageIO.cxx itkFITSImageIOFactory.cxx
  itkFITSWCSTransform.cxx
  )

TARGET_LINK_LIBRARIES(fits2itk da_util itkcfitsio itkwcs ITKIO ITKCommon)
