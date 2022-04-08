FetchContent_Declare(
  cbc
  URL https://www.coin-or.org/download/source/Cbc/Cbc-2.10.5.tgz
  URL_HASH MD5=46277180c0fc67f750e2de1836333189
)
include(ProcessorCount)
ProcessorCount(N)
if (NOT N)
  set(N 1)
endif()
FetchContent_GetProperties(cbc)
if(NOT cbc_POPULATED)
  FetchContent_Populate(cbc)
  message(STATUS "Building CBC library from source.")
  include(ExternalProject)
  ExternalProject_Add(cbc
      SOURCE_DIR ${cbc_SOURCE_DIR}
      CONFIGURE_COMMAND ${cbc_SOURCE_DIR}/configure --enable-cbc-parallel --enable-openmp --disable-zlib --disable-bzlib --without-blas --without-lapack --enable-static --disable-shared --with-pic --prefix=${cbc_BINARY_DIR}
      BUILD_BYPRODUCTS ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Cgl${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Clp${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Cbc${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Osi${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}OsiCbc${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}OsiClp${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}CoinUtils${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ClpSolver${CMAKE_STATIC_LIBRARY_SUFFIX} 
        ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}CbcSolver${CMAKE_STATIC_LIBRARY_SUFFIX} 
      BUILD_COMMAND make -j${N}
      BUILD_IN_SOURCE 1
      )
  set(cbc_LIBRARIES
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}OsiCbc${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}CbcSolver${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Cbc${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Cgl${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}OsiClp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ClpSolver${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Clp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Osi${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${cbc_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}CoinUtils${CMAKE_STATIC_LIBRARY_SUFFIX}
    )
endif()
