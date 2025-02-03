# linalg.cmake
# This CMake script configures the linear algebra libraries (BLAS and LAPACK) for the project.
# It supports multiple libraries including Intel MKL, OpenBLAS, and user-defined libraries.
# Note that for BLAS and LAPACK, this script can find vendor-specific versions.
# However, in some builds, it may not locate parallel versions, defaulting to the serial version in OpenBLAS.
# Refer to the CMake documentation for more information: 
# https://cmake.org/cmake/help/v3.15/module/FindBLAS.html

# Options to select which linear algebra library to use
option(OPENBLAS "Use OpenBLAS" OFF)
option(AMDLINALG "Use AMD libraries for the linear algebra" OFF)
option(OTHERLINALG "Use other linear algebra library" OFF)

# Check for conflicts
set(LINALG_COUNT 0)
if(OPENBLAS)
  math(EXPR LINALG_COUNT "${LINALG_COUNT} + 1")
endif()
if(AMDLINALG)
  math(EXPR LINALG_COUNT "${LINALG_COUNT} + 1")
endif()
if(OTHERLINALG)
  math(EXPR LINALG_COUNT "${LINALG_COUNT} + 1")
endif()
if(MKL)
  math(EXPR LINALG_COUNT "${LINALG_COUNT} + 1")
endif()
if(LINALG_COUNT GREATER 1)
	message(FATAL_ERROR "You cannot set OPENBLAS, Intel MKL, AMDLINALG or other linear 
	algebra library at the same time (selected: ${LINALG_COUNT}).")
endif()

if(MKL)
    # Check if MKL is already found before attempting to locate it
    if(NOT MKL_FOUND)
	set(MKL_INTERFACE "lp64")
	if (NOT ${OMP})
            set(MKL_THREADING "sequential")
	endif()
	find_package(MKL CONFIG REQUIRED)
        message(STATUS "MKL found (DIR)  : ${MKL_DIR}")
    endif()
    
    # Set up MKL as the BLAS library
    set(BLAS_LIBRARIES MKL::MKL)

elseif(OPENBLAS)
    # Configure for OpenBLAS
    if(OMP)
        find_library(OpenBLAS_LIBRARY NAMES openblasp openblas_mp openblas_omp openblas)
    else()
        find_library(OpenBLAS_LIBRARY NAMES openblas)
    endif()
    set(BLAS_LIBRARIES "${OpenBLAS_LIBRARY}")
    message(STATUS "OpenBLAS found : ${OpenBLAS_LIBRARY}")

elseif(AMDLINALG)
   # Configure for AMD libraries
   if(OMP)
   	find_library(BLIS_LIBRARY NAMES blis-mt blis)
   else()
	find_library(BLIS_LIBRARY NAMES blis)
   endif()
   find_library(FLAME_LIBRARY NAMES flame)

   if (NOT FLAME_LIBRARY OR NOT BLIS_LIBRARY)
   	message(FATAL_ERROR "AMD linear algebra libraries (BLIS and libFLAME) were not found.")
   endif()

   set(BLAS_LIBRARIES "${BLIS_LIBRARY}")
   list(APPEND BLAS_LIBRARIES ${FLAME_LIBRARY})
   message(STATUS "AMD linear algebra libraries found : ${BLAS_LIBRARIES}")

elseif(OTHERLINALG)
    # Support for user-defined linear algebra libraries (non-CMake vendors like libsci from Cray)
    set(LINALGLIB "None" CACHE STRING "Linear algebra library")
    
    if(LINALGLIB STREQUAL "None")
        message(FATAL_ERROR "For OTHERLINALG=ON please provide the library with its path using LINALGLIB variable")
    endif()
    
    set(BLAS_LIBRARIES "${LINALGLIB}")
    message(WARNING "Be aware that this library may be untested; please run the test cases (e.g. libsci fails in some cases)")
    message(STATUS "Linear algebra library found : ${LINALGLIB}")

else()
    # Default to using BLA_VENDOR to find BLAS and LAPACK
    message(STATUS "Using BLA_VENDOR option to search for the library.")
    message(WARNING "Falling back to default BLAS and LAPACK from the system. That might be unoptimal.")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    # Combine LAPACK libraries into BLAS libraries
    list(APPEND BLAS_LIBRARIES ${LAPACK_LIBRARIES})
    message(STATUS "Linear algebra library found : ${BLAS_LIBRARIES}")
endif()

