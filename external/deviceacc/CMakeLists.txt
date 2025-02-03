cmake_minimum_required(VERSION 3.27.0 FATAL_ERROR)

# Compile options
option(NVIDIA "NVIDIA GPU" OFF)
option(AMD "AMD GPU" OFF)
option(AMD_HIPSETVALIDDEVICE_SUPPORTED "Turn ON if besides defined, hipSetValidDevices is supported by the driver" OFF)
option(INTEL "INTEL GPU" OFF)
if(NVIDIA OR AMD OR INTEL)
    option(CPUBACKEND "If enabled activates the CPU backend so GPU accelerated routines are replaced by their CPU backend" OFF)
else()
    option(CPUBACKEND "If enabled activates the CPU backend so GPU accelerated routines are replaced by their CPU backend" ON)
endif()

project(deviceacc Fortran CXX)

add_compile_definitions(DEVICEACC)

# In case we want device supported compilation
# find and properly set libraries
if(NOT CPUBACKEND)
  # Set a general compiler option to protect
  add_compile_definitions(DEVICEOFFLOAD)
  # Now select the vendor
  if(NVIDIA)
    message(STATUS "GPU :  NVIDIA")
    add_compile_definitions(NVIDIAGPU)
    set(NVIDIAARCH "89" CACHE STRING "NVIDIA architecture (e.g. 80, 89)")
    message(STATUS "NVIDIA architecture : ${NVIDIAARCH}")
  elseif(AMD)
    message(STATUS "GPU :  AMD")
    add_compile_definitions(AMDGPU)
    set(AMDTARGET "None" CACHE STRING "Target march for AMD")
    if(AMDTARGET MATCHES "None")
        message(FATAL_ERROR "For AMD GPU the march of the target device is needed.")
    endif()
    message(STATUS "AMD target : ${AMDTARGET}")
  elseif(INTEL)
    message(STATUS "GPU :  INTEL")
    add_compile_definitions(INTELGPU)
  else()
    message(FATAL_ERROR "GPU vendor has not been selected, supported GPUs are: AMD, NVIDIA, and INTEL")
  endif()

  if(NVIDIA OR AMD)
    # Setting MAGMA library for the NVIDIA and AMD case
    set(MAGMA_F90_DIR ${PROJECT_SOURCE_DIR}/external/magma)
    set(SRC_MAGMA_F90
          ${MAGMA_F90_DIR}/magma2.F90
          ${MAGMA_F90_DIR}/magma2_common.F90
          ${MAGMA_F90_DIR}/magma2_sfortran.F90
          ${MAGMA_F90_DIR}/magma2_cfortran.F90
          ${MAGMA_F90_DIR}/magma2_dfortran.F90
          ${MAGMA_F90_DIR}/magma2_zfortran.F90)
    set(MAGMA_DIR "None" CACHE PATH "Path to magma library")
    message(STATUS "MAGMA library path : ${MAGMA_DIR}")

    # Checking version
    set(MAJOR_SEARCH_STRING "MAGMA_VERSION_MAJOR")
    set(MINOR_SEARCH_STRING "MAGMA_VERSION_MINOR")
    set(MICRO_SEARCH_STRING "MAGMA_VERSION_MICRO")
    file(STRINGS ${MAGMA_DIR}/include/magma_types.h FILE_CONTENTS)
    set(MAGMA_MAJOR 0)
    set(MAGMA_MINOR 0)
    set(MAGMA_MICRO 0)

    foreach(LINE IN LISTS FILE_CONTENTS)
      if("${LINE}" MATCHES "${MAJOR_SEARCH_STRING}")
        string(REPLACE " " ";" LINE_ELEMENTS "${LINE}")
        list(GET LINE_ELEMENTS 2 THIRD_ELEMENT)
        set(MAGMA_MAJOR ${THIRD_ELEMENT})
      endif()
      if("${LINE}" MATCHES "${MINOR_SEARCH_STRING}")
        string(REPLACE " " ";" LINE_ELEMENTS "${LINE}")
        list(GET LINE_ELEMENTS 2 THIRD_ELEMENT)
        set(MAGMA_MINOR ${THIRD_ELEMENT})
      endif()
      if("${LINE}" MATCHES "${MICRO_SEARCH_STRING}")
        string(REPLACE " " ";" LINE_ELEMENTS "${LINE}")
        list(GET LINE_ELEMENTS 2 THIRD_ELEMENT)
        set(MAGMA_MICRO ${THIRD_ELEMENT})
      endif()
    endforeach()

    if ((MAGMA_MAJOR EQUAL 0) AND (MAGMA_MINOR EQUAL 0) AND (MAGMA_MICRO EQUAL 0))
      message(FATAL_ERROR "MAGMA version not found")
    endif()

    set(MAGMA_VERSION_NUM ${MAGMA_MAJOR}.${MAGMA_MINOR}${MAGMA_MICRO})

    if (MAGMA_VERSION_NUM LESS 2.72)
      message(FATAL_ERROR "MAGMA version needs to be at least 2.7.2")
    endif()

    message(STATUS "MAGMA version : ${MAGMA_MAJOR}.${MAGMA_MINOR}.${MAGMA_MICRO}")

  else()
    # In INTEL case one uses the MKL offloading to the GPU through openMP
    set(MAGMA_F90_DIR "")
    set(SRC_MAGMA_F90 "")
  endif()

  # Find MPI
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_libraries(${MPI_Fortran_LIBRARIES})

  if(NVIDIA)

    # Find CUDA
    find_package(CUDAToolkit REQUIRED)

    if(CUDAToolkit_FOUND)
        message(STATUS "CUDAToolkit found. Include directories: ${CUDAToolkit_INCLUDE_DIRS}")
        message(STATUS "CUDAToolkit libraries: ${CUDAToolkit_LIBRARIES}")
    else()
        message(FATAL_ERROR "CUDAToolkit not found.")
    endif()

  elseif(AMD)
    # Find rocFFT
    find_package(rocfft REQUIRED)
    set(rocfftlib ${ROCFFT_LIBRARIES})
    message(STATUS "rocFFT found : ${rocfftlib}")
    
    message(STATUS "ROCM include directory: ${ROCFFT_INCLUDE_DIRS}")
    file(READ ${ROCFFT_INCLUDE_DIRS}/rocm-core/rocm_version.h FILE_CONTENT)

    # Use regular expressions to extract version numbers
    string(REGEX MATCH "#define ROCM_VERSION_MAJOR[ \t]+([0-9]+)" _major_match "${FILE_CONTENT}")
    string(REGEX MATCH "#define ROCM_VERSION_MINOR[ \t]+([0-9]+)" _minor_match "${FILE_CONTENT}")
    string(REGEX MATCH "#define ROCM_VERSION_PATCH[ \t]+([0-9]+)" _patch_match "${FILE_CONTENT}")

    # If not found, default to 0
    if(_major_match)
        string(REGEX REPLACE "#define ROCM_VERSION_MAJOR[ \t]+([0-9]+)" "\\1" ROCM_VERSION_MAJOR "${_major_match}")
    else()
        set(ROCM_VERSION_MAJOR 0)
    endif()

    if(_minor_match)
        string(REGEX REPLACE "#define ROCM_VERSION_MINOR[ \t]+([0-9]+)" "\\1" ROCM_VERSION_MINOR "${_minor_match}")
    else()
        set(ROCM_VERSION_MINOR 0)
    endif()

    if(_patch_match)
        string(REGEX REPLACE "#define ROCM_VERSION_PATCH[ \t]+([0-9]+)" "\\1" ROCM_VERSION_PATCH "${_patch_match}")
    else()
        set(ROCM_VERSION_PATCH 0)
    endif()

    if ((ROCM_VERSION_MAJOR EQUAL 0) AND (ROCM_VERSION_MINOR EQUAL 0) AND (ROCM_VERSION_PATCH EQUAL 0))
        message(FATAL_ERROR "ROCm version not found")
    else()
	message(STATUS "ROCM version : ${ROCM_VERSION_MAJOR}.${ROCM_VERSION_MINOR}.${ROCM_VERSION_PATCH}")
    endif()

    set(ROCM_VERSION_NUM ${ROCM_VERSION_MAJOR}.${ROCM_VERSION_MINOR}${ROCM_VERSION_PATCH})

    if (ROCM_VERSION_NUM GREATER_EQUAL 6.22)
	if (AMD_HIPSETVALIDDEVICE_SUPPORTED)
	    message(STATUS "ROCM defines and supports setting valid devices")
            add_compile_definitions(AMD_CAN_SET_VALID_DEVICES)
	else()
	    message(STATUS "ROCM defines but don't supports setting valid devices")
	endif()
    else()
	message(STATUS "ROCM does not support setting valid devices")
    endif()

  elseif(INTEL)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "Only NVIDIA, AMD, or INTEL GPUs are supported. This part of the code 
            is unreachable by construction. Consequently, this error indicates some
	    catastrophic CMake configuration step.")
  endif()

  # Few compiler options
  if (CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
      message(STATUS "Intel compilers")
      if (NOT INTEL)
          message(FATAL_ERROR "For non Intel cards please use GNU/Cray compilers. Exiting.")
      endif()

      set(CMAKE_CXX_COMPILER "icpx")
      set(CMAKE_Fortran_COMPILER "mpiifx")

      SET(CMAKE_AR ${CMAKE_CXX_COMPILER_AR} CACHE PATH "AR" FORCE)
      SET(CMAKE_RANLIB ${CMAKE_CXX_COMPILER_RANLIB} CACHE PATH "RANLIB" FORCE)

      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -fiopenmp \
                           -fopenmp-targets=spir64 -qmkl=parallel -fsycl")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fpp -O3 -fPIC \
                              -no-wrap-margin -fp-model source -fiopenmp \
                              -fopenmp-targets=spir64 -qmkl=parallel -fsycl -fpp -free")
      set(devacc_link_libs "${MKL_LIBRARIES}")

  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      message(STATUS "GNU compilers")
      set(CMAKE_CXX_COMPILER "g++")
      set(CMAKE_Fortran_COMPILER "mpif90")
      if(NOT NVIDIA AND NOT AMD)
        message(FATAL_ERROR "For Intel cards please use Intel compilers. Exiting.")
      endif()
      if(AMD)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -fopenmp -foffload-options=-march=${AMDTARGET}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -fPIC \
                              -ffree-line-length-512 -fopenmp -fopenmp-allocators  \
                              -foffload-options=-march=${AMDTARGET}")
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -fopenmp \
                              -foffload-options=-lm ")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -fPIC \
                                  -ffree-line-length-512 -fopenmp -fopenmp-allocators  \
                                  -foffload-options=-lgfortran -foffload-options=-lm")
      endif()

      # Find FFTW3 library
      find_library(FFTW3_LIBRARIES NAMES fftw3)
      find_library(FFTW3F_LIBRARIES NAMES fftw3f)
      find_library(FFTW3_LIBRARIES_OMP NAMES fftw3_omp)
      find_library(FFTW3F_LIBRARIES_OMP NAMES fftw3f_omp)
      find_path(FFTW3_INCLUDE_DIR NAMES "fftw3.h")
      include_directories(${FFTW3_INCLUDE_DIR})

      if (NOT FFTW3_LIBRARIES)
        message(FATAL_ERROR "FFTW3 not found. FTTW3 is required for GNU/Cray build of deviceacc.")
      endif()

      set(devacc_link_libs "${LINALG_LIBRARY};${FFTW3_LIBRARIES};${FFTW3F_LIBRARIES};${FFTW3_LIBRARIES_OMP};${FFTW3F_LIBRARIES_OMP}")

  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
      if(NOT NVIDIA AND NOT AMD)
        message(FATAL_ERROR "For Intel cards please use Intel compilers. Exiting.")
      endif()
      # Cray should find the proper targets with the acceleration modules
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fopenmp")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fopenmp")


      # Find FFTW3 library
      find_library(FFTW3_LIBRARIES NAMES fftw3)
      find_library(FFTW3F_LIBRARIES NAMES fftw3f)
      find_library(FFTW3_LIBRARIES_OMP NAMES fftw3_omp)
      find_library(FFTW3F_LIBRARIES_OMP NAMES fftw3f_omp)
      find_path(FFTW3_INCLUDE_DIR NAMES "fftw3.h")
      include_directories(${FFTW3_INCLUDE_DIR})

      if (NOT FFTW3_LIBRARIES)
        message(FATAL_ERROR "FFTW3 not found. FTTW3 is required for GNU/Cray build of deviceacc.")
      endif()

      set(devacc_link_libs "${LINALG_LIBRARY};${FFTW3_LIBRARIES};${FFTW3F_LIBRARIES};${FFTW3_LIBRARIES_OMP};${FFTW3F_LIBRARIES_OMP}")

  else()
      message(FATAL_ERROR "Compiler is not recognized: only GNU and Intel (ifx) compilers are supported")
  endif()
else()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
      message(STATUS "Intel compilers")
      set(CMAKE_CXX_COMPILER "icpx")
      set(CMAKE_Fortran_COMPILER "mpiifx")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -fiopenmp")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fpp -O3 -fPIC \
                              -no-wrap-margin -fp-model source -fiopenmp \
                              -qmkl=parallel -free")
      find_package(MKL REQUIRED)
      include_directories(${MKL_INCLUDE_DIRS})
      set(devacc_link_libs "${MKL_LIBRARIES}")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      message(STATUS "GNU compilers")
      set(CMAKE_CXX_COMPILER "g++")
      set(CMAKE_Fortran_COMPILER "mpif90")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -fPIC \
                                -ffree-line-length-512 -fopenmp")
      # Find linear algebra
      find_library(LINALG_LIBRARY REQUIRED NAMES openblasp openblas_mp openblas_omp openblas)
      if (NOT ${LINALG_LIBRARY})
        message(FATAL_ERROR "OpenBLAS not found. OpenBLAS is required for GNU/Cray build of deviceacc.")
      endif()

      # Find FFTW3 library
      find_library(FFTW3_LIBRARIES NAMES fftw3)
      find_library(FFTW3F_LIBRARIES NAMES fftw3f)
      find_library(FFTW3_LIBRARIES_OMP NAMES fftw3_omp)
      find_library(FFTW3F_LIBRARIES_OMP NAMES fftw3f_omp)
      find_path(FFTW3_INCLUDE_DIR NAMES "fftw3.h")
      include_directories(${FFTW3_INCLUDE_DIR})

      if (NOT FFTW3_LIBRARIES)
        message(FATAL_ERROR "FFTW3 not found. FTTW3 is required for GNU/Cray build of deviceacc.")
      endif()

      set(devacc_link_libs "${LINALG_LIBRARY};${FFTW3_LIBRARIES};${FFTW3F_LIBRARIES};${FFTW3_LIBRARIES_OMP};${FFTW3F_LIBRARIES_OMP}")

  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
      set(CMAKE_CXX_COMPILER "CC")
      set(CMAKE_Fortran_COMPILER "mpiftn")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fopenmp")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fopenmp")

      # Find FFTW3 library
      find_library(FFTW3_LIBRARIES NAMES fftw3)
      find_library(FFTW3F_LIBRARIES NAMES fftw3f)
      find_library(FFTW3_LIBRARIES_OMP NAMES fftw3_omp)
      find_library(FFTW3F_LIBRARIES_OMP NAMES fftw3f_omp)
      find_path(FFTW3_INCLUDE_DIR NAMES "fftw3.h")
      include_directories(${FFTW3_INCLUDE_DIR})

      if (NOT FFTW3_LIBRARIES)
        message(FATAL_ERROR "FFTW3 not found. FTTW3 is required for GNU/Cray build of deviceacc.")
      endif()

      set(devacc_link_libs "${LINALG_LIBRARY};${FFTW3_LIBRARIES};${FFTW3F_LIBRARIES};${FFTW3_LIBRARIES_OMP};${FFTW3F_LIBRARIES_OMP}")

  else()
      message(FATAL_ERROR "Compiler is not recognized: only GNU, Intel, Cray, and LLVMFlang compilers are supported")
  endif()

endif()

# Create INSTALL DIRECTORY
set(INSTALL_DIR "${PROJECT_SOURCE_DIR}/deviceacc")
message(STATUS "INSTALL DIRECTORY : ${PROJECT_SOURCE_DIR}")
file(MAKE_DIRECTORY ${INSTALL_DIR})
file(MAKE_DIRECTORY ${INSTALL_DIR}/lib)

# Set the components source
set(SRC "${SRC_MAGMA_F90}")
set(MODDIR "")

# Process source
add_subdirectory(src)

# Compile test
add_subdirectory(test)


