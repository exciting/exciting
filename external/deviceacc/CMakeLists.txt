cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

#Some options
option(NVIDIA "NVIDIA GPU" OFF)
option(AMD "AMD GPU" OFF)
option(INTEL "INTEL GPU" OFF)
option(CPUBACKEND "If enabled activates the CPU backend so GPU accelerated routines are replaced by their CPU backend" OFF)

project(deviceacc Fortran CXX)

if(NOT CPUBACKEND)
  # Now select the vendor
  if(NVIDIA)
    message(STATUS "GPU :  NVIDIA")
    add_compile_definitions(NVIDIAGPU)
  elseif(AMD)
    message(STATUS "GPU :  AMD")
    add_compile_definitions(AMDGPU)
    set(AMDTARGET "None" CACHE PATH "Target march for AMD")
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
    find_library(rocfftlib NAMES rocfft librocfft REQUIRED)
    message(STATUS "rocFFT found : ${rocfftlib}")
  elseif(INTEL)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "Should not be here")
  endif()

  # Select compilers: Intel or GNU compilers
  if(CMAKE_Fortran_COMPILER MATCHES "ifx")
      message(STATUS "Intel compilers")
      if (NOT INTEL)
        message(FATAL_ERROR "For non Intel cards please use GNU compilers. Exiting.")
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
  elseif(CMAKE_Fortran_COMPILER MATCHES "gfortran")
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
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto -O3 -fPIC -fopenmp \
                              -foffload-options=-lm ")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -flto -O3 -fPIC \
                                  -ffree-line-length-512 -fopenmp -fopenmp-allocators  \
                                  -foffload-options=-lgfortran -foffload-options=-lm")
      endif()
      set(devacc_link_libs "openblasp")
  else()
      message(FATAL_ERROR "Compiler is not recognized: only GNU and Intel (ifx) compilers are supported")
  endif()
else()
  if(CMAKE_Fortran_COMPILER MATCHES "ifx")
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
  elseif(CMAKE_Fortran_COMPILER MATCHES "gfortran")
      message(STATUS "GNU compilers")
      set(CMAKE_CXX_COMPILER "g++")
      set(CMAKE_Fortran_COMPILER "mpif90")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -flto  -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -flto -cpp -O3 -fPIC \
                                -ffree-line-length-512 -fopenmp")
      find_package(FFTW3 REQUIRED)
      include_directories(${FFTW3_INCLUDE_DIRS})                 
      set(devacc_link_libs "-lopenblasp -lfftw3 -lfftw3f")
  else()
      message(FATAL_ERROR "Compiler is not recognized: only GNU and Intel (ifx) compilers are supported")
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


