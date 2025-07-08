# Include this file in your build to use the tools.
# This is necessary because of the linker’s limitations 
# and offload table in GNU compilers.
# For ifx and Cray this problem does not exist, so technically a simple link should do the magic.
# TODO: Check new flang compiler once is usable and with OpenMP support

# Protect the file against double inclusion
include_guard(GLOBAL)

# Compile options
option(NVIDIA "NVIDIA GPU" OFF)
option(AMD "AMD GPU" OFF)
option(AMD_HIPSETVALIDDEVICE_SUPPORTED "Turn ON if besides defined, hipSetValidDevices is supported by the driver" OFF)
option(INTEL "INTEL GPU" OFF)
option(USM "Turn on if your GPU and CPU share physical memory")
if(NVIDIA OR AMD OR INTEL)
    option(CPUBACKEND "If enabled activates the CPU backend so GPU accelerated routines are replaced by their CPU backend" OFF)
else()
    option(CPUBACKEND "If enabled activates the CPU backend so GPU accelerated routines are replaced by their CPU backend" ON)
endif()

if(NOT DEVICEACC_SOURCE_DIR)
    message(FATAL_ERROR "Set DEVICEACC_SOURCE_DIR in your CMake to the path where this file is located")
endif()

# Classic intel compiler has small bug, here we tackle it
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
   execute_process(
    COMMAND ifort --version
    OUTPUT_VARIABLE IFORT_VERSION_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   # Parse the version number from the output (optional)
   string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" IFORT_VERSION "${IFORT_VERSION_OUTPUT}")

   # Check if the compiler version is 2021.4.0 or lower
   if(${IFORT_VERSION} VERSION_LESS_EQUAL 2021.4.0)
      message(STATUS "Correction for ifort 2021.4.0 bug on pure function with assumed ranks: ${IFORT_VERSION}")
      add_compile_definitions(IFORT_2021_PURE_ASSUMED_RANK_BUG)
   endif()
endif()

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

  if (USM)
    message(STATUS "GPU-CPU (Unified Shared Memory) : ON")
    add_compile_definitions(_USM_)
  else()
    message(STATUS "GPU-CPU (Unified Shared Memory) : OFF")
  endif()

  if(NVIDIA OR AMD)
    # Setting MAGMA library for the NVIDIA and AMD case
    set(MAGMA_F90_DIR ${DEVICEACC_SOURCE_DIR}/external/magma)
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

      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fPIC -fiopenmp -fopenmp-targets=spir64 -qmkl=parallel -fsycl")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fopenmp-targets=spir64=\"-fp-model=precise\" -fsycl -fpp -free")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ax -fPIC -fiopenmp -fopenmp-targets=spir64 -qmkl=parallel -fsycl")
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fopenmp-targets=spir64=\"-fp-model=precise\" -fsycl -fpp -free")

      set(devacc_link_libs "${MKL_LIBRARIES}")

  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      message(STATUS "GNU compilers")
      if(NOT NVIDIA AND NOT AMD)
        message(FATAL_ERROR "For Intel cards please use Intel compilers. Exiting.")
      endif()

      if(AMD)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fPIC -fopenmp-allocators -fopenmp -foffload-options=-march=${AMDTARGET}")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fopenmp-allocators -foffload-options=-march=${AMDTARGET}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fPIC -fopenmp-allocators -fopenmp -foffload-options=-march=${AMDTARGET}")
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fopenmp-allocators -foffload-options=-march=${AMDTARGET}")
      else()
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fPIC -fopenmp -foffload-options=-lm -foffload=nvptx-none -fopenmp-allocators -foffload-options=-march-map=sm_${NVIDIAARCH}")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}  -foffload=nvptx-none -foffload-options=-lgfortran -fopenmp-allocators -foffload-options=-lm -foffload-options=-march-map=sm_${NVIDIAARCH}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fPIC -fopenmp -foffload-options=-lm -fopenmp-allocators -foffload=nvptx-none -foffload-options=-march-map=sm_${NVIDIAARCH}")
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -foffload=nvptx-none -foffload-options=-lgfortran -fopenmp-allocators -foffload-options=-lm -foffload-options=-march-map=sm_${NVIDIAARCH}")
      endif()
      set(devacc_link_libs "${BLAS_LIBRARIES}")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
      if(NOT NVIDIA AND NOT AMD)
        message(FATAL_ERROR "For Intel cards please use Intel compilers. Exiting.")
      endif()
      # For Cray add this instruction whenever USM is required
      if (USM)
        set(USM_FLAGS "-fopenmp-force-usm")
      endif()
      # Cray should find the proper targets with the acceleration modules
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fopenmp -h acc_model=auto_async_none:no_fast_addr:deep_copy ${USM_FLAGS}")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fPIC -fopenmp")
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fopenmp -h acc_model=auto_async_none:no_fast_addr:deep_copy ${USM_FLAGS}")
  else()
      message(FATAL_ERROR "Compiler is not recognized: only GNU, Intel (ifx) and Cray compilers are supported")
  endif()
else()
   if(NOT FFTW3_LIBRARIES AND NOT MKL)
       message(FATAL_ERROR "Deviceacc requires FFTW3. Library not found.")
   endif()
   if(NOT BLAS_LIBRARIES AND NOT MKL)
       message(FATAL_ERROR "Deviceacc requires BLAS and LAPACK routines. Library not found.")
   endif()
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fPIC")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fPIC")
   set(devacc_link_libs "${BLAS_LIBRARIES};${FFTW3_LIBRARIES}")
endif()

# Set source files
set(SRC_DEVICEACC "${SRC_MAGMA_F90}")
set(DEVICEACC_SRC_DIR "${DEVICEACC_SOURCE_DIR}/src")

# Depend on the backend
if(NOT CPUBACKEND)
    # Control
    set(SRC_DEVICEACC_CONTROL ${DEVICEACC_SRC_DIR}/control/device/device_world_t.f90)
    #Linear algebra
    set(SRC_DEVICEACC_LINALG ${DEVICEACC_SRC_DIR}/linalg/device/linalg_device_common.f90)
    # FFT
    if (AMD)
        set(SRC_DEVICEACC_FFT ${DEVICEACC_SOURCE_DIR}/external/hipfort/hipfort_rocfft_enums.f90
                              ${DEVICEACC_SOURCE_DIR}/external/hipfort/hipfort_rocfft.f90)
    endif()
    set(SRC_DEVICEACC_FFT ${SRC_DEVICEACC_FFT} ${DEVICEACC_SRC_DIR}/fft/device/fft_device_t.f90)
    # Allocation
    set(SRC_DEVICEACC_MEMORY ${DEVICEACC_SRC_DIR}/memory/common/memory_device.f90
	                     ${DEVICEACC_SRC_DIR}/memory/common/memory_device.cpp)
    set(SRC_DEVICEACC_MEMORY ${SRC_DEVICEACC_MEMORY} ${DEVICEACC_SRC_DIR}/memory/device/s_memory_device.f90)
    # Macros
    set(DEVICEACC_MACROS ${DEVICEACC_SRC_DIR}/macros/device/offload.fpp)
    include_directories(${DEVICEACC_SRC_DIR}/macros/device/)
else()
    set(SRC_DEVICEACC_CONTROL ${DEVICEACC_SRC_DIR}/control/host/device_world_t.f90)
    set(SRC_DEVICEACC_LINALG  ${DEVICEACC_SRC_DIR}/linalg/host/linalg_device_common.f90)
    set(SRC_DEVICEACC_FFT     ${SRC_DEVICEACC_FFT} ${DEVICEACC_SRC_DIR}/fft/host/fft_device_t.f90)
    set(SRC_DEVICEACC_MEMORY  ${DEVICEACC_SRC_DIR}/memory/common/memory_device.f90)
    set(SRC_DEVICEACC_MEMORY  ${SRC_DEVICEACC_MEMORY} ${DEVICEACC_SRC_DIR}/memory/common/memory_device.cpp)
    set(SRC_DEVICEACC_MEMORY  ${SRC_DEVICEACC_MEMORY} ${DEVICEACC_SRC_DIR}/memory/host/s_memory_device.f90)
    set(DEVICEACC_MACROS      ${DEVICEACC_SRC_DIR}/macros/host/offload.fpp)
    include_directories(${DEVICEACC_SRC_DIR}/macros/host/)
endif()

# Create deviceacc object
set(SRC_DEVICEACC ${SRC_MAGMA_F90}
                  ${SRC_DEVICEACC_CONTROL}
                  ${SRC_DEVICEACC_LINALG}
                  ${SRC_DEVICEACC_FFT}
                  ${SRC_DEVICEACC_MEMORY})

if(NOT CPUBACKEND)
    # Create a link variable for the link against it and not its dependencies
    if (NVIDIA)
        set(LINKS_DEVICEACC ${devacc_link_libs} ${MAGMA_DIR}/lib/libmagma.so CUDA::cufft CUDA::cudart)
    elseif(AMD)
        set(LINKS_DEVICEACC ${devacc_link_libs} ${MAGMA_DIR}/lib/libmagma.so ${rocfftlib})
    elseif(INTEL)
        set(LINKS_DEVICEACC ${devacc_link_libs})
    else()
        message(FATAL_ERROR "Should not be here")
    endif()

else()
    set(LINKS_DEVICEACC ${devacc_link_libs})
endif()
