# This CMake installs IDieL

option(IDIEL "Install IDieL" ON)

if (IDIEL)
    add_compile_definitions(IDIEL)
    message(STATUS "IDieL : ON")

    set(IDieLInstallDir "${CMAKE_BINARY_DIR}/INTERNAL_IDieL_install")
    file(MAKE_DIRECTORY ${IDieLInstallDir})
    include(ExternalProject)

    ExternalProject_Add(INTERNAL_IDieL
      SOURCE_DIR    "${CMAKE_SOURCE_DIR}/external/IDieL"
      BINARY_DIR    "${CMAKE_BINARY_DIR}/INTERNAL_IDieL_build"
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -S <SOURCE_DIR> -B <BINARY_DIR>
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DMPI=${MPI}
        -DOMP=${OMP}
        -DHDF5=${HDF5}
        -DMKL=${MKL}
	-DTESTS=OFF
        -DOPENBLAS=${OPENBLAS}
        -DAMDLINALG=${AMDLINALG}
        -DCRAYLIBSCI=${CRAYLIBSCI}
        -DOTHERLINALG=${OTHERLINALG}
        -DLINALGLIB=${LINALGLIB}
        -DFFTW3_ROOT=${FFTW3_ROOT}
        -DNVIDIA=${NVIDIA}
        -DNVIDIAARCH=${NVIDIAARCH}
        -DAMD=${AMD}
        -DAMDTARGET=${AMDTARGET}
        -DAMD_HIPSETVALIDDEVICE_SUPPORTED=${AMD_HIPSETVALIDDEVICE_SUPPORTED}
        -DINTEL=${INTEL}
        -DUSM=${USM}
        -DCPUBACKEND=${CPUBACKEND}
        -DMAGMA_DIR=${MAGMA_DIR}
        -DCMAKE_INSTALL_PREFIX=${IDieLInstallDir}
      BUILD_COMMAND $(MAKE)
      INSTALL_COMMAND $(MAKE) install
    )

    add_library(libIDieL SHARED IMPORTED GLOBAL)
    set_target_properties(libIDieL PROPERTIES
       IMPORTED_LOCATION "${IDieLInstallDir}/lib/libIDieL.so"
    )
    include_directories(${IDieLInstallDir}/include/)
    add_dependencies(libIDieL INTERNAL_SPGLIB)

    # We need to install it
    install(FILES "${IDieLInstallDir}/lib/libIDieL.so"
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )


else()
    message(STATUS "IDieL : OFF")
    add_library(libIDieL INTERFACE)
endif()
