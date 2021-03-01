cmake_minimum_required(VERSION 3.13)

project(equation_of_state VERSION 1.4.0
                          DESCRIPTION "Equation of state (EOS) program for fitting energy-volume data"
		                  LANGUAGES Fortran)

set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

add_executable(eos)

set_target_properties(eos
        PROPERTIES
        RUNTIME_OUTPUT_NAME eos.exe
        VERSION 1.4.0
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY}
        )

target_sources(eos
        PRIVATE
        src/eos.f90       
        src/eveos.f90     
        src/fitdata.f90   
        src/fmin_nm.f90   
        src/getedata.f90  
        src/minf_nm.f90   
        src/modmain.f90   
        src/output.f90    
        src/pveos.f90     
        src/readinput.f90
        )
