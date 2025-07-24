##############################
# Create objects using CMake #
##############################

# Command to execute XSLT transformation for input.xsd

find_program(XSLTPROC_EXECUTABLE xsltproc)
if (XSLTPROC_EXECUTABLE)
    message(STATUS "xsltproc found: ${XSLTPROC_EXECUTABLE}")
else ()
    message(FATAL_ERROR "xsltproc not found. Please install it.")
endif ()

execute_process(COMMAND ${XSLTPROC_EXECUTABLE} ${CMAKE_SOURCE_DIR}/xml/schema/schemaexpand.xsl ${CMAKE_SOURCE_DIR}/xml/schema/input.xsd
                OUTPUT_FILE ${CMAKE_SOURCE_DIR}/xml/excitinginput.xsd)

# Command to execute XSLT transformation for inputmodules.f90
execute_process(COMMAND ${XSLTPROC_EXECUTABLE} ${CMAKE_SOURCE_DIR}/xml/schematofortran.xsl ${CMAKE_SOURCE_DIR}/xml/excitinginput.xsd
                OUTPUT_FILE ${CMAKE_SOURCE_DIR}/src/src_inputparser/inputmodules.f90)

# Command to execute XSLT transformation for speciesmodules.f90
execute_process(COMMAND ${XSLTPROC_EXECUTABLE} ${CMAKE_SOURCE_DIR}/xml/schematofortran.xsl ${CMAKE_SOURCE_DIR}/xml/species.xsd
                OUTPUT_FILE ${CMAKE_SOURCE_DIR}/src/src_inputparser/speciesmodules.f90)

# Get src/version.inc
# Function to get compiler version
function(GetCompilerVersion outVar)
    # Execute the compiler command and capture the output
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version OUTPUT_VARIABLE compiler_output)
     string(REGEX REPLACE "\n.*" "" compiler_output_first_line "${compiler_output}")
    set(${outVar} "${compiler_output_first_line}" PARENT_SCOPE)
endfunction()

# Get Git information
set(GITHASH  "sodium")
set(GITHASH2 "alpha")

# Get compiler version
GetCompilerVersion(COMPILERVERSION)

# Get date
string(TIMESTAMP CURRENT_DATE "%y,%m,%d")

# Write version information to file
file(WRITE  ${CMAKE_SOURCE_DIR}/src/version.inc "#define GITHASH \"${GITHASH}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.inc "#define GITHASH2 \"${GITHASH2}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.inc "#define COMPILERVERSION \"${COMPILERVERSION}\"\n")
file(APPEND ${CMAKE_SOURCE_DIR}/src/version.inc "#define VERSIONFROMDATE /${CURRENT_DATE}/\n")

###################################################################################
# In Cray 18 there is an incompatibility with what is generated and CMake expects #
###################################################################################

if (CMAKE_VERSION VERSION_LESS_EQUAL "3.31.3" AND CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    string(REGEX MATCH "^([0-9]+)" COMPILER_MAJOR_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
    message(STATUS "Correcting ${CMAKE_Fortran_COMPILER_VERSION}")
    if (COMPILER_MAJOR_VERSION AND COMPILER_MAJOR_VERSION VERSION_GREATER_EQUAL "18")
	# Create the modules folder
	file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/)
	# Create the symlink
	execute_process(
		COMMAND ln -s ${CMAKE_BINARY_DIR}/modules/m_memory_device.SM_MEMORY_DEVICE.smod
	                      ${CMAKE_BINARY_DIR}/modules/sm_memory_device.mod
        )
    endif()
endif()

