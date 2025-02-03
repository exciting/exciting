# libbzint
#
# This cmake files compiles the libBZINT library
#

# Specify the folder containing the source files for the libBZINT library
set(FOLDER_LIBBZINT   "${CMAKE_SOURCE_DIR}/src/src_libbzint")

# Recursively find all Fortran source files in the specified folder.
# This includes files with extensions .f90, .F90, and .f.
file(GLOB_RECURSE LIBBZINT_SRC_FILES
    "${FOLDER_LIBBZINT}/*.f90"
    "${FOLDER_LIBBZINT}/*.F90"
    "${FOLDER_LIBBZINT}/*.f"
)

# Create a static library named 'bzint' using the Fortran source files found above.
add_library(bzint STATIC ${LIBBZINT_SRC_FILES})
