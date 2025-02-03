# manage_exciting_source.cmake
#
# This cmake manages the exciting source files
# and set appropiate properties for them.
# exciting files are then

# Define the current folder
set(EXCITINGLIB_SRC_FOLDERS ${CMAKE_SOURCE_DIR}/src)
set(EXCITINGLIB_SRC_FILES "")
set(EXCITINGLIB_SRC_FILES_FIXED "")

# Specify folders to ignore
set(FOLDER_SPACEGROUP "${CMAKE_SOURCE_DIR}/src/spacegroup")
set(FOLDER_LIBBZINT   "${CMAKE_SOURCE_DIR}/src/src_libbzint")

# Specify files not inside ignored folders to ignore
# TODO: Remove or check those files
set(EXCLUDE_SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/species/species.f90"
    # Add more files to exclude as needed
)

# We need special treatment of files containing INQUIRE so that it compiles using the Cray compiler
# Developers, please add here anycase of INQUIRE using not the default integer
set(EXCITING_INQUIRE_FILES "")
if (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")

    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_SOURCE_DIR}/src/src_xs/putgeteps0.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_SOURCE_DIR}/src/src_xs/src_rttddft/rttddft_io_serial.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_SOURCE_DIR}/src/src_gw/calcpmatgw.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_SOURCE_DIR}/src/src_gw/src_tests/checkmbrot.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_SOURCE_DIR}/src/src_hybrids/putvxnl.f90")

    # We create an INQUIRE_FLAGS
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(INQUIRE_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG}")
    else()
	set(INQUIRE_FLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
    endif()
    string(REPLACE "integer32" "integer64" INQUIRE_FLAGS "${INQUIRE_FLAGS}")
    set_source_files_properties(${EXCITING_INQUIRE_FILES} PROPERTIES COMPILE_FLAGS ${INQUIRE_FLAGS})
endif()

# Iterate over files and purge main
foreach(folder ${EXCITINGLIB_SRC_FOLDERS})
    # Search for F90 files in the current folder
    file(GLOB_RECURSE F90_FILES_IN_FOLDER "${folder}/*.f90" "${folder}/*.F90")

    # Exclude folder as this
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_SPACEGROUP}/.*$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_LIBBZINT}/.*$")

    # Iterate over files and exclude those in EXCLUDE_FILES list
    foreach(file ${EXCLUDE_SRC_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()

    # Iterate over files and exclude those in INQUIRE_FILES list
    foreach(file ${EXCITING_INQUIRE_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()

    # Exclude main.f90 and main.F90
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.f90$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.F90$")

    # Append the filtered files to the list
    list(APPEND EXCITINGLIB_SRC_FILES ${F90_FILES_IN_FOLDER})
endforeach()

foreach(folder ${EXCITINGLIB_SRC_FOLDERS})
    # Search for F90 files in the current folder
    file(GLOB_RECURSE F90_FILES_IN_FOLDER "${folder}/*.f")

    # Exclude folder as this
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_SPACEGROUP}/.*$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_LIBBZINT}/.*$")

    # Iterate over files and exclude those in EXCLUDE_FILES list
    foreach(file ${EXCLUDE_SRC_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()

    # Iterate over files and exclude those in INQUIRE_FILES list
    foreach(file ${EXCITING_INQUIRE_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()


    # Exclude main.f90 and main.F90
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.f90$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.F90$")

    # Append the filtered files to the list
    list(APPEND EXCITINGLIB_SRC_FILES_FIXED ${F90_FILES_IN_FOLDER})
endforeach()

# New Intel compiler fails for some reason to automatically do it (2024.02)
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        foreach(source_file ${EXCITINGLIB_SRC_FILES_FIXED})
                set_source_files_properties(${source_file} PROPERTIES
                        COMPILE_OPTIONS "-fixed"
                )
        endforeach()
endif()

# Fuse all sources
set(EXCITING_SOURCE_FILES ${EXCITINGLIB_SRC_FILES} ${EXCITINGLIB_SRC_FILES_FIXED} ${EXCITING_INQUIRE_FILES})

# Set the exciting main source
set(EXCITING_MAIN "${CMAKE_SOURCE_DIR}/src/mainxml/main.f90")
