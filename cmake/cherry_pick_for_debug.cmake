# cherry_pick_for_debug.cmake
#
# This CMake file implements a function to apply debug flags to cherry-picked Fortran files
# Arguments:
#   SRC - list of source files
#   CHERRY_PICK_FILES_FOR_DEBUG - semicolon-separated list of filenames (with paths)
#   CHERRY_PICK_FILES_DEBUG_FLAGS - flags to apply to the cherry picked files.
function(apply_cherry_pick_debug SRC CHERRY_PICK_FILES_FOR_DEBUG CHERRY_PICK_FILES_DEBUG_FLAGS)

    # Do nothing if no cherry-pick files are provided
    if("${CHERRY_PICK_FILES_FOR_DEBUG}" STREQUAL "")
        return()
    endif()

    # Do nothing if no debug flags are provided
    if("${CHERRY_PICK_FILES_DEBUG_FLAGS}" STREQUAL "")
        message(STATUS "No debug flags provided for cherry-picked files")
        return()
    endif()

    # Set the cherry picking debug flags
    set(cherry_flags "${CHERRY_PICK_FILES_DEBUG_FLAGS}")

    message(STATUS "Applying debug flags to cherry-picked files")
    message(STATUS "Debug flags for cherry-picked files: ${cherry_flags}")

    # Loop over all source files
    foreach(src_file ${SRC})
        foreach(cherry_file ${CHERRY_PICK_FILES_FOR_DEBUG})
            if("${src_file}" STREQUAL "${CMAKE_SOURCE_DIR}/src/${cherry_file}")
                set_property(
                    SOURCE "${src_file}"
                    APPEND_STRING
                    PROPERTY COMPILE_FLAGS " ${cherry_flags}"
                )
                message(STATUS "Applying cherry-pick debug flags to ${src_file}")
            endif()
        endforeach()
    endforeach()
endfunction()
