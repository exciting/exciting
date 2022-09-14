# Find Python
find_package(Python3 3.7 COMPONENTS Interpreter Development)
if(Python3_FOUND)
    message("-- Python 3 interpreter version: " ${Python3_VERSION})
else()
    message("-- Python 3 interpreter not found")
endif()

# Python helper functions

function(find_python_module module)
    # Find if a Python module is installed. Copied from
    # https://github.com/ivansafrin/Polycode/blob/master/CMake/FindPythonModule.cmake

    string(TOUPPER ${module} module_upper)
    if(NOT PY_${module_upper})
        if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
            set(${module}_FIND_REQUIRED TRUE)
        endif()
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND "${Python3_EXECUTABLE}" "-c"
                "import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
                RESULT_VARIABLE _${module}_status
                OUTPUT_VARIABLE _${module}_location
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT _${module}_status)
            set(PY_${module_upper} ${_${module}_location} CACHE STRING
                    "Location of Python module ${module}")
        endif(NOT _${module}_status)
    endif(NOT PY_${module_upper})
    find_package_handle_standard_args(PY_${module} DEFAULT_MSG PY_${module_upper})
endfunction(find_python_module)

function(find_pythonhome)
    if ($ENV{PYTHONHOME})
        set(_PYTHONHOME $ENV{PYTHONHOME})
    else()
        execute_process(
                COMMAND "${Python3_EXECUTABLE}" "-c"
                "import sys; print(sys.prefix + ':' + sys.exec_prefix)"
                RESULT_VARIABLE _PYTHONHOME_FAILED
                OUTPUT_VARIABLE _PYTHONHOME
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(_PYTHONHOME_FAILED)
            message(FATAL_ERROR "Could not determine PYTHONHOME. Error:" ${_PYTHONHOME_FAILED})
        endif()
    endif()
    message("    Found PYTHONHOME: " ${_PYTHONHOME})
    set(PYTHONHOME ${_PYTHONHOME} PARENT_SCOPE)
endfunction(find_pythonhome)
