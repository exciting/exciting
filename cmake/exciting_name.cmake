# exiting_name.cmake
# 
# This CMake file sets the appropiate name for exciting executable
# using set_exciting_name function 
#

# This CMake function sets the program name based on the provided flags for
# OpenMP (OMP) and MPI (MPI). It defines the name of the program as 
# follows:
# - If MPI is enabled:
#   - If OpenMP is enabled: "exciting_mpismp"
#   - If OpenMP is disabled: "exciting_purempi"
# - If MPI is disabled:
#   - If OpenMP is enabled: "exciting_smp"
#   - If OpenMP is disabled: "exciting_serial"

function(set_exciting_name IS_MPI IS_OMP)

    # Check the conditions based on OMP and MPI flags
    if(IS_MPI)
	if(IS_OMP)
	    set(exciting_name "exciting_mpismp" PARENT_SCOPE)
        else()
	    set(exciting_name "exciting_purempi" PARENT_SCOPE)
        endif()
    else()
	if(IS_OMP)
	    set(exciting_name "exciting_smp" PARENT_SCOPE)
        else()
	    set(exciting_name "exciting_serial" PARENT_SCOPE)
        endif()
    endif()

endfunction()

# Run the name selector and print the name
set_exciting_name(${MPI} ${OMP})
message(STATUS "Set exciting executable to: ${exciting_name}")


