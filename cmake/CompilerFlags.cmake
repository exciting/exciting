# Set compiler flags
# CMake will append CMAKE_Fortran_FLAGS with CMAKE_Fortran_FLAGS_BUILDTYPE
# CMAKE_Fortran_FLAGS_BUILDTYPE may also have predefined values, hence initialise it

# Standard flags to use in all cases
set(STD_FFLAGS
    -std='f2008'            # Fortran standard set to 2008
    -fimplicit-none         # Specify that no implicit typing is allowed
    -ffree-line-length-0    # No fixed line length
    -march=native           # Produces code optimized for the local machine under the constraints of the
                            # selected instruction set (hence the result might not run on different machines).
   )
string(REPLACE ";" " " STD_FFLAGS "${STD_FFLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${STD_FFLAGS}")

# GCC
set(GCC_DEBUG
     -g               # Generate symbols
     -fbacktrace      # symbolic stack traceback
     -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
     -finit-real=nan  #  All real scalars are initialised to NaN
     -fcheck=all      # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
    )
# More debug flags to consider:
# -finit-integer=2147483647 -finit-real=snan \
# -frecord-gcc-switches -finit-character=42 -finit-logical=true -fdump-core -fstack-protector-all -pipe


set(GCC_RELEASE -O3)  # Level 3 optimisation. Could also consider -fpack-derived

# Intel
set(INTEL_DEBUG
    -g            # Generate symbols
    -traceback    # symbolic stack traceback
    -fp           # Disables the ebp register in optimizations and sets the ebp register to be used as the frame pointer.
    -check all    # Checks for all runtime failures.
    -check bounds # Generates code to perform runtime checks on array subscript and character substring expressions.
    -check-uninit #  Enables runtime checking for uninitialized variables.
    -ftrapuv      #  Set unassigned scalars as a very large integer or an invalid address
    -fpe3         # control over floating-point exception (divide by zero, overflow, invalid operation, underflow, denormalized number, positive infinity, negative infinity or a NaN)
    )

# -g -O0 -DUSE_ASSERT -debug all -implicitnone -warn unused \
#   -fp-stack-check -heap-arrays -ftrapuv -check pointers \
#   -check bounds -check all -check noarg_temp_created -traceback -I"${MKLROOT}/include"

set(INTEL_RELEASE
    -O3           # Optimsation level 3
    -no-prec-div  # Heurisitics to improves precision of floating-point divides: enables optimizations that give slightly less precise results than full IEEE division.
    -fp-model fast=2 # Semantics of floating-point calculations: Enables more aggressive optimizations on floating-point data.
    -foptimize-sibling-calls # Optimise tail recursive calls.
   )

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set(FF_DEBUG ${GCC_DEBUG})
   set(FF_RELEASE ${GCC_RELEASE})

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set(FF_DEBUG ${INTEL_DEBUG})
   set(FF_RELEASE ${INTEL_RELEASE})

else ()
     message(SEND_ERROR "flags have not been defined for this compiler: \
            ${CMAKE_Fortran_COMPILER_ID}")
endif()

string(REPLACE ";" " " FF_DEBUG "${FF_DEBUG}")
string(REPLACE ";" " " FF_RELEASE "${FF_RELEASE}")

# Initialise BUILDTYPE flags so we completely define/control
# the compiler settings
set(CMAKE_Fortran_FLAGS_DEBUG "${FF_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${FF_RELEASE}")
