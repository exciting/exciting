# Set compiler flags
# CMake will append CMAKE_Fortran_FLAGS with CMAKE_Fortran_FLAGS_BUILDTYPE
# CMAKE_Fortran_FLAGS_BUILDTYPE may also have predefined values, hence initialise it

# Function to check if the CPU is from Intel by parsing /proc/cpuinfo or 
# for macOS use sysctl
function(CheckIfIntelCPU result_var)
    if(APPLE)
        # macOS: Use sysctl to get CPU brand string
        execute_process(
            COMMAND sysctl -n machdep.cpu.brand_string
            OUTPUT_VARIABLE CPUINFO_CONTENTS
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    elseif(UNIX)
        # Linux: Read /proc/cpuinfo
        file(READ "/proc/cpuinfo" CPUINFO_CONTENTS)
    else()
	message(WARNING "CheckIfIntelCPU: Unsupported platform. Assuming no Intel CPU.")
        set(${result_var} FALSE PARENT_SCOPE)
        return()
    endif()

    # Check for "Intel" in CPU info
    if(CPUINFO_CONTENTS MATCHES "Intel")
        set(${result_var} TRUE PARENT_SCOPE)
    else()
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()


# Function to check if AVX-512 is supported by parsing /proc/cpuinfo
# for macOS sysctl does not provide information about the 
# supported instruction set for M1/M2/M3. Therefore, for macOS we return FALSE.
function(CheckForAVX512Support result_var)
    if(APPLE)
        # macOS (Intel or Apple Silicon): no reliable way to check AVX-512 via sysctl,
        # and Apple Silicon does not support AVX-512.
        set(${result_var} FALSE PARENT_SCOPE)
    elseif(UNIX)
        # Linux: read /proc/cpuinfo and check for avx512 flags
        file(READ "/proc/cpuinfo" CPUINFO_CONTENTS)
        if(CPUINFO_CONTENTS MATCHES "avx512")
            set(${result_var} TRUE PARENT_SCOPE)
        else()
            set(${result_var} FALSE PARENT_SCOPE)
        endif()
    else()
        message(WARNING "CheckForAVX512Support: Unsupported platform, assuming no AVX-512.")
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()

# Like CheckForAVX512Support but for AVX-2 
function(CheckForAVX2Support result_var)
    if(APPLE)
        # macOS does not expose CPU flags like AVX2 via sysctl on Apple Silicon or Intel.
        # Apple Silicon does NOT support AVX2.
        set(${result_var} FALSE PARENT_SCOPE)
    elseif(UNIX)
        # Linux: read /proc/cpuinfo and check for avx2 flag
        file(READ "/proc/cpuinfo" CPUINFO_CONTENTS)
        if(CPUINFO_CONTENTS MATCHES "avx2")
            set(${result_var} TRUE PARENT_SCOPE)
        else()
            set(${result_var} FALSE PARENT_SCOPE)
        endif()
    else()
        message(WARNING "CheckForAVX2Support: Unsupported platform, assuming no AVX2.")
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()

# Function to check if the Compiler supports pointer remapping
function(CheckPointerRemapping)
   file(WRITE ${CMAKE_BINARY_DIR}/config_tests/test_pointer_remapping.f90 "
       program check_pointer_remapping
           implicit none

           real, allocatable, target :: a(:)
           real, pointer, contiguous :: b(:,:)

           integer :: i, j

           allocate(a(100), source = 4.0)

           b(1:10,1:10) => a(:)
           b(0:9,0:9) => b(:,:)

           do i = 0,9
              do j = 0, 9
                 if (b(j,i) /= 4.0) error stop 'Bad value'
              end do
           end do

       end program check_pointer_remapping
   ")
   # Create a log file for the config test
    set(REMAPPING_CONFIG_TEST_LOG "${CMAKE_BINARY_DIR}/config_tests/remapping_test_log.txt")
    file(WRITE "${REMAPPING_CONFIG_TEST_LOG}" "Compilation Log:\n")
    set(REMAPPING_CONFIG_TEST_LOG_STR "")

    # Try to compile and run the Fortran test program during configuration.
    execute_process(
            COMMAND ${CMAKE_Fortran_COMPILER} -o test_pointer_remapping test_pointer_remapping.f90
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/config_tests
            RESULT_VARIABLE REMAPPING_TEST_COMPILE_RESULT
            OUTPUT_VARIABLE REMAPPING_CONFIG_TEST_LOG_STR
            ERROR_VARIABLE  REMAPPING_CONFIG_TEST_LOG_STR
    )

    # Append the compilation log to the log file
    file(APPEND "${REMAPPING_CONFIG_TEST_LOG}" "${REMAPPING_CONFIG_TEST_LOG_STR}\n")

    execute_process(
            COMMAND ./test_pointer_remapping
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/config_tests
            RESULT_VARIABLE REMAPPING_TEST_RUN_RESULT
            OUTPUT_VARIABLE REMAPPING_CONFIG_TEST_LOG_STR
            ERROR_VARIABLE  REMAPPING_CONFIG_TEST_LOG_STR
    )

    # Append the compilation log to the log file
    file(APPEND "${REMAPPING_CONFIG_TEST_LOG}" "${REMAPPING_CONFIG_TEST_LOG_STR}\n")

    if(NOT REMAPPING_TEST_COMPILE_RESULT EQUAL 0)
       message(FATAL_ERROR "Your Fortran compiler does not support pointer remapping. This is required.")
    endif()

    if(NOT REMAPPING_TEST_RUN_RESULT EQUAL 0)
        message(STATUS "Fortran compiler check (pointer remapping) : fails")
        message(STATUS "Pointer remapping is bugged; setting REMAPPING_BUG to activate workarounds.")
        add_compile_definitions(REMAPPING_BUG)
        if(AMD OR NVIDIA OR INTEL)
            message(FATAL_ERROR "GPU compilation requires pointer remapping.")
        endif()
    else()
        message(STATUS "Fortran compiler check (pointer remapping) : works")
    endif()

endfunction()

# Intel compilers might fail to detect non Intel machines. Here we perform few tests.
set(ARCH_MARCH "native")
CheckIfIntelCPU(IS_INTEL_CPU)
CheckForAVX512Support(AVX512_SUPPORTED)
CheckForAVX2Support(AVX2_SUPPORTED)
set(INTEL_CODE_NAME "Host" CACHE STRING "Code name for Intel processors")
message(STATUS "Compiling for an Intel CPU: ${IS_INTEL_CPU} (Intel code name: ${INTEL_CODE_NAME})")
message(STATUS "CPU instruction support for AVX2: ${AVX2_SUPPORTED}")
message(STATUS "CPU instruction support for AVX512: ${AVX512_SUPPORTED}")

# GCC
set(GCC_COMMON -march=${ARCH_MARCH} -fPIC -cpp -ffree-line-length-0 -fallow-argument-mismatch -fallow-invalid-boz )

set(GCC_DEBUG
     -g               # Generate symbols
     -fbacktrace      # symbolic stack traceback
     -fPIC            # Code independent position
     -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
     -finit-real=nan  #  All real scalars are initialised to NaN
     -fcheck=all      # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
     ${GCC_COMMON}
    )
# More debug flags to consider:
# -finit-integer=2147483647 -finit-real=snan \
# -frecord-gcc-switches -finit-character=42 -finit-logical=true -fdump-core -fstack-protector-all -pipe

set(GCC_RELEASE -g -O3 ${GCC_COMMON})  # Level 3 optimisation. Could also consider -fpack-derived

# Intel
set(INTEL_DEBUG
    -O0
    -g
    -fPIC
    -fpp                   # Use preprocessor
    -allow nofpp_comments  # Preprocessor option
    -g            # Generate symbols
    -traceback    # symbolic stack traceback
    -fp           # Disables the ebp register in optimizations and sets the ebp register to be used as the frame pointer.
    -check all    # Checks for all runtime failures.
    -check bounds # Generates code to perform runtime checks on array subscript and character substring expressions.
    -check-uninit #  Enables runtime checking for uninitialized variables.
    -ftrapuv      #  Set unassigned scalars as a very large integer or an invalid address
    -fpe3         # control over floating-point exception (divide by zero, overflow, invalid operation, underflow, denormalized number, positive infinity, negative infinity or a NaN)
    )

set(INTEL_RELEASE -O3 -g -fPIC -fpp -allow nofpp_comments -save-temps -no-wrap-margin -fp-model source )

# Cray compiler (add -hlist=m -h keepfiles to save the temporary files)
set(CRAY_RELEASE -O2 -ef -g -craype-verbose -e Z -dC -s real64 -s integer32 -fPIC -hipa0 -h flex_mp=strict -hnopattern)
set(CRAY_DEBUG   -O0 -ef -g -fsanitize=thread -craype-verbose -e Z -dC -s real64 -s integer32 -fPIC -hipa0 -h flex_mp=strict -hnopattern)

# Flang compiler
set(FLANG_RELEASE    -cpp -fno-fast-math -flto)
set(FLANG_DEBUG   -g -cpp -fno-fast-math -flto)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set(FF_DEBUG ${GCC_DEBUG})
   set(FF_RELEASE ${GCC_RELEASE})

   if(OMP)
      set(FF_DEBUG "${GCC_DEBUG} ${OpenMP_Fortran_FLAGS}")
      set(FF_RELEASE "${GCC_RELEASE} ${OpenMP_Fortran_FLAGS}")
   endif()

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

    # In case of Intel also add
    add_compile_definitions(IFORT)


   # It is little tricky for Intel compilers to work with non-Intel CPUs
   # Here we select the proper instruction set based on the supported instruction sets
   if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
      if(OMP)
         set(IFX_IPO "-ipo")
      else()
	 set(IFX_IPO "")
      endif()
      if (IS_INTEL_CPU AND AVX512_SUPPORTED)
	 set(INTEL_OPTIMIZATION "-x${INTEL_CODE_NAME} ${IFX_IPO} -mprefer-vector-width=512")
      elseif(IS_INTEL_CPU AND NOT AVX512_SUPPORTED)
         set(INTEL_OPTIMIZATION "-x${INTEL_CODE_NAME} ${IFX_IPO}")
      elseif(NOT IS_INTEL_CPU AND AVX512_SUPPORTED)
         set(INTEL_OPTIMIZATION "-march=x86-64-v4")
      elseif(NOT IS_INTEL_CPU AND NOT AVX512_SUPPORTED AND AVX2_SUPPORTED)
         set(INTEL_OPTIMIZATION "-march=x86-64-v3")
      else()
         set(INTEL_OPTIMIZATION "-march=x86-64-v2")
      endif()
   else()
      if (IS_INTEL_CPU)
         set(INTEL_OPTIMIZATION "-march=native")
      elseif(NOT IS_INTEL_CPU AND AVX2_SUPPORTED)
         set(INTEL_OPTIMIZATION "-march=core-avx2")
      else()
         set(INTEL_OPTIMIZATION "-march=sse4.2")
      endif()
   endif()

   set(FF_DEBUG ${INTEL_DEBUG} ${INTEL_OPTIMIZATION})
   set(FF_RELEASE ${INTEL_RELEASE} ${INTEL_OPTIMIZATION})

   if(OMP)
      set(FF_DEBUG "${INTEL_DEBUG} ${OpenMP_Fortran_FLAGS} -qmkl=parallel")
      set(FF_RELEASE "${INTEL_RELEASE} ${OpenMP_Fortran_FLAGS} -qmkl=parallel")
   endif()

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")

   set(FF_DEBUG ${CRAY_DEBUG})
   set(FF_RELEASE ${CRAY_RELEASE})

   set(ENV{FORMAT_TYPE_CHECKING} "RELAXED")

   if(OMP)
      set(FF_DEBUG "${CRAY_DEBUG} ${OpenMP_Fortran_FLAGS}")
      set(FF_RELEASE "${CRAY_RELEASE} ${OpenMP_Fortran_FLAGS}")
   endif()

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")

   set(FF_DEBUG ${FLANG_DEBUG})
   set(FF_RELEASE ${FLANG_RELEASE})

   if(OMP)
      set(FF_DEBUG "${FLANG_DEBUG} ${OpenMP_Fortran_FLAGS} -fopenmp-version=51")
      set(FF_RELEASE "${FLANG_RELEASE} ${OpenMP_Fortran_FLAGS} -fopenmp-version=51")
   endif()

else ()
     message(STATUS "Unrecognized compiler: only Intel, GNU, Cray, Flang (new) compilers are supported")
     message(FATAL_ERROR "Compiler id: ${CMAKE_Fortran_COMPILER_ID}")
endif()

string(REPLACE ";" " " FF_DEBUG "${FF_DEBUG}")
string(REPLACE ";" " " FF_RELEASE "${FF_RELEASE}")
string(REPLACE ";" " " STD_FFLAGS "${STD_FFLAGS}")

# Initialise BUILDTYPE flags so we completely define/control
# the compiler settings
# Note, these flags are GLOBALS and will apply to ALL libs/executables built by CMake
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} ${FF_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} ${FF_RELEASE}")

# Here test compiler bugs
CheckPointerRemapping()
