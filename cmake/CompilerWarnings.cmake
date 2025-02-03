# Add compiler warnings to project

  set(GCC_WARNINGS
      -Wall  # Contains aliasing, apersand, conversion, surprising, binding-type,
             # intrinsics-std, tabs, intrinsic-shadow, line-truncation,
	     # target-lifetime, integer-division, real-q-constant, unused and undefined-do-loop
      -Warray-bounds
      -Wcharacter-truncation # Warn when a character assignment will truncate the assigned string.
      -Wfunction-elimination # Warn if any calls to impure functions are eliminated by the optimizations
                             # enabled by the -ffrontend-optimize option
      -Wimplicit-interface   # Warn if a procedure is called without an explicit interface.
      -Wimplicit-procedure   # Warn if a procedure is called that has neither an explicit interface nor has been declared as EXTERNAL.
      -Wuse-without-only     # Warn if modules are used without only keyword
      -Wunderflow            # Produce a warning when numerical constant expressions are encountered, which yield an UNDERFLOW during compilation.
      -Wrealloc-lhs          # I assume this warns if the shape of the L.H.S. and R.H.S of an allocatable assignment differ.
                             # Fortran 2003: The L.H.S.is assumed to be allocated with the correct shape to hold the right-hand side. If it is not, incorrect behavior will occu. Intel's default is when the L.H.S of an assignment is an allocatable object, it should be reallocated to the shape of the right-hand side of the assignment before the assignment occurs. This is the current Fortran Standard definition. This feature may cause extra overhead at run time.
      -Wrealloc-lhs-all      #
      -Wfrontend-loop-interchange # Warn when using -ffrontend-loop-interchange for performing loop interchanges.
      -pedantic)    # Issue warnings for uses of extensions to Fortran i.e. C's #include
                    # This should be used in conjunction with -std=f95 to f2018

    # Removed warnings
    # Wimplicit-interface and Wimplicit-procedure are not usable with MPI.
    # See: https://stackoverflow.com/questions/41938663/why-does-mpich-3-0s-mpi-module-omit-explicit-interfaces-for-some-procedures

   # See https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/compiler-diagnostic-options/warn.html
   set(INTEL_WARNINGS
       -warn all     # Enables all warning messages except errors and stderrors
      )

   # This sets compiler output to: Error, Warning, Caution, Note, and Comment   
   set(CRAY_WARNINGS
      -m 0)

   set(FLANG_WARNINGS
       -Wall)

   option(WARNINGS_AS_ERRORS "Treat compiler warnings as error" FALSE)
   if (WARNINGS_AS_ERRORS)
     set(GCC_WARNINGS ${GCC_WARNINGS} -Werror)
     set(INTEL_WARNINGS ${INTEL_WARNINGS} "-warn errors")
   endif()

   # List to string
   string(REPLACE ";" " " GCC_WARNINGS "${GCC_WARNINGS}")
   string(REPLACE ";" " " INTEL_WARNINGS "${INTEL_WARNINGS}")

   if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set(COMPILER_WARNINGS ${GCC_WARNINGS})
   elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
      set(COMPILER_WARNINGS ${INTEL_WARNINGS})
   elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
      set(COMPILER_WARNINGS ${CRAY_WARNINGS})
   elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
      set(COMPILER_WARNINGS ${FLANG_WARNINGS})
   else ()
     message(SEND_ERROR "Warning flags have not been defined for this compiler: \
            ${CMAKE_Fortran_COMPILER_ID}")
   endif()
