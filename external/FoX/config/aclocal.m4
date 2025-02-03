dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006

dnl @synopsis TW_CHECK_ASSOCIATED_BUG([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler has a bug
dnl forbidding the use of ASSOCIATED in restricted expressions
dnl (versions of gfortran, pathscale & xlf have been found which do)
dnl 
dnl 
dnl @version 1.0
dnl @author <tow@uszla.me.uk>
dnl
AC_DEFUN([TW_CHECK_ASSOCIATED_BUG],[
dnl
AC_MSG_CHECKING([for ASSOCIATED in restricted expression bug])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      function test_bug(a) result(b)
      integer, pointer :: a
      integer, dimension(merge(1, 2, associated(a))) :: b
      b = 0
      end function test_bug
   ]]),
   [AC_MSG_RESULT([no])
    m4_default([$1],[:])
   ],
   [AC_MSG_RESULT([yes])
    m4_default([$2],[:] 
               [AC_MSG_ERROR([This Fortran compiler does not understand ASSOCIATED in restricted expressions.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006

dnl @synopsis TW_CHECK_FC_90([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with Fortran 90 (ISO/IEC-1539:1991)
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl It currently tests for:
dnl
dnl Modules
dnl Private 
dnl New-style variable declarations.
dnl
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_90],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_MSG_CHECKING([for Fortran 90 compliance])
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Module test_module

      Implicit None
      Private

      Contains

      Function test_function() Result(out)
      Integer :: out
      out = 0
      End Function test_function

      End Module test_module
   ]]),
   [AC_MSG_RESULT([yes])
    m4_default([$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_default([$2],
               [AC_MSG_ERROR([ A fully Fortran-90-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006

dnl @synopsis TW_CHECK_FC_95([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with Fortran 95 (ISO-IEC 1539-1:1997)
dnl 
dnl It currently tests for:
dnl
dnl Named End Interface
dnl Derived type initialization
dnl The Null() intrinsic
dnl The Forall statement 
dnl The Cpu_Time intrinsic
dnl Pure functions
dnl Elemental functions
dnl 
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_95],[
dnl
AC_MSG_CHECKING([for Fortran 95 compliance])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_f95

!      Interface test_interface
!      End Interface test_interface

      Type test_type
        Integer :: i = 1
      End Type test_type

      Integer, Pointer :: j => Null()

      Integer :: i
      Real :: a

      Forall (i=1:50)
      End Forall

      Call CPU_TIME(a)

      Contains

      Pure Integer Function test_pure()
        test_pure = 0
      End Function test_pure

      Elemental Integer Function test_elemental(in)
        Integer, Intent(In) :: in
        test_elemental = 0
      End Function test_elemental

      End Program test_f95
   ]]),
   [AC_MSG_RESULT([yes])
    m4_default([$1],[:])
   ],
   [AC_MSG_RESULT([no])
    m4_default([$2], 
               [AC_MSG_ERROR([A fully Fortran-95-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl @synopsis TW_CHECK_FC_FPP([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler supports
dnl cpp-like functionality when called on a suitable fixed-format file.
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_FPP],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_FC_SRCEXT(F)
dnl
AC_MSG_CHECKING([whether $FC has an integrated Fortran cpp-style preprocessor for fixed-form source])
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_cpp
#if 1
      Integer i
#else
      Integer j
#endif
      End Program
   ]]),
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A Fortran compiler with integrated cpp-style preprocessor for fixed-form source is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl @synopsis TW_CHECK_FC_FPP_90([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler supports
dnl cpp-like functionality when called on a suitable fixed-format file.
dnl If so, ACTION_IF_TRUE is performed; if not, ACTION_IF_FALSE
dnl 
dnl @version 1.0
dnl @author Toby White <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_FPP_90],[
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_FPP_SRCEXT(F90)
ac_ext=F90
dnl
AC_MSG_CHECKING([whether $FC has an integrated Fortran cpp-style preprocessor for free-form source])
dnl
FCFLAGS_save=$FCFLAGS
FCFLAGS="$FCFLAGS $FPPFLAGS_F90 $FCFLAGS_free"
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
Program test_cpp
#if 1
  Integer i
#else
  Integer j
#endif
End Program
   ]]),
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A Fortran compiler with integrated cpp-style preprocessor for free-from source is required.])])
   ]
)
AC_LANG_POP(Fortran)

FCFLAGS=$FCFLAGS_save
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl @synopsis TW_CHECK_FC_TR15580([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with the Fortran 95 Floating Point Exception Handling
dnl Extension, ISO TR15580.
dnl
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_TR15580],[
dnl
AC_MSG_CHECKING([$FC for compliance to the Floating Point Exception Handling Extension])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_tr15580

      Use, Intrinsic :: IEEE_Arithmetic
      Use, Intrinsic :: IEEE_Exceptions
      Use, Intrinsic :: IEEE_Features

      End Program test_tr15580
   ]]),   
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A fully TR15580-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl @synopsis TW_CHECK_FC_TR15581([ACTION_IF_TRUE],[ACTION_IF_FALSE])
dnl
dnl Checks whether the currently selected Fortran compiler is fully
dnl compliant with the Fortran 95 Enhanced Datatype Facilities
dnl Extension, ISO TR15581.
dnl
dnl @version 1.0
dnl @author <tow21@cam.ac.uk>
dnl
AC_DEFUN([TW_CHECK_FC_TR15581],[
dnl
AC_MSG_CHECKING([$FC for compliance to the Enhanced Datatype Facilities Extension])
dnl
AC_LANG_PUSH(Fortran)
dnl
AC_COMPILE_IFELSE(
dnl The program is written in fixed-form source to avoid worrying
dnl about filename extensions.
  AC_LANG_SOURCE([[
      Program test_tr15581

      Type test_type
        Integer, Allocatable :: array(:)
      End Type test_type

      End Program test_tr15581

      Function test_function
        Integer, Allocatable :: test_function(:)

        Allocate(test_function(5))
      End Function test_function

      Subroutine test_subroutine(array)
        Integer, Allocatable :: array(:)

        Allocate(array(5))
      End Subroutine test_subroutine
   ]]),   
   [AC_MSG_RESULT([yes])
    m4_ifval([$1],[$1],[])
   ],
   [AC_MSG_RESULT([no])
    m4_ifval([$2],[$2],
                  [AC_MSG_ERROR([A fully TR15581-compliant compiler is required.])])
   ]
)
AC_LANG_POP(Fortran)
dnl
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl Check how to get at the abort intrinsic.

AC_DEFUN([_TW_TRY_ABORT_BARE],
[
      PROGRAM TESTABORT
      CALL ABORT
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_NAG],
[
      PROGRAM TESTABORT
      USE F90_UNIX_PROC, ONLY:ABORT
      CALL ABORT
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_INTEL],
[
      PROGRAM TESTABORT
      CALL ABORT("")
      END PROGRAM TESTABORT
])
AC_DEFUN([_TW_TRY_ABORT_XLF],
[
      PROGRAM TESTABORT
      CALL ABORT_
      END PROGRAM TESTABORT
])

AC_DEFUN([TW_FC_CHECK_ABORT], [
AC_REQUIRE([AC_PROG_FC])dnl
dnl
AC_MSG_CHECKING([how to compile a call to ABORT])
dnl
dnl Try first with nothing
dnl
tw_abort_ok=no
dnl
dnl First check with one arg (this will fail if no args are necessary; testing
dnl in the opposite order will succeed when it shouldnt)
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_INTEL])],
    [tw_abort_ok=yes; tw_method="with argument";echo "#define FC_HAVE_ABORT 1" >> confdefs.h; echo "#define FC_ABORT_ARG 1" >> confdefs.h],
    [])
fi
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_XLF])],
    [tw_abort_ok=yes; tw_method="with underscore";echo "#define FC_HAVE_ABORT 1" >> confdefs.h; echo "#define FC_ABORT_UNDERSCORE 1" >> confdefs.h],
    [])
fi
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_BARE])],
    [tw_abort_ok=yes; tw_method=default;echo "#define FC_HAVE_ABORT 1" >> confdefs.h],
    [])
fi
dnl
if test $tw_abort_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_NAG])],
    [tw_abort_ok=yes; tw_method="with f90_unix_proc";echo "#define FC_HAVE_ABORT 1" >> confdefs.h],
    [])
fi
dnl
dnl Cant get it to compile alone - need a compiler flag.
dnl Now try with -Vaxlib for intel:
dnl
if test $tw_abort_ok = no; then
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS -Vaxlib"
   AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_ABORT_BARE])],
    [tw_abort_ok=yes; tw_method="with -Vaxlib";echo "#define FC_HAVE_ABORT 1" >> confdefs.h],
    [])
   if test $tw_abort_ok = no; then
      LDFLAGS=$save_LDFLAGS
   fi
fi
AC_MSG_RESULT([$tw_method])
dnl
AS_IF([test $tw_abort_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile call to ABORT ])])]
     )
dnl
])# TW_FC_CHECK_ABORT
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     


# TW_FC_CHECK_EOL()
# -------------------
#
# This macro checks what kind of line ending the Fortran compiler
# writes out, and expects to read in. (Checking reading explicitly really doesnt work)


AC_DEFUN([TW_FC_CHECK_EOL], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG_PUSH([Fortran])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
AC_MSG_CHECKING([for EOR character used by $FC])
AC_RUN_IFELSE([
       program output_eol
       open(unit=10, file="conf_eol.out")
       write(10,"(a)") ""
       end program
],
[
if ! test -f conf_eol.out; then
AC_MSG_ERROR([Could not find test output])
elif od -b conf_eol.out | grep 5015 >/dev/null; then
  ac_cv_FC_check_output_eol=CRLF
elif od -b conf_eol.out | grep 15 >/dev/null; then
  ac_cv_FC_check_output_eol=CR
elif od -b conf_eol.out | grep 12 >/dev/null; then
  ac_cv_FC_check_output_eol=LF
else
  ac_cv_FC_check_output_eol=UNKNOWN
fi
rm -rf conf_eol.out
if test $ac_cv_FC_check_output_eol = UNKNOWN; then
  AC_MSG_ERROR([Could not determine line-ending convention])
fi
], 
[AC_MSG_ERROR([Could not execute compiled program])],
[ac_cv_FC_check_output_eol=EOL_CR] dnl take a wild guess at Unix if x-compiling
)
AC_MSG_RESULT([$ac_cv_FC_check_output_eol])
])
dnl output EOL
dnl
AC_DEFUN([_TW_FC_CHECK_INPUT_EOL], [
AC_MSG_CHECKING([for input EOL character produced by $FC])
touch conf_empty.txt
# Surely there must be a better way to create a CRLF file than this!
echo 12 | tr -d "\012\015" | tr 12 "\015\012" > conf_crlf.txt
AC_RUN_IFELSE([
       program input_eol
       integer :: i, io_eof
       integer :: s1, result
       character :: c
       open(unit=10, file="conf_empty.txt")
       open(unit=11, file="conf_crlf.txt")
       open(unit=12, file="conf_result.txt")
       ! Pick up eof first
       read(10, "(a1)", iostat=io_eof) c
       i = 0
       n = 1
       s = 0
       result = 0
       read(11, "(a1)", iostat=i) c
       ! If we are on an LF-EOL machine,
       ! then we should get CR followed by EOR
       if (i==0) then
         if (iachar(c)==13) then
           s1 = 13
         elseif (iachar(c)==32) then
           ! some compilers translate it into a space, unhelpfully
           s1 = 32
         else
           write(12, *) "UNKNOWN"
           stop
         endif
       else
         s1 = -1 ! End of Record, we assume
       endif
       read(11, "(a1)", iostat=i) c
       if (i==0) then
         if (iachar(c)==10.and.s1==-1) then
	   ! Sequence was EOR, LF, therefore EOR=CR.
           ! Next must be EOF
           read(11, "(a1)", iostat=i)
           if (i==io_eof) result = 1 ! EOR_CR
         endif
       elseif (i==io_eof) then
         if (s1==-1) then
           ! Sequence was EOR, EOF, therefore EOR=CRLF
           result = 2 ! EOR_CRLF
         elseif (s1==32) then
           ! Sequence was SPACE, EOF. Empirically, this seems to happen on PPC Macs, so:
           result = 3 ! EOR_LF
         endif
       elseif (s1==13) then
         ! We assume this non-zero iostat is EOR
         ! Sequence was CR, EOR, therefore EOR=LF
         ! Next must be EOF
         read(11, "(a1)", iostat=i)
         if (i==io_eof) result = 3 ! EOR_LF
       endif
       select case(result)
       case (1)
         write(12,*) "CR"
       case (2)
         write(12,*) "CRLF"
       case (3)
         write(12,*) "LF"
       case default
         write(12,*) "UNKNOWN"
       end select
       end program
],
[
rm -f conf_empty.txt conf_crlf.txt
if ! test -f conf_result.txt; then
  AC_MSG_ERROR([Could not find test output])
elif grep CRLF conf_result.txt >/dev/null; then
  ac_cv_FC_check_input_eol=CRLF
elif grep CR conf_result.txt >/dev/null; then
  ac_cv_FC_check_input_eol=CR
elif grep LF conf_result.txt >/dev/null; then
  ac_cv_FC_check_input_eol=LF
else
  ac_cv_FC_check_input_eol=UNKNOWN
fi
dnl rm -f conf_result.txt
if test $ac_cv_FC_check_input_eol = UNKNOWN; then
  AC_MSG_ERROR([Could not determine input line-ending convention])
fi
], 
[AC_MSG_ERROR([Could not execute compiled program])],
[ac_cv_FC_check_input_eol=CR] dnl take a wild guess at Unix if x-compiling
)
AC_MSG_RESULT([$ac_cv_FC_check_output_eol])
dnl check input_eol
dnl
FCFLAGS="$FCFLAGS_save"
])])dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl Check how to get at the flush intrinsic.

AC_DEFUN([_TW_TRY_FLUSH_BARE],
[
      PROGRAM TESTFLUSH
      PRINT*
      CALL FLUSH(5)
      END PROGRAM TESTFLUSH
])
AC_DEFUN([_TW_TRY_FLUSH_NAG],
[
      PROGRAM TESTFLUSH
      USE F90_UNIX_IO, ONLY:FLUSH
      PRINT*
      CALL FLUSH(5)
      END PROGRAM TESTFLUSH
])
AC_DEFUN([_TW_TRY_FLUSH_XLF],
[
      PROGRAM TESTFLUSH
      PRINT*
      CALL FLUSH_(5)
      END PROGRAM TESTFLUSH
])

AC_DEFUN([TW_FC_CHECK_FLUSH], [
AC_REQUIRE([AC_PROG_FC])dnl
dnl
AC_MSG_CHECKING([how to compile a call to FLUSH])
dnl
dnl Try first with nothing
dnl
tw_flush_ok=no
dnl
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_BARE])],
    [tw_flush_ok=yes; TW_FLUSH=bare;tw_method=default;echo "#define FC_HAVE_FLUSH 1" >> confdefs.h],
    [])
if test $tw_flush_ok = no; then
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS -Vaxlib"
   AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_BARE])],
    [tw_flush_ok=yes; TW_FLUSH=INTEL;tw_method="with -Vaxlib";echo "$define FC_HAVE_FLUSH 1" >> confdefs.h],
    [])
   if test $tw_flush_ok = no; then
      LDFLAGS=$save_LDFLAGS
   fi
fi
if test $tw_flush_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_NAG])],
    [tw_flush_ok=yes; TW_FLUSH=NAG;tw_method="with f90_unix_io";echo "#define FC_HAVE_FLUSH 1" >> confdefs.h],
    [])
fi
if test $tw_flush_ok = no; then
  AC_LINK_IFELSE(
   [AC_LANG_SOURCE([_TW_TRY_FLUSH_XLF])],
    [tw_flush_ok=yes; TW_FLUSH=XLF;tw_method="with underscore"],
    [])
fi
AC_MSG_RESULT([$tw_method])
dnl
AS_IF([test $tw_flush_ok = yes],
      [$1],
      [m4_default([$2],[AC_MSG_ERROR([Cannot compile FLUSH statement])])]
     )
dnl
])# TW_FC_CHECK_FLUSH
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl A macro to determine which compiler is being used, in order that
dnl different flags can be set

AC_DEFUN([TW_FC_ID], [
AC_REQUIRE([AC_PROG_FC])

FC_ID=

dnl Firstly go by compiler name.

case $FC in 
   
   g77*)
      FC_ID=G77
      ;;

   g95*)
      FC_ID=G95
      ;;

   gfortran*)
      FC_ID=Gfortran
      ;;

   if*)
      FC_ID=Intel
      ;;

   lf9*)
      FC_ID=Lahey
      ;;
   
   path*)
      FC_ID=Pathscale
      ;;

   pgf*)
      FC_ID=Portland
      ;;

   xlf*)
      FC_ID=Xlf 

esac

dnl then try and disambiguate all f77, f90, and f95 types.
dnl We should have a choice between
dnl nag. absoft. sun. sgi. digital. hp. cray. ...?

if test x$FC_ID = x; then
   tw_fc_v_output=$($FC -V 2>&1 )
   if test $?; then
      case $tw_fc_v_output in
         *NAG*)
            FC_ID=Nag
            ;;
         *Sun*)
            FC_ID=Sun # there's more than one compiler here ...
            ;;
         *Absoft*)
            FC_ID=Absoft # there's more than one compiler here ...
            ;;
      esac
   fi
fi
 if test x$FC_ID = x; then
   tw_fc_v_output=$($FC -version 2>&1)
   if test $?; then
      case $tw_fc_v_output in
         *Compaq*)
            FC_ID=Digital
            ;;
         *Digital*)
            FC_ID=Digital
            ;;
         *SGI*)
            FC_ID=SGI
            ;;
      esac
   fi
fi   
   
AS_IF([test x$FC_ID != x],
      [AC_MSG_NOTICE([$FC seems to be a $FC_ID compiler])],
      [FC_ID=unknown; AC_MSG_NOTICE([Could not determine type of compiler])])

dnl for more fun, try and get the version number now ...


])# TW_FC_ID
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl A macro to set various compiler-dependent things that can't be sensibly
dnl deduced.

AC_DEFUN([TW_FC_ID_FLAGS], [
AC_REQUIRE([TW_FC_ID])

case $FC_ID in

  Absoft)
     FFLAGS_DEBUG="-et -g -Rb -Rc -Rp -Rs"
     ;;

  Digital)
     FFLAGS_FAST=-O2
     FFLAGS_DEBUG="-g -Rabc -ei"
     ;;

  G77)
     ;;

  G95)
    FFLAGS_FAST=-O2
    FFLAGS_DEBUG="-ggdb3 -ftrace=full -fbounds-check -flogical=false -freal=nan -fpointer=invalid -finteger=-1"
    ;;

  Gfortran)
     ;;
 
  Intel)
     FFLAGS_DEBUG="-C -g -inline_debug_info"
     ;;

  Lahey)
     FFLAGS_DEBUG="--chk aesux --chkglobal -g --trace"
     FFLAGS_FAST="-O --warn --quiet --tpp --ntrace"
     ;;

  Nag)
     FFLAGS_DEBUG="-C=all -g -gline -nan"
     echo "#define __NAG__ 1" >> confdefs.h
     ;;
  
  Pathscale)
     ;;

  Portland)
     FFLAGS_DEBUG="-g -Mbounds"
     FFLAGS_FAST="-fast"
     echo "#define PGF90 1" >> confdefs.h
     ;;

  SGI)
     FFLAGS_DEBUG="-g -O0"
     FFLAGS_FAST="-O3 -OPT:Olimit=0"
     ;;

  Sun)
     FFLAGS_DEBUG="-C -g"
     FFLAGS_FAST="-fast"
     ;;

  Xlf)
     FFLAGS_DEBUG="-g -C -qinitauto -qsave -qmaxmem=16000 -qnolm"
     FFLAGS_FAST="-O3 -qarch=auto -qtune=auto -qcache=auto -qnolm"
     ;;

esac

AC_SUBST(FFLAGS_MPI)

])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     


# AC_FC_REAL_KIND([KIND_DECLARATION], [VARIABLE_SUFFIX], [ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro checks what integer is produced by the kind 
# declaration KIND_DECLARATION. This integer is placed in 
# AC_FC_KIND_<VARIABLE_SUFFIX>. If we successfully find a
# kind integer, ACTION_IF_SUCCESS is performed; otherwise
# ACTION_IF_FAIL.

AC_DEFUN([AC_FC_REAL_KIND], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for kind number produced by $1], 
                ac_cv_fc_real_kind_[]$2[],
[dnl
AC_LANG_PUSH([Fortran])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
ac_fc_kind_test=1
ac_fc_kind_found=no
while test $ac_fc_kind_test -lt 100
do
  cat > conftest.$ac_ext << ACEOF
dnl The program below will fail to compile if 
dnl sp != mysp; ie if the kind produced by the 
dnl supplied kind declaration ($1) is not the same
dnl same as $ac_fc_kind_test. This is because Fortran
dnl pointers & targets must be of the same kind. 
dnl All conforming compilers must fail to compile the
dnl subroutine otherwise.
dnl
dnl This approach is taken since it enables us to
dnl test for kind numbers at compile time rather
dnl than run time, which means the macro will support
dnl crosss-compilation.
dnl
dnl However, note that kind numbers can theoretically
dnl be anything from 1 to the largest default integer
dnl supported by the compiler. Here we only test up to
dnl 99, which is more than enough on all compilers tried
dnl so far
dnl
      subroutine kind_explorer
      integer, parameter :: sp = $1
      integer, parameter :: mysp = $ac_fc_kind_test
      real(kind=sp), target :: x
      real(kind=mysp), pointer :: y
      y=>x
      end subroutine kind_explorer
ACEOF
dnl
  if eval $ac_compile 2>&5
  then
    ac_fc_kind_found=yes
    break
  fi
  ac_fc_kind_test=`expr $ac_fc_kind_test + 1`
done
if test "$ac_fc_kind_found" = yes; then
  ac_cv_fc_real_kind_[]$2[]=$ac_fc_kind_test
else 
  ac_cv_fc_real_kind_[]$2[]=none
fi

FCFLAGS="$FCFLAGS_save"
AC_LANG_POP([Fortran])
])
AS_IF([test $ac_cv_fc_real_kind_[]$2[] != no],
      [ac_fc_real_kind_[]$2[]=$ac_cv_fc_real_kind_[]$2[]; m4_default([$3],[:])],
      [m4_default([$4],[AC_MSG_ERROR([Could not find Fortran real kind number for $1])])]
     )
]) # AC_FC_REAL_KIND

# AC_FC_INT_KIND([KIND_DECLARATION], [VARIABLE_SUFFIX], [ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro checks what integer is produced by the kind 
# declaration KIND_DECLARATION. This integer is placed in 
# AC_FC_KIND_<VARIABLE_SUFFIX>. If we successfully find a
# kind integer, ACTION_IF_SUCCESS is performed; otherwise
# ACTION_IF_FAIL.

AC_DEFUN([AC_FC_INT_KIND], [dnl
AC_REQUIRE([AC_PROG_FC])
AC_CACHE_CHECK([for kind number produced by $1], 
                ac_cv_fc_int_kind_[]$2[],
[dnl
AC_LANG_PUSH([Fortran])
FCFLAGS_save="$FCFLAGS"
FCFLAGS="$FCFLAGS $FCFLAGS_free_f90"
ac_fc_kind_test=1
ac_fc_kind_found=no
while test $ac_fc_kind_test -lt 100
do
  cat > conftest.f90 << ACEOF
dnl The program below will fail to compile if 
dnl sp != mysp; ie if the kind produced by the 
dnl supplied kind declaration ($1) is not the same
dnl same as $ac_fc_kind_test. This is because Fortran
dnl pointers & targets must be of the same kind. 
dnl All conforming compilers must fail to compile the
dnl subroutine otherwise.
dnl
dnl This approach is taken since it enables us to
dnl test for kind numbers at compile time rather
dnl than run time, which means the macro will support
dnl crosss-compilation.
dnl
dnl However, note that kind numbers can theoretically
dnl be anything from 1 to the largest default integer
dnl supported by the compiler. Here we only test up to
dnl 99, which is more than enough on all compilers tried
dnl so far
dnl
      subroutine kind_explorer
      integer, parameter :: sp = $1
      integer, parameter :: mysp = $ac_fc_kind_test
      integer(kind=sp), target :: x
      integer(kind=mysp), pointer :: y
      y=>x
      end subroutine kind_explorer
ACEOF
dnl
  if eval $ac_compile 2>&5
  then
    ac_fc_kind_found=yes
    break
  fi
  ac_fc_kind_test=`expr $ac_fc_kind_test + 1`
done
if test "$ac_fc_kind_found" = yes; then
  ac_cv_fc_int_kind_[]$2[]=$ac_fc_kind_test
else 
  ac_cv_fc_int_kind_[]$2[]=none
fi

FCFLAGS="$FCFLAGS_save"
AC_LANG_POP([Fortran])
])
AS_IF([test $ac_cv_fc_int_kind_[]$2[] != no],
      [ac_fc_int_kind_[]$2[]=$ac_cv_fc_int_kind_[]$2[]; m4_default([$3],[:])],
      [m4_default([$4],[AC_MSG_ERROR([Could not find Fortran integer kind number for $1])])]
     )
]) # AC_FC_INT_KIND


# AC_FC_GET_REAL_KINDS([ACTION_IF_SUCCESS], [ACTION_IF_FAIL])
# -------------------
#
# This macro attempts to find the Fortran compiler's kinds
# for the following four types of real number:
#  Compiler default (single) precision
#  Compiler double precision
#  IEEE single precision
#  IEEE double precision
# The first two are guaranteed to exist; the second two may
# or may not.
# If all 4 are succesfully found,. ACTION_IF_SUCCESS is
# performed.
# Otherwise, ACTION_IF_FAIL is performed
#
AC_DEFUN([AC_FC_GET_REAL_KINDS], [dnl
AC_REQUIRE([AC_PROG_FC])

ac_fc_got_kinds=yes

AC_FC_REAL_KIND([[kind(1.0)]], [sp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[kind(1.0d0)]], [dp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[selected_real_kind(6,34)]], [ieee_sp], 
                [], [ac_got_kinds=no])
AC_FC_REAL_KIND([[selected_real_kind(15,300)]], [ieee_dp], 
                [], [ac_got_kinds=no])

AS_IF([test $ac_fc_got_kinds != no],
      [m4_default([$1],[:])],
      [m4_default([$2],[AC_MSG_ERROR([Could not find all Fortran real kinds])])]
      )
]) dnl AC_FC_GET_REAL_KINDS

dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006     

dnl autoconf macros for detecting NetCDF (fortan implementation only)
dnl
AC_DEFUN([_TW_TRY_NETCDF], [
ac_ext=f
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([[
      PROGRAM NETCDF
      CALL NF_CLOSE()
      END PROGRAM
   ]])],
   [m4_ifval([$1],[$1],[])],
   [m4_ifval([$2],[$2],[])]
)
])
dnl
AC_DEFUN([TW_PATH_NETCDF],[
tw_netcdf_ok=no
dnl
AC_LANG_PUSH([Fortran])
dnl
case $with_netcdf in
  yes | "") ;;
  no) tw_netcdf_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) NETCDF_LIBS="$with_netcdf" ;;
   *) NETCDF_LIBS="-l$with_netcdf" ;;
esac
dnl
if test $tw_netcdf_ok != disable; then
  if test "x$NETCDF_LIBS" = x; then
     NETCDF_LIBS="-lnetcdf"
  fi
dnl
  AC_MSG_CHECKING([for netcdf])
  save_LIBS="$LIBS"
  LIBS="$LIBS $NETCDF_LIBS"
  _TW_TRY_NETCDF([tw_netcdf_ok=yes],[])
  AC_MSG_RESULT([$tw_netcdf_ok])
  LIBS="$save_LIBS"
fi
dnl
dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test $tw_netcdf_ok = yes],
      [ifelse([$1],,AC_DEFINE(HAVE_NETCDF,1,[Define if you have NetCDF library.]),[$1])],
      [NETCDF_LIBS="";tw_netcdf_ok=no;$2])
AC_LANG_POP
])
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl As a special exception, the Free Software Foundation gives unlimited
dnl permission to copy, distribute and modify the configure scripts that
dnl are the output of Autoconf.  You need not follow the terms of the GNU
dnl General Public License when using or distributing such scripts, even
dnl though portions of the text of Autoconf appear in them.  The GNU
dnl General Public License (GPL) does govern all other use of the material
dnl that constitutes the Autoconf program.
dnl
dnl Certain portions of the Autoconf source text are designed to be copied
dnl (in certain cases, depending on the input) into the output of
dnl Autoconf.  We call these the "data" portions.  The rest of the Autoconf
dnl source text consists of comments plus executable code that decides which
dnl of the data portions to output in any given case.  We call these
dnl comments and executable code the "non-data" portions.  Autoconf never
dnl copies any of the non-data portions into its output.
dnl
dnl This special exception to the GPL applies to versions of Autoconf
dnl released by the Free Software Foundation.  When you make and
dnl distribute a modified version of Autoconf, you may extend this special
dnl exception to the GPL to apply to your modified version as well, *unless*
dnl your modified version has the potential to copy into its output some
dnl of the text that was the non-data portion of the version that you started
dnl with.  (In other words, unless your change moves or copies text from
dnl the non-data portions to the data portions.)  If your modification has
dnl such potential, you must delete any notice of this special exception
dnl to the GPL from your modified version.
dnl
dnl Copyright Toby White <tow21@cam.ac.uk>  2004-2006 

dnl Macro to check for cygpath PATH-substitution utility.
dnl compilation with BLACS. Only works with Fortran at the moment.

AC_DEFUN([TW_PROG_CYGPATH_W], [

 # test whether we have cygpath 	 
 if test -z "$CYGPATH_W"; then 	 
   if (cygpath --version) >/dev/null 2>/dev/null; then 	 
     CYGPATH_W='cygpath -w' 	 
   else 	 
     CYGPATH_W=echo 	 
   fi 	 
 fi
 AC_SUBST([CYGPATH_W])
])
