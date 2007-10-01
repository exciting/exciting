
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! main routine for the EXCITING code
program main
use modmain
implicit none
! local variables
integer itask
! read input files
call readinput
! perform the appropriate task
do itask=1,ntasks
  task=tasks(itask)
  select case(task)
  case(-1)
    write(*,*)
    write(*,'("EXCITING version ",I1.1,".",I2.2,".",I3.3)') version
    write(*,*)
    stop
  case(0,1,2,3)
    call gndstate
  case(10)
    call dos
  case(15)
    call writelsj
  case(20,21)
    call bandstr
  case(25)
    call effmass
  case(31,32,33)
    call rhoplot
  case(41,42,43)
    call potplot
  case(51,52,53)
    call elfplot
  case(61,62,63,162)
    call wfplot
  case(72,73,82,83,142,143)
    call vecplot
  case(91,92,93)
    call dbxcplot
  case(100,101)
    call fermisurf
  case(110)
    call mossbauer
  case(115)
    call writeefg
  case(120)
    call writepmat
  case(121)
    call linopt
  case(200)
    call phonon
  case(210)
    call phdos
  case(220)
    call phdisp
  case(230)
    call writephn
  case default
    write(*,*)
    write(*,'("Error(main): task not defined : ",I8)') task
    write(*,*)
    stop
  end select
end do
stop
end program

!BOI
! !TITLE: The EXCITING Code Developers' Guide\\ Version 0.9.74
! !AUTHORS: J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl
! !AFFILIATION:
! !INTRODUCTION: Introduction
!   Welcome to the {\sf EXCITING} code developers' manual! If you're reading
!   this then you would presumably like to understand, change or contribute to
!   the code. This automatically generated manual is intended to make that task
!   as easy as possible.
!   \section{Contributing to {\sf EXCITING}}
!   Please bear in mind when writing code for the {\sf EXCITING} project that
!   it should be an exercise in physics and not software engineering. All code
!   should therefore be kept as simple and concise as possible, and above all it
!   should be easy for anyone to locate and follow the {\sf FORTRAN}
!   representation of the original mathematics. We would also appreciate the
!   following conventions being adhered to:
!   \begin{itemize}
!   \item Strict {\sf FORTRAN} 90 should be used. Features which are marked as
!    obsolescent in F90/95 should be avoided. These include assigned format
!    specifiers, labeled do-loops, computed goto statements and statement
!    functions.
!   \item Modules should be used in place of common blocks for declaring
!    global variables. Use the existing modules to declare new global variables.
!   \item Any code should be written in lower-case free form style, starting
!    from column one. Try and keep the length of each line to fewer than 80
!    characters using the \& character for line continuation.
!   \item Every function or subroutine, no matter how small, should be in its
!    own file named {\tt routine.f90}, where {\tt routine} is the function or
!    subroutine name. It is recommended that the routines are named so as to
!    make their purpose apparent from the name alone.
!   \item Use of {\tt implicit none} is mandatory. Remember also to define the
!    {\tt intent} of any passed arguments.
!   \item Local allocatable arrays must be deallocated on exit of the routine to
!    prevent memory leakage. Use of automatic arrays should be limited to arrays
!    of small size.
!   \item Every function or subroutine must be documented with the {\sf Protex}
!    source code documentation system. This should include a short \LaTeX\
!    description of the algorithms and methods involved. Equations which need to
!    be referenced should be labeled with {\tt routine\_1}, {\tt routine\_2}
!    etc. The authorship of each new piece of code or modification should be
!    indicated in the {\tt REVISION HISTORY} part of the header. See the
!    {\sf Protex} documentation for details.
!   \item Ensure as much as possible that a routine will terminate the program
!    when given improper input instead of continuing with erroneous results.
!    Specifically, functions should have a well-defined domain for which they
!    return accurate results. Input outside that domain should result in an
!    error message and termination.
!   \item Report errors prior to termination with a short description, for
!    example:
!    \begin{verbatim}
!     write(*,*)
!     write(*,'("Error(readparam): invalid spnst : ",I8)') spnst(is)
!     write(*,'(" for species ",I4)') is
!     write(*,*)
!    \end{verbatim}
!   \item Wherever possible, real numbers outputted as ASCII data should be
!    formatted with the {\tt G18.10} specifier.
!   \item Avoid redundant or repeated code: check to see if the routine you need
!    already exists, before writing a new one.
!   \item All reading in of ASCII data should be done in the subroutine
!    {\tt readinput}. For binary data, separate routines for reading and writing
!    should be used (for example, {\tt writestate} and {\tt readstate}).
!   \item Input file names should be in lowercase and have the extension
!    {\tt .in} . All output file names should be in uppercase with the extension
!    {\tt .OUT} .
!   \item All internal units should be atomic. Input and output units should be
!    atomic by default and clearly stated otherwise. Rydbergs should not be used
!    under any circumstances.
!   \end{itemize}
!   \section{Licensing}
!    Routines which constitute the main part of the code are released under the
!    GNU General Public License (GPL). Library routines are released under the
!    less restrictive GNU Lesser General Public License (LGPL). Both licenses
!    are contained in the file {\tt COPYING}. Any contribution to the code must
!    be licensed at the authors' discretion under either the GPL or LGPL.
!    Author(s) of the code retain the copyrights. Copyright and (L)GPL
!    information must be included at the beginning of every file, and no code
!    will be accepted without this.
!   \\\\
!   J. Kay Dewhurst\\
!   Edinburgh, 2006
!
!EOI

