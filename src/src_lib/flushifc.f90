!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: flushifc
! !INTERFACE:
!
!
Subroutine flushifc (fnum)
! !INPUT/OUTPUT PARAMETERS:
      Use modinput
!   fnum : unit specifier for file (in,integer)
! !DESCRIPTION:
!   Interface to the Fortran {\tt flush} statement. Some compilers do not
!   support the {\tt flush} command, which is very useful for keeping small
!   formatted files up-to-date on the disk. The routine implimented below is a
!   machine-independent emulation of {\tt flush}, but may be replaced with the
!   intrinsic command if preferred.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
! local variables
      Logical :: named_, opened_
      Character (32) :: action_, blank_, delim_, form_
      Character (1024) :: name_
      Inquire (fnum, Action=action_, Blank=blank_, Delim=delim_, &
     & Form=form_, Name=name_, Named=named_, Opened=opened_)
!
      If ((adjustl(action_) .Ne. 'WRITE') .And. (adjustl(action_) .Ne. &
     & 'READWRITE')) Then
         Write (*,*)
         Write (*, '("Error(flushifc): unit ", I4, " is read-only")') &
        & fnum
         Write (*,*)
  !stop
      End If
      If (adjustl(form_) .Ne. 'FORMATTED') Then
         Write (*,*)
         Write (*, '("Error(flushifc): unit ", I4, " is not a formatted&
        & file")') fnum
         Write (*,*)
  !stop
      End If
      If ( .Not. named_) Then
         Write (*,*)
         Write (*, '("Error(flushifc): unit ", I4, " is not named")') &
        & fnum
         Write (*,*)
  !stop
      End If
      If ( .Not. opened_) Then
         Write (*,*)
         Write (*, '("Error(flushifc): unit ", I4, " is not connected")&
        &') fnum
         Write (*,*)
  !stop
      End If
! close and re-open file
      Close (fnum)
      Open (fnum, Action=action_, Blank=blank_, Delim=delim_, &
     & Form=form_, File=trim(name_), Position='APPEND')
!
      Return
End Subroutine
!EOC
