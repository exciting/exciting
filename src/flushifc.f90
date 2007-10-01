
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: flushifc
! !INTERFACE:
subroutine flushifc(fnum)
! !INPUT/OUTPUT PARAMETERS:
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
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
logical named_,opened_
character(32) action_,blank_,delim_,form_
character(1024) name_
inquire(fnum,action=action_,blank=blank_,delim=delim_,form=form_,name=name_, &
 named=named_,opened=opened_)
if ((adjustl(action_).ne.'WRITE').and.(adjustl(action_).ne.'READWRITE')) then
  write(*,*)
  write(*,'("Error(flushifc): unit ",I4," is read-only")') fnum
  write(*,*)
  stop
end if
if (adjustl(form_).ne.'FORMATTED') then
  write(*,*)
  write(*,'("Error(flushifc): unit ",I4," is not a formatted file")') fnum
  write(*,*)
  stop
end if
if (.not.named_) then
  write(*,*)
  write(*,'("Error(flushifc): unit ",I4," is not named")') fnum
  write(*,*)
  stop
end if
if (.not.opened_) then
  write(*,*)
  write(*,'("Error(flushifc): unit ",I4," is not connected")') fnum
  write(*,*)
  stop
end if
! close and re-open file
close(fnum)
open(fnum,action=action_,blank=blank_,delim=delim_,form=form_, &
 file=trim(name_),position='APPEND')
return
end subroutine
!EOC
