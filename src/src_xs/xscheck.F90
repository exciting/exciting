
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xscheck
  use modmain, only: task,reducek,spinsprl,spinpol,version
  use modxs
  use modtetra
  implicit none
  ! local variables
  character(*), parameter :: thisnam='xscheck'
  integer errc, warnc
  errc = 0
  warnc = 0
  ! no spin-spirals
  if ( spinsprl ) then
     write(*,*) 'Error('//thisnam//'): not working for spin-spirals'
     errc = errc + 1
  end if
  ! warn for spin polarized calculations
  if ( spinpol ) then
     write(*,*) 'Warning('//thisnam//'): calculation is spin-polarized'
     warnc = warnc + 1
  end if
 !!$ ! type of response functions
 !!$ if (rsptype.eq.'tord') then
 !!$    write(*,'(a,2i8)') 'Error('//thisnam//'): code limitation - only &
 !!$         &retarded response functions implemented'
 !!$    errc = errc + 1
 !!$ end if
  ! tetrahedron method not implemented for analytic continuation
  if (tetra.and.acont) then
     write(*,*) 'Error('//thisnam//'): tetrahedron method does not work &
          &together with analytic continuation'
     errc = errc + 1
  end if
  ! stop on errors
  if ( errc .gt. 0 ) then
     write(*,*) '  Errors occurred - abort'
     call terminate
  end if
  ! warn on warnings
  if ( warnc .gt. 0 ) then
     write(*,*) '  Warnings occurred:', warnc
  end if
  tscreen=.false.
  if ((task.ge.400).and.(task.le.499)) tscreen=.true.
end subroutine xscheck
