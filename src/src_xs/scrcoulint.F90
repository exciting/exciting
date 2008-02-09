
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrcoulint
  use modmain
  use modxs
  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  integer :: iknr,jknr,iq,j
  real(8) :: vklofft(3)
  integer :: ngridkt(3)
  logical :: nosymt,reducekt
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  ! map variables for screened Coulomb interaction
  call initbse
  nosym=nosymscr
  ! no symmetries implemented for screened Coulomb interaction
  reducek=.false.
  ! q-point set of screening corresponds to (k,kp)-pairs
  ngridk(:)=ngridq(:)
  vkloff(:)=vkloffbse(:)
  ! check number of empty states
  if (nemptyscr.lt.nempty) then
     write(*,*)
     write(*,'("Error(",a,"): too few empty states in screening eigenvector &
          &file - this is strange anyway (BSE/screening)",2i8)') trim(thisnam),&
          nempty,nemptyscr
     write(*,*)
     call terminate
  end if
  ! loop over non-reduced number of k-points
  do iknr=1,nkptnr     
     do iq=1,nqpt
        ! determine set of k-points as classes associated to the q-point
        


        !do j=1,nkq
        !jknr=


     end do
  end do
  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine scrcoulint
