
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine idf
  use modmain
  use modmpi
  use modxs
  use modfxcifc
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'idf'
  integer :: iq
  ! initialise universal variables
  call init0
  call init1
  ! save Gamma-point variables
  call xssave0
  ! initialize q-point set
  call init2xs
  call readfermi
  ! w-point parallelization for dielectric function
  call genparidxran('w',nwdf)
  write(unitout,'("Exchange-correlation kernel type :",i4)') fxctype
  write(unitout,'("  ",a)') trim(fxcdescr)
  ! loop over q-points
  do iq=1,nqpt
     ! call for q-point
     call idfq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): inverse dielectric &
          &function finished for q-point:',iq
  end do
  call barrier
  if ((procs.gt.1).and.(rank.eq.0)) then
     call idfgather
        write(unitout,'(a)') 'Info('//thisnam//'): inverse dielectric &
             &function gathered for q-point:'
  end if
  call barrier
  if (rank.eq.0) then
     do iq=1,nqpt
        ! call for q-point
        call xslinopt(iq)
        write(unitout,'(a,i8)') 'Info('//thisnam//'): TDDFT linear optics &
             &finished for q-point:',iq
     end do
  end if
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): TDDFT linear optics &
       &finished"
  call genfilname(setfilext=.true.)
end subroutine idf
