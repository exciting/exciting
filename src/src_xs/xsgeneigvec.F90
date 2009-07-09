


! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine xsgeneigvec
  use modmain
use modinput
  use modmpi
  use modxs
  use m_writegqpts
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='xsgeneigvec'
  logical, parameter :: tq0ev=.true.
  real(8) :: vqlt(3)
  integer :: iq, qi, qf
  logical, external :: tqgamma
  ! initialize universal variables
  call init0
  call init1
  call init2
  ! SCF allready parallelized for k-point set
  qi=1
  qf=nqpt
  ! add extra q-point for if files for q=0 are to be calculated

  ! if first Q-point is Gamma-point we copy files
  if (tqgamma(1)) qi=1
  if (tscreen) then
     qi=0
     qf=0
  end if
  ! write q-points
  if (rank.eq.0) call writeqpts
  ! calculate eigenvectors for each q-point (k+q point set)
  do iq=qi, qf
     if (.not.tscreen) &
	  call genfilname(iqmt=max(0, iq), setfilext=.true.)
     vqlt(:)=0.d0
     if ((iq.ne.0).and.(rank.eq.0)) then
	vqlt(:)=vql(:, iq)
	call writegqpts(iq, filext)
     end if
     ! write eigenvectors, -values, occupancies and contracted MT coefficients
     call writeevec(vqlt, qvkloff(1, iq), filext)
     if (.not.tscreen) then
	write(unitout, '("Info(", a, "): eigenvectors generated for Q-point (iq, vql below)")') thisnam
	write(unitout, '(i6, 3g18.10)') iq, vqlt(:)
     else
	write(unitout, '(a)') 'Info('//thisnam//'): eigenvectors &
	     &generated for associated(input%xs%screening)/screened interaction:'
     end if
     if (rank.eq.0) then
        ! safely remove unnecessary files
	call filedel('EQATOMS'//trim(filext))
	call filedel('EVALCORE'//trim(filext))
	call filedel('FERMIDOS'//trim(filext))
	call filedel('GEOMETRY'//trim(filext))
	call filedel('associated(input%structure%symmetries%lattice)'//trim(filext))
	call filedel('IADIST'//trim(filext))
	call filedel('LINENGY'//trim(filext))
	call filedel('SYMCRYS'//trim(filext))
	call filedel('SYMLAT'//trim(filext))
	call filedel('SYMSITE'//trim(filext))
	call filedel('TOTENERGY'//trim(filext))
	call filedel('EVALFV'//trim(filext))
	call filedel('RMSDVEFF'//trim(filext))
	call filedel('SYMGENR'//trim(filext))
	call filedel('SYMINV'//trim(filext))
	call filedel('SYMMULT'//trim(filext))
	call filedel('SYMMULT_TABLE'//trim(filext))
	call filedel('SYMT2'//trim(filext))
     end if
     ! end loop over q-points
  end do
  if ((rank.eq.0).and.tqgamma(1).and.(.not.tscreen)) then
     write(unitout, '("Info(", a, "): first Q-point is Gamma-point - copying &
       &relevant files")') thisnam
     ! write files again one by one
     call copyfilesq0
  end if
  call barrier
  write(unitout, '("Info(", a, "): generation of eigenvectors finished")') thisnam
  call genfilname(setfilext=.true.)
end subroutine xsgeneigvec
