! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine xsgeneigvec
  use mod_qpoint, only: nqpt, vql
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi
  use modxs, only: tscreen, skipgnd, qvkloff, unitout
  use m_writegqpts, only: writegqpts
  use m_filedel, only: filedel
  use m_genfilname, only: genfilname
  implicit none

  ! Local variables
  character(*), parameter :: thisnam = 'xsgeneigvec'
  real(8) :: vqlt(3)
  integer(4) :: iq, qi, qf

  ! External functions
  logical, external :: tqgamma

  ! Initialize universal variables
  call init0
  call init1
  call init2

  ! SCF already parallelized for k-point set
  ! Add gamma q-point
  qi = 0
  qf = nqpt

  ! If first q-point is gamma-point we copy files
  if(tqgamma(1)) qi = 1

  ! For screening there is just the gamma q-point
  if(tscreen) then
    qi = 0
    qf = 0
  end if

  ! Write q-points
  if(rank .eq. 0) call writeqpts

  ! Calculate eigenvectors for each q-point (k+q point set)
  qloop: do iq = qi, qf

    if(.not. tscreen) call genfilname(iqmt=max(0, iq), setfilext=.true.)

    if(skipgnd) filext = '.OUT'

    vqlt(:) = 0.d0
    if(iq .ne. 0) then
      vqlt(:) = vql(:, iq)
      if(rank .eq. 0) call writegqpts(iq, filext)
    end if

    ! Write eigenvectors, -values, occupancies and contracted mt coefficients
    call writeevec(vqlt, qvkloff(1, iq), filext)

    if(.not. tscreen) then
      write(unitout, &
        & '("Info(", a, "): eigenvectors generated for Q-point (iq, vql below)")')&
        & thisnam
      write(unitout, '(i6, 3g18.10)') iq, vqlt(:)
    else
      write(unitout, '(a)') 'Info(' // thisnam // '):&
        & eigenvectors generated for&
        & associated(input%xs%screening)/screened interaction:'
    end if

    if(rank .Eq. 0) Then
      ! Safely remove unnecessary files
      call filedel('EQATOMS'//trim(filext))
      call filedel('EVALCORE'//trim(filext))
      call filedel('FERMIDOS'//trim(filext))
      call filedel('GEOMETRY'//trim(filext))
      call filedel('LATTICE'//trim(filext))
      call filedel('IADIST'//trim(filext))
      call filedel('LINENGY'//trim(filext))
      call filedel('SYMCRYS'//trim(filext))
      call filedel('SYMLAT'//trim(filext))
      call filedel('SYMSITE'//trim(filext))
      call filedel('TOTENERGY'//trim(filext))
      call filedel('EVALFV'//trim(filext))
      call filedel('RMSDVEFF'//trim(filext))
      call filedel('DTOTENERGY'//trim(filext))
      if(input%groundstate%tforce) call filedel('DFORCEMAX'//trim(filext))
      call filedel('CHGDIST'//trim(filext))
      call filedel('SYMGENR'//trim(filext))
      call filedel('SYMINV'//trim(filext))
      call filedel('SYMMULT'//trim(filext))
      call filedel('SYMMULT_TABLE'//trim(filext))
      call filedel('SYMT2'//trim(filext))
    end if

  ! End loop over q-points
  end do qloop

  if((rank .eq. 0 .or. (.not. input%sharedfs .and. firstinnode))&
   & .and. tqgamma(1) .and. ( .not. tscreen)) then

    write(unitout, '("Info(", a, "): first Q-point is Gamma-point&
     & - copying relevant files")') thisnam

    ! Write files again one by one
    call copyfilesq0

  end if

  call barrier

  write(unitout, '("Info(", a, "): generation of eigenvectors finished")') thisnam

  call genfilname(setfilext=.true.)

end subroutine xsgeneigvec
