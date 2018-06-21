! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: xsgeneigvec
! !INTERFACE:
subroutine xsgeneigvec(qi, qf, nqpts, vql, qvkloff, tscr, tmqmt)
! !USES:
  use modmpi
  use modinput, only: input
  use modxs, only: unitout, vqlmt, vqcmt
  use mod_misc, only: filext
  use m_filedel, only: filedel
  use m_genfilname, only: genfilname
! !INPUT/OUTPUT PARAMETERS:
! In:
!   integer(4) :: qi, qf ! Range qi:qf of then nqpts momentum transfer Qs form Q-point list
!   integer(4) :: nqpts  ! for which to do one-shot GS runs 
!   real(8)    :: vql(3,nqpts)     ! Q-point vectors form list
!   real(8)    :: qvkloff(3,nqpts) ! Offsets of corresponding k-grids
!   logical    :: tscr  ! If true use screening file extension
!   logical    :: tmqmt ! If true use file extension indicating that -Q/2 was used
! 
! !DESCRIPTION:
!   The routine generates the one-shot GS quantities used in the BSE.
!   The routine writes the files {\tt APWCMT\_*.OUT, EIGVAL\_*.OUT, EVALSV\_*.OUT,
!   EVECFV\_*.OUT, EVECSV\_*.OUT, BONDLENGTH\_*.OUT, EFERMI\_*.OUT,
!   LOCMT\_*.OUT, OCCSV\_*.OUT, geometry\_*.xml}. The file extension contains
!   the number of the considered Q-point {\tt \_QMTxyz}, the information whether 
!   the +Q/2 or -Q/2 shift was used {\tt \_m}. If the screening parameter were
!   used the extension {\tt \_SCR} is added.
!
! !REVISION HISTORY:
!   Based on xsgeneigvec
!   Created. 2016 (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer(4), intent(in) :: qi, qf, nqpts
  real(8), intent(in) :: vql(3,nqpts), qvkloff(3,nqpts)
  logical, intent(in) :: tscr, tmqmt
  real(8), parameter :: epslat=1.d-6

  ! Local variables
  character(*), parameter :: thisname = 'xsgeneigvec'
  integer(4) :: iq

  ! Check 
  if(qi > nqpts .or. qi < 1) then 
    write(*,*) "ERROR(xsgeneigvec): qi invalid"
    call terminate
  end if
  if(qf > nqpts .or. qf < 1 .or. qf < qi) then 
    write(*,*) "ERROR(xsgeneigvec): qf invalid"
    call terminate
  end if
  if(tscr .and. (qi /= 1 .or. qf /= 1)) then
    write(*,'("Error(",a,"):&
      & Screening eigenvectors only needed on (QMT=1) original k-grid.")')&
      & thisname
    call terminate()
  end if

  ! Ground state SCF already parallelized for k-point set
  ! Calculate eigenvectors for each qmt-point using the passed offsets qvkloff
  do iq = qi, qf

    ! Execute a GS run with an offset derived from qmt and 
    ! write eigenvectors, -values and occupancies.
    !   Writes: EVECSV_*, EVECFV_*, EVALSV_*, OCCSV_*
    !
    ! For later calculation of plane wave matrix elements:
    ! The product of APW matching coefficients with APW eigenvector coefficients
    ! are computed and stored. Also the eigenvector coefficients pertaining to the 
    ! local orbitals are written.
    !   Writes: APWCMT_*, LOCMT_*
    !
    ! The ending scheme: _SCR to indicate that the screening GS parameter were used
    !                    _QMTXYZ to indicate which qmt list point was used
    !                    _m to indicate that k-qmt/2 was used instead of k+qmt/2 
    !                    
    ! 
    ! Set file extension
    if(.not. tscr) then
      if(.not. tmqmt) then 
        call genfilname(iqmt=iq, setfilext=.true.)
      else
        call genfilname(iqmt=iq, auxtype="m", setfilext=.true.)
      end if
    else 
      if(.not. tmqmt) then 
        call genfilname(scrtype='', setfilext=.true.)
      else
        write(*,'("Error(",a,"):&
          & Screening eigenvectors only needed on original k-grid.")')&
          & thisname
        call terminate()
      end if
    end if

    write(unitout, &
      & '("Info(", a, "): Generating eigenvectors for Q-point ", i6)')&
      & thisname, iq
    if(tscr) then 
      write(unitout, '("Using screening GS parameters.")')
    end if
    if (input%xs%BSE%outputlevelnumber == 1) then
      write(unitout, '("vqlmt = ", 3g18.10)') vqlmt(1:3,iq)
      write(unitout, '("vqcmt = ", 3g18.10)') vqcmt(1:3,iq)
      write(unitout, '("Norm2 = ", 1g18.10)') norm2(vqcmt(1:3,iq))
      if(norm2(vqcmt(1:3,iq)) > 1.0d-8) then 
        write(unitout, '("vqcmt/Norm2 = ", 3g18.10)') vqcmt(1:3,iq)/norm2(vqcmt(1:3,iq))
      end if
      write(unitout, '("k-grid offset derived form qmt and k offset:")')
      write(unitout, '("qvkloff = ", 3g18.10)') qvkloff(1:3,iq)
      write(unitout, '("qvkloff/ngridk = ", 3g18.10)') qvkloff(1:3,iq)/input%xs%ngridk
      if(iq /= 1 .and. all(abs(qvkloff(1:3,iq)-qvkloff(1:3,1)) < epslat)) then
        write(unitout, '("Info(", a, "): Q-point ", i4,&
        &  " has the same k-grid offset as the Q=0 qmt point, skipping GS calculation.")')&
        &  thisname, iq
    end if
    call printline(unitout, "-")
      cycle
    else
      call printline(unitout, "-")
    end if

    ! Call groundstate and write results
    call writeevec(vql(1:3,iq), qvkloff(1:3, iq), filext)

    if(rank == 0) then
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
  end do

  call barrier(callername=trim(thisname))

  if(rank == 0) then 
    write(unitout, '("Info(", a, "): Generation of eigenvectors finished")') thisname
    write(unitout, *)
  end if

end subroutine xsgeneigvec
!EOC
