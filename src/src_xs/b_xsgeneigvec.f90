! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine b_xsgeneigvec(qi, nqpts, vql, qvkloff, tscr, tmqmt, tminus)
  use modmpi
  use modinput, only: input
  use modxs, only: unitout
  use mod_misc, only: filext
  use m_filedel, only: filedel
  use m_genfilname, only: genfilname

  implicit none

  ! Arguments
  integer(4), intent(in) :: qi, nqpts
  real(8), intent(in) :: vql(3,nqpts), qvkloff(3,nqpts)
  logical, intent(in) :: tscr, tmqmt, tminus
  real(8), parameter :: epslat=1.d-6

  ! Local variables
  character(*), parameter :: thisnam = 'b_xsgeneigvec'
  integer(4) :: iq

  ! Check 
  if(qi > nqpts) then 
    write(*,*) "ERROR(b_xsgeneigvec): qi > nqpts"
    call terminate
  end if

  write(*,*) "b_xsgeneigvec here at rank", rank

  ! Ground state SCF already parallelized for k-point set
  ! Calculate eigenvectors for each q-point (k+q point set)
  do iq = qi, nqpts

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
    !                    _mqmt to indicate that k-qmt was used instead of k+qmt 
    !                    _m denotes the -(k+-qmt)-grid was used (differs for vkloffs /= 0)
    !                    
    ! 
    ! Set file extension
    if(.not. tscr) then
      if(.not. tmqmt) then 
        if(.not. tminus) then 
          call genfilname(iqmt=iq, setfilext=.true.)
        else
          call genfilname(iqmt=iq, auxtype="m", setfilext=.true.)
        end if
      else
        if(.not. tminus) then
          call genfilname(iqmt=iq, auxtype="mqmt", setfilext=.true.)
        else
          call genfilname(iqmt=iq, auxtype="mqmt_m", setfilext=.true.)
        end if
      end if
    else 
      if(.not. tmqmt) then 
        if(.not. tminus) then 
          call genfilname(iqmt=iq, scrtype='', setfilext=.true.)
        else
          call genfilname(iqmt=iq, scrtype='', auxtype="m", setfilext=.true.)
        end if
      else
        if(.not. tminus) then
          call genfilname(iqmt=iq, scrtype='', auxtype="mqmt", setfilext=.true.)
        else
          call genfilname(iqmt=iq, scrtype='', auxtype="mqmt_m", setfilext=.true.)
        end if
      end if
    end if

    !if(iq /= 1 .and. all(qvkloff(1:3,iq) == qvkloff(1:3,1))) then
    !  if(rank == 0) then 
    !    write(unitout, '("Info(", a, "): Q-point ", i4,&
    !    &  " has the same k-grid offset as the Q=0 qmt point, skipping GS calculation.")')&
    !    &  thisnam, iq
    !  end if
    !  cycle
    !end if

    ! Call groundstate and write results
    call writeevec(vql(1:3,iq), qvkloff(1:3, iq), filext)

    if(rank == 0) then 
      write(unitout, &
        & '("Info(", a, "): Eigenvectors generated for Q-point ", i6)')&
        & thisnam, iq
      if(tscr) then 
        write(unitout, '("Using screening GS parameters.")')
      end if
      write(unitout, '("vql = ", 3g18.10)') vql(1:3,iq)
      if(tmqmt) then 
        if(tminus) then 
          write(unitout, '("-(k-qmt)-grid offset derived form qmt and k offset:")')
        else
          write(unitout, '("(k-qmt)-grid offset derived form qmt and k offset:")')
        end if
      else
        if(tminus) then 
          write(unitout, '("-(k+qmt)-grid offset derived form qmt and k offset:")')
        else
          write(unitout, '("(k+qmt)-grid offset derived form qmt and k offset:")')
        end if
      end if
      write(unitout, '("qvkloff = ", 3g18.10)') qvkloff(1:3,iq)
      write(unitout, '("qvkloff/ngridk = ", 3g18.10)') qvkloff(1:3,iq)/input%xs%ngridk
      write(unitout, *)
    end if

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

  call barrier

  if(rank == 0) then 
    write(unitout, '("Info(", a, "): Generation of eigenvectors finished")') thisnam
    write(unitout, *)
  end if

end subroutine b_xsgeneigvec
