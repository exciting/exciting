
! Copyright (C) 2005-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmatxs
! !INTERFACE:
subroutine writepmatasc
! !USES:
  use modinput, only: input
  use modmpi, only: procs, rank, firstofset, lastofset, barrier
  use mod_misc, only: task, filext
  use mod_kpoint, only: nkpt, vkl
  use mod_Gkvector, only: ngkmax, vgkl, ngk, gkc, tpgkc, sfacgk, igkig, vgkc
  use mod_APW_LO, only: apwordmax, nlotot, nlomax, lolmax
  use mod_muffin_tin, only: lmmaxapw
  use mod_atoms, only: natmtot
  use mod_eigensystem, only: nmatmax 
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use modxas, only: ncg
  use modxs, only: tscreen, fnpmat, fnpmat_t, kpari,&
                  & kparf, hybridhf, ripaa, ripalo,&
                  & riploa, riplolo, apwcmt, locmt,&
                  & unitout, iqmtgamma, fhdf5
  use m_putpmat
  use m_genfilname
  use mod_hdf5
  use m_getunit, only: getunit

! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}, {\tt PMAT\_XS.OUT} or
!   {\tt PMAT\_SCR.OUT} depending on the context of the execution.
!
! !REVISION HISTORY:
!   Created 2006 (S. Sagmeister)
!   Modifications, August 2010 (S. Sagmeister)
!   Re-purposed for HDF5 output in BSE module, July 2017 (C. Vorwerk)
!EOP
!BOC

  implicit none

  ! Local variables
  integer :: ik, reclen
  character(32) :: fnam
  complex(8), allocatable :: apwalmt(:, :, :, :)
  complex(8), allocatable :: evecfvt(:, :)
  complex(8), allocatable :: evecsvt(:, :)
  complex(8), allocatable :: pmat(:, :, :,:)
  character(256) :: string
  logical :: fast
  character(*), parameter :: thisname="writepmatxs"
  character(128) :: gname, cik
  ! External functions
  logical, external :: tqgamma

  integer :: ist, ist1, ist2, oct, un

  !write(*,*) "writepmatxs here at rank", rank

  ! Initialise universal variables
  call init0
  ! k-point setup
  ! Also allocated the radial functions (mod_APW_LO)
  call init1
  if(task .ne. 120) then 
    ! q-point and qmt-point setup
    !   Init 2 sets up (task 320/420):
    !   * A list of momentum transfer vectors form the q-point list 
    !     (modxs::vqmtl and mod_qpoint::vql)
    !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
    !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
    !   * G+qmt quantities (modxs)
    !   * The square root of the Coulomb potential for the G+qmt points
    !   * Reads STATE.OUT
    !   * Generates radial functions (mod_APW_LO)
    call init2
  end if

  ! Check if fast (default) version of matrix elements is used
  fast=.false.

  ! Check if first Q-point in list is the gamma point,
  ! if not return.
  if(.not. tqgamma(1)) then
    if(rank == 0) then 
      write(unitout,*) "Warning(writepmatxs): First Q pont not gamma, returning"
    end if
    return
  end if

  ! Generate index ranges for parallel execution
  call genparidxran('k', nkpt)

  allocate(apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
  allocate(evecfvt(nmatmax, nstfv))
  allocate(evecsvt(nstsv, nstsv))

  ! Allocate the momentum matrix elements array
  if (input%xs%bse%xas) then ! Allocation for xas calculation
    allocate(pmat(3, ncg, nstsv,nkpt))
  else
    allocate(pmat(3, nstsv, nstsv,nkpt))
  end if  

  ! Get eigenvectors for qmt=0, i.e. set file extension to _QMT001
  call genfilname(iqmt=iqmtgamma, setfilext=.true.)

    if(allocated(apwcmt)) deallocate(apwcmt)
    allocate(apwcmt(nstfv, apwordmax, lmmaxapw, natmtot))
    if(allocated(ripaa)) deallocate(ripaa)
    allocate(ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw,natmtot, 3))
    if(nlotot .gt. 0) then
      if(allocated(locmt)) deallocate(locmt)
      allocate(locmt(nstfv, nlomax,-lolmax:lolmax, natmtot))
      if(allocated(ripalo)) deallocate(ripalo)
      allocate(ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, natmtot, 3))
      if(allocated(riploa)) deallocate(riploa)
      allocate(riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, natmtot, 3))
      if(allocated(riplolo)) deallocate(riplolo)
      allocate(riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, natmtot, 3))
    end if

    ! Calculate gradient of radial functions times spherical harmonics
    call pmatrad
#ifdef _HDF5_
  if (rank == 0) then
    if (.not. hdf5_exist_group(fhdf5, "/", "pmat")) then
      call hdf5_create_group(fhdf5,"/", "pmat")
    end if
  end if
#else
  call getunit(un)
  open(un, File="PMAT_XS_ASC.OUT",Action="write")
#endif
  kloop: do ik = kpari, kparf
    if(task .ne. 120) call chkpt(2, (/ task, ik /), 'ematqk:&
     & task, k - point index; momentum matrix elements')

    ! Get the eigenvectors and values from file
    call getevecfv(vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
    call getevecsv(vkl(1, ik), evecsvt)

    ! Find the matching coefficients
    call match(ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik),&
      & sfacgk(1, 1, 1, ik), apwalmt)

    ! Generate apw expansion coefficients for muffin-tin
     call genapwcmt(input%groundstate%lmaxapw, ngk(1, ik), 1,&
       & nstfv, apwalmt, evecfvt, apwcmt)

    ! Generate local orbital expansion coefficients for muffin-tin
    if(nlotot .gt. 0) & 
      & call genlocmt(ngk(1, ik), 1, nstfv, evecfvt, locmt)

    ! Calculate the momentum matrix elements
    if((input%xs%bse%xas) .and. (task .le. 400)) then
      call genpmatcorxs(ik, ngk(1, ik), apwalmt, evecfvt, evecsvt, pmat(:,:,:,ik))
    else
      call genpmatxs(ngk(1, ik), igkig(1, 1, ik),&
       & vgkc(1, 1, 1, ik), evecfvt, evecsvt, pmat(:,:,:,ik))
    end if
#ifndef _HDF5_
    if (.NOT. input%xs%bse%xas) then
      Do ist1 = 1, nstsv
        Do ist2 = 1, nstsv
          Do oct = 1, 3
            Write (un, '(3i8, i4, 2g18.10)') ik, ist1, ist2, oct, &
              & pmat (oct, ist1, ist2, ik)
          End Do
        End Do
      End Do
    else
      Do ist1 = 1, ncg
        Do ist2 = 1, nstsv
          Do oct = 1, 3
            Write (un, '(3i8, i4, 2g18.10)') ik, ist1, ist2, oct, &
              & pmat (oct, ist1, ist2, ik)
          End Do
        End Do
      End Do
    end if
#endif
  end do kloop
#ifdef _HDF5_
  if (rank == 0) then
    do ik=1,nkpt 
      write(cik, '(I4.4)') ik
      if (.not. hdf5_exist_group(fhdf5, "/pmat/", trim(adjustl(cik)))) then
        call hdf5_create_group(fhdf5,"/pmat/", trim(adjustl(cik)))
      end if
      gname="/pmat/"//trim(adjustl(cik))
      ! Write hdf5
      !write(*,*) 'shape(pmat)=', shape(pmat)
      !write(*,*) "pmat(1,1,1)=", pmat(1,1,1)
      call hdf5_write(fhdf5,gname,'pmat', pmat(1,1,1,1), shape(pmat(:,:,:,ik)))
    end do
  end if
#endif
  call barrier(callername=trim(thisname))
#ifndef _HDF5_
  close(un)
#endif
  inquire(iolength=reclen) vkl(:, ik), nstsv, pmat
  deallocate(apwalmt, evecfvt, evecsvt, pmat)

  deallocate(apwcmt)
  deallocate(ripaa)
  if(nlotot .gt. 0) then
    deallocate(locmt)
    deallocate(ripalo, riploa, riplolo)
  end if

  call barrier(callername=trim(thisname))

  if(.not. input%sharedfs) call cpfiletonodes(trim(fnpmat))

  if(task .eq. 120) then
    close(50)
    if(rank==0) then
      write(*,*)
      write(*, '("Info(writepmatxs):")')
      write(*, '(" Momentum matrix elements written to file PMAT.OUT")')
      write(*,*)
    end if
  else
    write(unitout, '(a)') "Info(writepmatxs): Momentum matrix elements finished"
  end if

  ! Reset global file extension to default
  call genfilname(setfilext=.true.)

end subroutine writepmatasc
!EOC

