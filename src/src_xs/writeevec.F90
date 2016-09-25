! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine writeevec(vq, voff, filxt)
  use mod_constants, only: zzero
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_spin, only: nspnfv 
  use mod_Gkvector, only: ngkmax, vgkl, ngk, gkc, tpgkc, sfacgk
  use mod_APW_LO, only: apwordmax, nlomax, lolmax
  use mod_muffin_tin, only: lmmaxapw
  use mod_atoms, only: natmtot
  use mod_kpoint, only: nkpt, vkl
  use modinput, only: input
  use modmpi, only: rank, procs, barrier, endloopbarrier
  use modxs, only: isreadstate0, evecfv, apwcmt, locmt,&
                 & kpari, kparf 
  use m_gndstateq, only: gndstateq
  use m_filedel, only: filedel

  implicit none

  ! Arguments
  real(8), intent(in) :: vq(3), voff(3)
  character(*), intent(in) :: filxt
  ! Local variables
  integer :: ik, j
  complex(8), allocatable :: apwalm(:, :, :, :)

  ! Read from STATE.OUT exclusively
  isreadstate0 = .true.

  ! SCF calculation with one cycle
  call gndstateq(voff, filxt)
  
  ! Allocate first variational eigenvectors 
  if(allocated(evecfv)) deallocate(evecfv)
  allocate(evecfv(nmatmax, nstfv, nspnfv))

  ! Local array
  allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))

  ! Allocate space for expansion coefficients of the APWs
  if(allocated(apwcmt)) deallocate(apwcmt)
  allocate(apwcmt(nstfv, apwordmax, lmmaxapw, natmtot))

  ! Allocate expansion coefficients of the LOs
  if(allocated(locmt)) deallocate(locmt)
  allocate(locmt(nstfv, nlomax, -lolmax:lolmax, natmtot))

  ! Delete existing coefficients files
  if(rank .eq. 0) call filedel('APWCMT'//trim(filxt))
  if(rank .eq. 0) call filedel('LOCMT'//trim(filxt))

  ! Get the k-point indices for this rank
  call genparidxran('k', nkpt)

  do ik = kpari, kparf

    apwcmt(:, :, :, :) = zzero
    locmt(:, :, :, :) = zzero

    ! Read the first variational eigenvectors for k-point ik
    ! form file.
    call getevecfv(vkl(1, ik), vgkl(1, 1, 1, ik), evecfv)

    ! Compute the matching coefficients of the APWs
    call match(ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik),&
      & sfacgk(1, 1, 1, ik), apwalm)

    ! Extract the expansion coefficients corresponding to the
    ! APWs form evecfv.
    call genapwcmt(input%groundstate%lmaxapw, ngk(1, ik), 1,&
      & nstfv, apwalm, evecfv, apwcmt)

    ! Extract the expansion coefficients corresponding to the
    ! LOs form evecfv.
    call genlocmt(ngk(1, ik), 1, nstfv, evecfv, locmt)

    do j = 0, procs - 1
      if(rank .eq. j) then
         call putapwcmt('APWCMT'//trim(filxt), ik, vkl(1, ik), vq, apwcmt)
         call putlocmt('LOCMT'//trim(filxt), ik, vkl(1, ik), vq, locmt)
      end if
      call barrier
    end do

  end do

  call endloopbarrier(nkpt, procs)

  isreadstate0 = .false.

  deallocate(evecfv, apwalm, apwcmt, locmt)

end subroutine writeevec
