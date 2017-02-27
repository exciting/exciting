! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematqalloc
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use mod_spin, only: nspnfv
  use mod_kpoint, only: nkpt
  use mod_APW_LO, only: nlotot, apwordmax, nlomax, lolmax
  use mod_atoms, only: natmtot
  use modxs, only: evecfv, evecfvb, evecfv0, evecfv0b, evalsv0, nmatmax0,&
                 & xih, apwcmt, apwcmtb, lmmaxapwwf, apwcmt0, apwcmt0b,&
                 & locmt, locmt0

  implicit none

  ! Allocate eigenvalue and eigenvector arrays
  if(allocated(evecfv)) deallocate(evecfv)
  allocate(evecfv(nmatmax, nstfv, nspnfv))
  if(allocated(evecfvb)) deallocate(evecfvb)
  allocate(evecfvb(nmatmax, nstfv, nspnfv))

  if(allocated(evecfv0)) deallocate(evecfv0)
  allocate(evecfv0(nmatmax0, nstfv, nspnfv))
  if(allocated(evecfv0b)) deallocate(evecfv0b)
  allocate(evecfv0b(nmatmax0, nstfv, nspnfv))

  if(allocated(evalsv0)) deallocate(evalsv0)
  allocate(evalsv0(nstsv, nkpt))
  
  ! Allocate helper matrix
  if(allocated(xih)) deallocate(xih)
  allocate(xih(nlotot, nlotot))

  ! Allocate contracted coefficients array for apw-part
  if(allocated(apwcmt)) deallocate(apwcmt)
  allocate(apwcmt(nstfv, apwordmax, lmmaxapwwf, natmtot))

  if(allocated(apwcmt0)) deallocate(apwcmt0)
  allocate(apwcmt0(nstfv, apwordmax, lmmaxapwwf, natmtot))

  ! Allocate contracted coefficients array for local orbitals-part
  if(allocated(locmt)) deallocate(locmt)
  allocate(locmt(nstfv, nlomax,-lolmax:lolmax, natmtot))
  if(allocated(locmt0)) deallocate(locmt0)
  allocate(locmt0(nstfv, nlomax,-lolmax:lolmax, natmtot))

end subroutine ematqalloc
