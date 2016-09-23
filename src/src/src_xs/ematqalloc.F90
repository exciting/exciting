!
!
!
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqalloc
      Use modmain
      Use modmpi
      Use modxs
      Implicit None
  ! allocate eigenvalue and eigenvector arrays
      If (allocated(evecfv)) deallocate (evecfv)
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      If (allocated(evecfv0)) deallocate (evecfv0)
      Allocate (evecfv0(nmatmax0, nstfv, nspnfv))
      If (allocated(evalsv0)) deallocate (evalsv0)
      Allocate (evalsv0(nstsv, nkpt))
  ! allocate helper matrix
      If (allocated(xih)) deallocate (xih)
      Allocate (xih(nlotot, nlotot))
  ! allocate contracted coefficients array for APW-part
      If (allocated(apwcmt)) deallocate (apwcmt)
      Allocate (apwcmt(nstfv, apwordmax, lmmaxapwwf, natmtot))
      If (allocated(apwcmt0)) deallocate (apwcmt0)
      Allocate (apwcmt0(nstfv, apwordmax, lmmaxapwwf, natmtot))
  ! allocate contracted coefficients array for local orbitals-part
      If (allocated(locmt)) deallocate (locmt)
      Allocate (locmt(nstfv, nlomax,-lolmax:lolmax, natmtot))
      If (allocated(locmt0)) deallocate (locmt0)
      Allocate (locmt0(nstfv, nlomax,-lolmax:lolmax, natmtot))
End Subroutine ematqalloc
