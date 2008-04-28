
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqalloc
  use modmain
  use modmpi
  use modxs
  implicit none
  ! allocate eigenvalue and eigenvector arrays
  if (allocated(evecfv)) deallocate(evecfv)
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  if (allocated(evecfv0)) deallocate(evecfv0)
  allocate(evecfv0(nmatmax0,nstfv,nspnfv))
  if (allocated(evalsv0)) deallocate(evalsv0)
  allocate(evalsv0(nstsv,nkpt))
  ! allocate helper matrix
  if (allocated(xih)) deallocate(xih)
  allocate(xih(nlotot,nlotot))
  ! allocate contracted coefficients array for APW-part
  if (allocated(apwdlm)) deallocate(apwdlm)
  allocate(apwdlm(nstsv,apwordmax,lmmaxapwtd,natmtot))
  if (allocated(apwdlm0)) deallocate(apwdlm0)
  allocate(apwdlm0(nstsv,apwordmax,lmmaxapwtd,natmtot))
  ! allocate contracted coefficients array for local orbitals-part
  if (allocated(lodlm)) deallocate(lodlm)
  allocate(lodlm(nstsv,nlomax,-lolmax:lolmax,natmtot))
  if (allocated(lodlm0)) deallocate(lodlm0)
  allocate(lodlm0(nstsv,nlomax,-lolmax:lolmax,natmtot))
end subroutine ematqalloc
