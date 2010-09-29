!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: xssave0
! !INTERFACE:
!
!
Subroutine xssave0
! !USES:
      Use modmain
      Use modxs
! !DESCRIPTION:
!   This routine should be called after init0, init1 and init2 in order to save
!   variables related to the k-point set for ${\bf q}=0$.
!
! !REVISION HISTORY:
!   Created March 2005 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! allocate the k-point arrays
      If (allocated(vkl0)) deallocate (vkl0)
      Allocate (vkl0(3, nkptnr))
  ! allocate the G+k-point arrays
      If (allocated(ngk0)) deallocate (ngk0)
      If (allocated(igkig0)) deallocate (igkig0)
      If (allocated(vgkl0)) deallocate (vgkl0)
      If (allocated(vgkc0)) deallocate (vgkc0)
      If (allocated(gkc0)) deallocate (gkc0)
      If (allocated(tpgkc0)) deallocate (tpgkc0)
      If (allocated(sfacgk0)) deallocate (sfacgk0)
      Allocate (ngk0(nspnfv, nkpt))
      Allocate (igkig0(ngkmax, nspnfv, nkpt))
      Allocate (vgkl0(3, ngkmax, nspnfv, nkpt))
      Allocate (vgkc0(3, ngkmax, nspnfv, nkpt))
      Allocate (gkc0(ngkmax, nspnfv, nkpt))
      Allocate (tpgkc0(2, ngkmax, nspnfv, nkpt))
      Allocate (sfacgk0(ngkmax, natmtot, nspnfv, nkpt))
  ! overlap and Hamiltonian matrix sizes
      If (allocated(nmat0)) deallocate (nmat0)
      Allocate (nmat0(nspnfv, nkpt))
  ! save variables for k-vectors
      nkpt0 = nkpt
      vkl0 (:, :) = vkl (:, :)
  ! save variables for G+k-vectors
      ngkmax0 = ngkmax
      ngk0 (:, :) = ngk (:, :)
      igkig0 (:, :, :) = igkig (:, :, :)
      vgkl0 (:, :, :, :) = vgkl (:, :, :, :)
      vgkc0 (:, :, :, :) = vgkc (:, :, :, :)
      gkc0 (:, :, :) = gkc (:, :, :)
      tpgkc0 (:, :, :, :) = tpgkc (:, :, :, :)
      sfacgk0 (:, :, :, :) = sfacgk (:, :, :, :)
  ! save variables for overlap and Hamiltonian matrix sizes
      nmatmax0 = nmatmax
      nmat0 (:, :) = nmat (:, :)
End Subroutine xssave0
!EOC
