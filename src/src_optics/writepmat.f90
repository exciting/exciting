!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writepmat
! !INTERFACE:
!
!
Subroutine writepmat
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, recl
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: pmat (:, :, :)
      Complex (8), Allocatable :: apwcmt (:, :, :, :)
      Complex (8), Allocatable :: locmt (:, :, :, :)
      Real (8), Allocatable :: ripaa (:, :, :, :, :, :)
      Real (8), Allocatable :: ripalo (:, :, :, :, :, :)
      Real (8), Allocatable :: riploa (:, :, :, :, :, :)
      Real (8), Allocatable :: riplolo (:, :, :, :, :, :)
! initialise universal variables
      Call init0
      Call init1
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
! allocate the momentum matrix elements array
      Allocate (pmat(3, nstsv, nstsv))
! read in the density and potentials from file
      Call readstate
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
      Allocate (ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw, &
     & natmtot, 3))
      Allocate (apwcmt(nstsv, apwordmax, lmmaxapw, natmtot))
      If (nlotot .Gt. 0) Then
         Allocate (ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, &
        & natmtot, 3))
         Allocate (riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, &
        & natmtot, 3))
         Allocate (riplolo(nlomax,-lolmax:lolmax, &
        & nlomax,-lolmax:lolmax, natmtot, 3))
         Allocate (locmt(nstsv, nlomax,-lolmax:lolmax, natmtot))
      End If
! calculate gradient of radial functions times spherical harmonics
      Call pmatrad (ripaa, ripalo, riploa, riplolo)
! find the record length
      Inquire (IoLength=Recl) pmat
      Open (50, File='PMAT.OUT', Action='WRITE', Form='UNFORMATTED', &
     & Access='DIRECT', Status='REPLACE', Recl=Recl)
      Do ik = 1, nkpt
! get the eigenvectors from file
         Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(:, ik), evecsv)
! find the matching coefficients
         Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
        & sfacgk(:, :, 1, ik), apwalm)
! generate APW expansion coefficients for muffin-tin
         Call genapwcmt (input%groundstate%lmaxapw, ngk(1, ik), 1, &
        & nstfv, apwalm, evecfv, apwcmt)
! generate local orbital expansion coefficients for muffin-tin
         If (nlotot .Gt. 0) Call genlocmt (ngk(1, ik), 1, nstfv, &
        & evecfv, locmt)
! calculate the momentum matrix elements
         Call genpmat2 (ngk(1, ik), igkig(:, 1, ik), vgkc(:, :, 1, ik), &
        & ripaa, ripalo, riploa, riplolo, apwcmt, locmt, evecfv, &
        & evecsv, pmat)
! calculate the momentum matrix elements
!!$  call genpmat(ngk(1,ik),igkig(:,1,ik),vgkc(:,:,1,ik),apwalm,evecfv,evecsv,pmat)
! write the matrix elements to direct-access file
         Write (50, Rec=ik) pmat
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writepmat):")')
      Write (*, '(" momentum matrix elements written to file PMAT.OUT")&
     &')
      Write (*,*)
      Deallocate (apwalm, evecfv, evecsv, pmat)
      Deallocate (ripaa, apwcmt)
      If (nlotot .Gt. 0) deallocate (ripalo, riploa, riplolo, locmt)
End Subroutine
!EOC
