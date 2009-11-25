!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genvmat (vmt, vir, vmat)
! generates potential matrix elements for all states and k-points
      Use modinput
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: vmt (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (In) :: vir (ngrtot)
      Complex (8), Intent (Out) :: vmat (nstsv, nstsv, nkpt)
! local variables
      Integer :: is, ia, ias, irc, ir
      Integer :: ik
! local arrays
      Real (8), Allocatable :: rfmt (:, :, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: wfmt (:, :, :, :, :)
      Complex (8), Allocatable :: wfir (:, :, :)
! allocate local arrays
      Allocate (rfmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (wfmt(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir(ngrtot, nspinor, nstsv))
! convert muffin-tin potential to spherical coordinates
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, vmt(:, ir, ias), 1, 0.d0, rfmt(:, irc, ias), &
              & 1)
            End Do
         End Do
      End Do
! loop over k-points
      Do ik = 1, nkpt
! get the eigenvectors and values from file
         Call getevalsv (vkl(:, ik), evalsv)
         Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(:, ik), evecsv)
! find the matching coefficients
         Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
        & sfacgk(:, :, 1, ik), apwalm)
! calculate the wavefunctions for all states
         Call genwfsv (.False., ngk(1, ik), igkig(:, 1, ik), evalsv, &
        & apwalm, evecfv, evecsv, wfmt, wfir)
         Call genvmatk (rfmt, vir, wfmt, wfir, vmat(:, :, ik))
      End Do
      Deallocate (apwalm, evecfv, evecsv, wfmt, wfir)
      Return
End Subroutine
