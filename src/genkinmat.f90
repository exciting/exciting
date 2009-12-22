!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genkinmat
! generates kinetic matrix elements for all states and k-points
      Use modinput
      Use modmain
      Implicit None
! local variables
      Integer :: is, ia, ias, idm
      Integer :: ik, ist, ir, irc
! allocatable arrays
      Real (8), Allocatable :: rfmt (:, :, :)
      Real (8), Allocatable :: rvfmt (:, :, :, :)
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: wfmt (:, :, :, :, :)
      Complex (8), Allocatable :: wfir (:, :, :)
      Complex (8), Allocatable :: vmat (:, :)
      Complex (8), Allocatable :: bmat (:, :)
      Complex (8), Allocatable :: c (:, :)
! allocate local arrays
      Allocate (rfmt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (evalfv(nstfv, nspnfv))
      If (associated(input%groundstate%spin)) allocate (rvfmt(lmmaxvr, &
     & nrcmtmax, natmtot, ndmag))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (wfmt(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir(ngrtot, nspinor, nstsv))
      Allocate (vmat(nstsv, nstsv))
      Allocate (bmat(nstsv, nstsv))
      Allocate (c(nstsv, nstsv))
! convert muffin-tin effective potential and magnetic field to spherical
! coordinates
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, veffmt(:, ir, ias), 1, 0.d0, rfmt(:, irc, &
              & ias), 1)
               Do idm = 1, ndmag
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 & lmmaxvr, bxcmt(:, ir, ias, idm), 1, 0.d0, rvfmt(:, &
                 & irc, ias, idm), 1)
               End Do
            End Do
         End Do
      End Do
! loop over k-points
      Do ik = 1, nkpt
! solve the first- and second-variational secular equations
         Call seceqn (ik, evalfv, evecfv, evecsv)
! write the first variational eigenvalues/vectors to file (this ensures the
! phase in eigenvectors is the same for subsequent matrix element evaluations)
         Call putevalfv (ik, evalfv)
         Call putevecfv (ik, evecfv)
! find the matching coefficients
         Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
        & sfacgk(:, :, 1, ik), apwalm)
! calculate the wavefunctions for all states of the input k-point
         Call genwfsv (.False., ngk(1, ik), igkig(:, 1, ik), evalsv(:, &
        & ik), apwalm, evecfv, evecsv, wfmt, wfir)
! compute effective potential matrix elements
         Call genvmatk (rfmt, veffir, wfmt, wfir, kinmatc(:, :, ik))
         kinmatc (:, :, ik) = - kinmatc (:, :, ik)
! add second-variational eigenvalues along the diagonal
         Do ist = 1, nstsv
            kinmatc (ist, ist, ik) = kinmatc (ist, ist, ik) + evalsv &
           & (ist, ik)
         End Do
! compute the exchange-correlation magnetic field matrix elements
         If (associated(input%groundstate%spin)) Then
            Call genbmatk (rvfmt, bxcir, wfmt, wfir, bmat)
            kinmatc (:, :, ik) = kinmatc (:, :, ik) - bmat (:, :)
         End If
! rotate kinetic matrix elements to Cartesian basis
         Call zgemm ('N', 'C', nstsv, nstsv, nstsv, zone, kinmatc(:, :, &
        & ik), nstsv, evecsv, nstsv, zzero, c, nstsv)
         Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, evecsv, &
        & nstsv, c, nstsv, zzero, kinmatc(:, :, ik), nstsv)
      End Do
      If (associated(input%groundstate%spin)) deallocate (rvfmt)
      Deallocate (rfmt, evalfv, apwalm, evecfv, evecsv)
      Deallocate (wfmt, wfir, vmat, bmat, c)
      Return
End Subroutine
