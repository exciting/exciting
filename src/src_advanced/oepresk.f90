!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepresk (ik, vnlcv, vnlvv, dvxmt, dvxir, dbxmt, dbxir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: vnlcv (ncrmax, natmtot, nstsv, nkpt)
      Complex (8), Intent (In) :: vnlvv (nstsv, nstsv, nkpt)
      Real (8), Intent (Inout) :: dvxmt (lmmaxvr, nrcmtmax, natmtot)
      Real (8), Intent (Inout) :: dvxir (ngrtot)
      Real (8), Intent (Inout) :: dbxmt (lmmaxvr, nrcmtmax, natmtot, &
     & ndmag)
      Real (8), Intent (Inout) :: dbxir (ngrtot, ndmag)
! local variables
      Integer :: is, ia, ias, ist, jst
      Integer :: nrc, ic, m, idm
      Real (8) :: de
      Complex (8) zde, zvnl, zrvx, zmbx, zt1, zt2
! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: wfmt (:, :, :, :, :)
      Complex (8), Allocatable :: wfir (:, :, :)
      Complex (8), Allocatable :: wfcr (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zmagmt (:, :, :, :)
      Complex (8), Allocatable :: zmagir (:, :)
      Complex (8), Allocatable :: zvfmt (:, :, :)
      Complex (8), Allocatable :: zfmt (:, :)
! external functions
      Complex (8) zfinp, zfmtinp
      External zfinp, zfmtinp
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (wfmt(lmmaxvr, nrcmtmax, natmtot, nspinor, nstsv))
      Allocate (wfir(ngrtot, nspinor, nstsv))
      Allocate (wfcr(lmmaxvr, nrcmtmax, 2))
      Allocate (zrhomt(lmmaxvr, nrcmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      If (associated(input%groundstate%spin)) Then
         Allocate (zmagmt(lmmaxvr, nrcmtmax, natmtot, ndmag))
         Allocate (zmagir(ngrtot, ndmag))
         Allocate (zvfmt(lmmaxvr, nrcmtmax, ndmag))
         Allocate (zfmt(lmmaxvr, nrcmtmax))
      End If
! get the eigenvalues/vectors from file for input k-point
      Call getevalsv (vkl(:, ik), evalsv(:, ik))
      Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
      Call getevecsv (vkl(:, ik), evecsv)
! find the matching coefficients
      Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
     & sfacgk(:, :, 1, ik), apwalm)
! calculate the wavefunctions for all states
      Call genwfsv (.False., ngk(1, ik), igkig(:, 1, ik), evalsv(:, &
     & ik), apwalm, evecfv, evecsv, wfmt, wfir)
!-----------------------------------------------------------!
!     core-conduction overlap density and magnetisation     !
!-----------------------------------------------------------!
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            ic = 0
            Do ist = 1, spnst (is)
               If (spcore(ist, is)) Then
                  Do m = - spk (ist, is), spk (ist, is) - 1
                     ic = ic + 1
! pass in m-1/2 to wavefcr
                     Call wavefcr (input%groundstate%lradstep, is, ia, &
                    & ist, m, nrcmtmax, wfcr)
                     Do jst = 1, nstsv
                        If (evalsv(jst, ik) .Gt. efermi) Then
                           de = evalcr (ist, ias) - evalsv (jst, ik)
                           zde = occmax * wkpt (ik) / &
                          & (de+zi*input%groundstate%swidth)
! calculate the complex overlap density in the muffin-tin
                           Call vnlrhomt (.False., is, wfcr(:, :, 1), &
                          & wfmt(:, :, ias, 1, jst), zrhomt(:, :, ias))
                           If (associated(input%groundstate%spin)) Then
                              Call vnlrhomt (.False., is, wfcr(:, :, &
                             & 2), wfmt(:, :, ias, 2, jst), zfmt)
                              zrhomt (:, 1:nrc, ias) = zrhomt (:, &
                             & 1:nrc, ias) + zfmt (:, 1:nrc)
                           End If
                           zvnl = conjg (vnlcv(ic, ias, jst, ik))
                           zrvx = zfmtinp (.False., &
                          & input%groundstate%lmaxvr, nrc, rcmt(:, is), &
                          & lmmaxvr, zrhomt(:, :, ias), zvxmt(:, :, &
                          & ias))
                           zt1 = zvnl - zrvx
! spin-polarised case
                           If (associated(input%groundstate%spin)) Then
                              Call oepmagmt (.False., is, wfcr(:, :, &
                             & 1), wfcr(:, :, 2), wfmt(:, :, ias, 1, &
                             & jst), wfmt(:, :, ias, 2, jst), zvfmt)
! integral of magnetisation dot exchange field
                              zmbx = 0.d0
                              Do idm = 1, ndmag
                                 zmbx = zmbx + zfmtinp (.False., &
                                & input%groundstate%lmaxvr, nrc, &
                                & rcmt(:, is), lmmaxvr, zvfmt(:, :, &
                                & idm), zbxmt(:, :, ias, idm))
                              End Do
                              zt1 = zt1 - zmbx
! end spin-polarised case
                           End If
                           zt2 = zde * zt1
! residuals for exchange potential and field
!$OMP CRITICAL
                           dvxmt (:, 1:nrc, ias) = dvxmt (:, 1:nrc, &
                          & ias) + dble (zt2*zrhomt(:, 1:nrc, ias))
                           Do idm = 1, ndmag
                              dbxmt (:, 1:nrc, ias, idm) = dbxmt (:, &
                             & 1:nrc, ias, idm) + dble (zt2*zvfmt(:, &
                             & 1:nrc, idm))
                           End Do
!$OMP END CRITICAL
! end loop over jst
                        End If
                     End Do
                  End Do
! end loop over ist
               End If
            End Do
! end loops over atoms and species
         End Do
      End Do
!--------------------------------------------------------------!
!     valence-conduction overlap density and magnetisation     !
!--------------------------------------------------------------!
      Do ist = 1, nstsv
         If (evalsv(ist, ik) .Lt. efermi) Then
            Do jst = 1, nstsv
               If (evalsv(jst, ik) .Gt. efermi) Then
! calculate the overlap density
                  Call vnlrho (.False., wfmt(:, :, :, :, ist), wfmt(:, &
                 & :, :, :, jst), wfir(:, :, ist), wfir(:, :, jst), &
                 & zrhomt, zrhoir)
                  de = evalsv (ist, ik) - evalsv (jst, ik)
                  zde = occmax * wkpt (ik) / &
                 & (de+zi*input%groundstate%swidth)
                  zvnl = conjg (vnlvv(ist, jst, ik))
                  zrvx = zfinp (.False., zrhomt, zvxmt, zrhoir, zvxir)
                  zt1 = zvnl - zrvx
! spin-polarised case
                  If (associated(input%groundstate%spin)) Then
                     Call oepmag (.False., wfmt(:, :, :, :, ist), &
                    & wfmt(:, :, :, :, jst), wfir(:, :, ist), wfir(:, &
                    & :, jst), zmagmt, zmagir)
! integral of magnetisation dot exchange field
                     zmbx = 0.d0
                     Do idm = 1, ndmag
                        zmbx = zmbx + zfinp (.False., zmagmt(:, :, :, &
                       & idm), zbxmt(:, :, :, idm), zmagir(:, idm), &
                       & zbxir(:, idm))
                     End Do
                     zt1 = zt1 - zmbx
                  End If
                  zt2 = zde * zt1
! residuals for exchange potential and field
!$OMP CRITICAL
                  Do is = 1, nspecies
                     nrc = nrcmt (is)
                     Do ia = 1, natoms (is)
                        ias = idxas (ia, is)
                        dvxmt (:, 1:nrc, ias) = dvxmt (:, 1:nrc, ias) + &
                       & dble (zt2*zrhomt(:, 1:nrc, ias))
                        Do idm = 1, ndmag
                           dbxmt (:, 1:nrc, ias, idm) = dbxmt (:, &
                          & 1:nrc, ias, idm) + dble (zt2*zmagmt(:, &
                          & 1:nrc, ias, idm))
                        End Do
                     End Do
                  End Do
                  dvxir (:) = dvxir (:) + dble (zt2*zrhoir(:))
                  Do idm = 1, ndmag
                     dbxir (:, idm) = dbxir (:, idm) + dble &
                    & (zt2*zmagir(:, idm))
                  End Do
!$OMP END CRITICAL
! end loop over jst
               End If
            End Do
! end loop over ist
         End If
      End Do
      Deallocate (apwalm, evecfv, evecsv)
      Deallocate (wfmt, wfir, wfcr, zrhomt, zrhoir)
      If (associated(input%groundstate%spin)) Then
         Deallocate (zmagmt, zmagir, zvfmt, zfmt)
      End If
      Return
End Subroutine
