!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! Modified March 2014 (UW)
Subroutine exxengy
      Use modmain
      Use modinput
      Use modmpi
      Implicit None
! local variables
      Integer :: is, ia, nrc, m1, m2
      Integer :: ik, ist, jst
      Real (8) :: evv, ecv, ecc
      Complex (8) zt1
! allocatable arrays
      Complex (8), Allocatable :: wfcr1 (:, :, :)
      Complex (8), Allocatable :: wfcr2 (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :)
      Complex (8), Allocatable :: zvclmt (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
! external functions
      Complex (8) zfmtinp
      External zfmtinp
      Allocate (wfcr1(lmmaxvr, nrcmtmax, 2))
      Allocate (wfcr2(lmmaxvr, nrcmtmax, 2))
      Allocate (zrhomt(lmmaxvr, nrcmtmax))
      Allocate (zvclmt(lmmaxvr, nrcmtmax))
      Allocate (zfmt(lmmaxvr, nrcmtmax))
      evv = 0.d0
      ecv = 0.d0
      ecc = 0.d0
#ifdef MPI
         Do ik = firstk (rank), lastk (rank)
            Write (*, '("Info(exxengy): ", I6, " of ", I6, " k-points on pr&
        &oc:", I6)') ik, nkpt, rank

#else
         Do ik = 1, nkpt
             Write (*, '("Info(exxengy): ", I6, " of ", I6, " k-points")') ik, nkpt
#endif       
             call exxengyk (evv, ecv, ik) 
         End Do
#ifdef MPI
       call MPI_ALLREDUCE(MPI_IN_PLACE, evv, 1,  MPI_DOUBLE_PRECISION,  MPI_SUM, MPI_COMM_WORLD, ierr)
       call MPI_ALLREDUCE(MPI_IN_PLACE, ecv, 1,  MPI_DOUBLE_PRECISION,  MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            Do jst = 1, spnst (is)
               If (spcore(jst, is)) Then
                  Do m2 = - spk (jst, is), spk (jst, is) - 1
                     Call wavefcr (input%groundstate%lradstep, is, ia, &
                    & jst, m2, nrcmtmax, wfcr2)
                     Do ist = 1, spnst (is)
                        If (spcore(ist, is)) Then
                           Do m1 = - spk (ist, is), spk (ist, is) - 1
                              Call wavefcr (input%groundstate%lradstep, &
                             & is, ia, ist, m1, nrcmtmax, wfcr1)
! calculate the complex overlap density
                              Call vnlrhomt (.True., is, wfcr1(:, :, &
                             & 1), wfcr2(:, :, 1), zrhomt)
                              Call vnlrhomt (.True., is, wfcr1(:, :, &
                             & 2), wfcr2(:, :, 2), zfmt)
                              zrhomt (:, 1:nrc) = zrhomt (:, 1:nrc) + &
                             & zfmt (:, 1:nrc)
! calculate the Coulomb potential
                              Call zpotclmt (input%groundstate%ptnucl, &
                             & input%groundstate%lmaxvr, nrc, rcmt(:, &
                             & is), 0.d0, lmmaxvr, zrhomt, zvclmt)
                              zt1 = zfmtinp (.True., &
                             & input%groundstate%lmaxvr, nrc, rcmt(:, &
                             & is), lmmaxvr, zrhomt, zvclmt)
                              ecc = ecc - 0.5d0 * dble (zt1)
                           End Do
! end loop over ist
                        End If
                     End Do
                  End Do
! end loop over jst
               End If
            End Do
! end loops over atoms and species
         End Do
      End Do
! total exchange energy
      engyx = engyx + ex_coef*(evv + ecv + ecc)
      Deallocate (wfcr1, wfcr2, zrhomt, zvclmt, zfmt)
      Return
End Subroutine
