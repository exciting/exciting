!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genwfsv
! !INTERFACE:
!
!
Subroutine genwfsv (tocc, ngp, igpig, evalsvp, apwalm, evecfv, evecsv, &
& wfmt, wfir)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tocc    : .true. if only occupied wavefunctions are required (in,logical)
!   ngp     : number of G+p-vectors (in,integer)
!   igpig   : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   evalsvp : second-variational eigenvalue for every state (in,real(nstsv))
!   apwalm  : APW matching coefficients
!             (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv  : first-variational eigenvectors (in,complex(nmatmax,nstfv))
!   evecsv  : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   wfmt    : muffin-tin part of the wavefunctions for every state in spherical
!             coordinates (out,complex(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
!   wfir    : interstitial part of the wavefunctions for every state
!             (out,complex(ngrtot,nspinor,nstsv))
! !DESCRIPTION:
!   Calculates the second-variational spinor wavefunctions in both the
!   muffin-tin and interstitial regions for every state of a particular
!   $k$-point. The wavefunctions in both regions are stored on a real-space
!   grid. A coarse radial mesh is assumed in the muffin-tins with with angular
!   momentum cut-off of {\tt lmaxvr}. If {\tt tocc} is {\tt .true.}, then only
!   the occupied states (those below the Fermi energy) are calculated.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tocc
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: evalsvp (nstsv)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Complex (8), Intent (Out) :: wfmt (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor, nstsv)
      Complex (8), Intent (Out) :: wfir (ngrtot, nspinor, nstsv)
! local variables
      Integer :: ispn, is, ia, ias
      Integer :: i, j, n, ist, igp, ifg
      Real (8) :: t1
      Complex (8) zt1
! allocatable arrays
      Logical, Allocatable :: done (:)
      Complex (8), Allocatable :: wfmt1 (:, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Allocate (done(nstfv))
      Allocate (wfmt1(lmmaxvr, nrcmtmax))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, nstfv))
!--------------------------------!
!     muffin-tin wavefunction    !
!--------------------------------!
      Do is = 1, nspecies
         n = lmmaxvr * nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            done (:) = .False.
            Do j = 1, nstsv
               If (( .Not. tocc) .Or. ((tocc) .And. (evalsvp(j) .Lt. &
              & efermi))) Then
                  If (input%groundstate%tevecsv) Then
! generate spinor wavefunction from second-variational eigenvectors
                     wfmt (:, :, ias, :, j) = 0.d0
                     i = 0
                     Do ispn = 1, nspinor
                        Do ist = 1, nstfv
                           i = i + 1
                           zt1 = evecsv (i, j)
                           If (Abs(dble(zt1))+Abs(aimag(zt1)) .Gt. &
                          & input%groundstate%epsocc) Then
                              If ( .Not. done(ist)) Then
                                 Call wavefmt &
                                & (input%groundstate%lradstep, &
                                & input%groundstate%lmaxvr, is, ia, &
                                & ngp, apwalm, evecfv(:, ist), lmmaxvr, &
                                & wfmt1)
! convert from spherical harmonics to spherical coordinates
                                 Call zgemm ('N', 'N', lmmaxvr, &
                                & nrcmt(is), lmmaxvr, zone, zbshtvr, &
                                & lmmaxvr, wfmt1, lmmaxvr, zzero, &
                                & wfmt2(:, :, ist), lmmaxvr)
                                 done (ist) = .True.
                              End If
! add to spinor wavefunction
                              Call zaxpy (n, zt1, wfmt2(:, :, ist), 1, &
                             & wfmt(:, :, ias, ispn, j), 1)
                           End If
! end loop over first-variational states
                        End Do
! end loop over spin
                     End Do
                  Else
! spin-unpolarised wavefunction
                     Call wavefmt (input%groundstate%lradstep, &
                    & input%groundstate%lmaxvr, is, ia, ngp, apwalm, &
                    & evecfv(:, j), lmmaxvr, wfmt1)
! convert from spherical harmonics to spherical coordinates
                     Call zgemm ('N', 'N', lmmaxvr, nrcmt(is), lmmaxvr, &
                    & zone, zbshtvr, lmmaxvr, wfmt1, lmmaxvr, zzero, &
                    & wfmt(:, :, ias, 1, j), lmmaxvr)
                  End If
               End If
! end loop over second-variational states
            End Do
! end loops over atoms and species
         End Do
      End Do
!-----------------------------------!
!     interstitial wavefunction     !
!-----------------------------------!
      t1 = 1.d0 / Sqrt (omega)
      Do j = 1, nstsv
         wfir (:, :, j) = 0.d0
         If (( .Not. tocc) .Or. ((tocc) .And. (evalsvp(j) .Lt. &
        & efermi))) Then
            If (input%groundstate%tevecsv) Then
! generate spinor wavefunction from second-variational eigenvectors
               i = 0
               Do ispn = 1, nspinor
                  Do ist = 1, nstfv
                     i = i + 1
                     zt1 = evecsv (i, j)
                     If (Abs(dble(zt1))+Abs(aimag(zt1)) .Gt. &
                    & input%groundstate%epsocc) Then
                        zt1 = t1 * zt1
                        Do igp = 1, ngp
                           ifg = igfft (igpig(igp))
                           wfir (ifg, ispn, j) = wfir (ifg, ispn, j) + &
                          & zt1 * evecfv (igp, ist)
                        End Do
                     End If
                  End Do
               End Do
            Else
! spin-unpolarised wavefunction
               Do igp = 1, ngp
                  ifg = igfft (igpig(igp))
                  wfir (ifg, 1, j) = t1 * evecfv (igp, j)
               End Do
            End If
! Fourier transform wavefunction to real-space
            Do ispn = 1, nspinor
               Call zfftifc (3, ngrid, 1, wfir(:, ispn, j))
            End Do
         End If
      End Do
      Deallocate (done, wfmt1, wfmt2)
      Return
End Subroutine
!EOC
