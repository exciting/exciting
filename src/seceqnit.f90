!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine seceqnit (nmatp, ngp, igpig, vpl, vgpl, vgpc, apwalm, &
& evalfv, evecfv)
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: nmatp
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (In) :: vgpl (3, ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
! local variables
      Integer :: is, ia, it, i
      Integer :: ist, jst
      Real (8) :: ts1, ts0
      Real (8) :: t1
      Complex (8) zt1
! allocatable arrays
      Complex (8), Allocatable :: h (:)
      Complex (8), Allocatable :: o (:, :)
! external functions
      Complex (8) zdotc
      External zdotc
      Call timesec (ts0)
      Allocate (o(nmatp, nstfv))
      If ((iscl .Ge. 2) .Or. (task .Eq. 1) .Or. (task .Eq. 3)) Then
! read in the eigenvalues/vectors from file
         Call getevalfv (vpl, evalfv)
         Call getevecfv (vpl, vgpl, evecfv)
      Else
! initialise the eigenvectors to canonical basis vectors
         evecfv (:, :) = 0.d0
         Do ist = 1, nstfv
            evecfv (ist, ist) = 1.d0
         End Do
      End If
! start iteration loop
      Do it = 1, nseqit
! begin parallel loop over states
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(h,is,ia,t1,i)
!$OMP DO
         Do ist = 1, nstfv
            Allocate (h(nmatp))
! operate with H and O on the current vector
            h (:) = 0.d0
            o (:, ist) = 0.d0
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  Call hmlaa (.True., is, ia, ngp, apwalm, evecfv(:, &
                 & ist), h)
                  Call hmlalo (.True., is, ia, ngp, apwalm, evecfv(:, &
                 & ist), h)
                  Call hmllolo (.True., is, ia, ngp, evecfv(:, ist), h)
                  Call olpaa (.True., is, ia, ngp, apwalm, evecfv(:, &
                 & ist), o(:, ist))
                  Call olpalo (.True., is, ia, ngp, apwalm, evecfv(:, &
                 & ist), o(:, ist))
                  Call olplolo (.True., is, ia, ngp, evecfv(:, ist), &
                 & o(:, ist))
               End Do
            End Do
            Call hmlistl (.True., ngp, igpig, vgpc, evecfv(:, ist), h)
            Call olpistl (.True., ngp, igpig, evecfv(:, ist), o(:, &
           & ist))
! normalise
            t1 = dble (zdotc(nmatp, evecfv(:, ist), 1, o(:, ist), 1))
            If (t1 .Gt. 0.d0) Then
               t1 = 1.d0 / Sqrt (t1)
               Do i = 1, nmatp
                  evecfv (i, ist) = t1 * evecfv (i, ist)
                  h (i) = t1 * h (i)
                  o (i, ist) = t1 * o (i, ist)
               End Do
            End If
! estimate the eigenvalue
            evalfv (ist) = dble (zdotc(nmatp, evecfv(:, ist), 1, h, 1))
! subtract the gradient of the Rayleigh quotient from the eigenvector
            t1 = evalfv (ist)
            Do i = 1, nmatp
               evecfv (i, ist) = evecfv (i, ist) - tauseq * &
              & (h(i)-t1*o(i, ist))
            End Do
! normalise
            o (:, ist) = 0.d0
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  Call olpaa (.True., is, ia, ngp, apwalm, evecfv(:, &
                 & ist), o(:, ist))
                  Call olpalo (.True., is, ia, ngp, apwalm, evecfv(:, &
                 & ist), o(:, ist))
                  Call olplolo (.True., is, ia, ngp, evecfv(:, ist), &
                 & o(:, ist))
               End Do
            End Do
            Call olpistl (.True., ngp, igpig, evecfv(:, ist), o(:, &
           & ist))
            t1 = dble (zdotc(nmatp, evecfv(:, ist), 1, o(:, ist), 1))
            If (t1 .Gt. 0.d0) Then
               t1 = 1.d0 / Sqrt (t1)
               Do i = 1, nmatp
                  evecfv (i, ist) = t1 * evecfv (i, ist)
                  o (i, ist) = t1 * o (i, ist)
               End Do
            End If
            Deallocate (h)
! end parallel loop over states
         End Do
!$OMP END DO
!$OMP END PARALLEL
! perform Gram-Schmidt orthonormalisation
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jst,zt1,t1,i)
!$OMP DO ORDERED
         Do ist = 1, nstfv
!$OMP ORDERED
            Do jst = 1, ist - 1
               zt1 = - zdotc (nmatp, evecfv(:, jst), 1, o(:, ist), 1)
               Call zaxpy (nmatp, zt1, evecfv(:, jst), 1, evecfv(:, &
              & ist), 1)
               Call zaxpy (nmatp, zt1, o(:, jst), 1, o(:, ist), 1)
            End Do
!$OMP END ORDERED
! normalise
            t1 = dble (zdotc(nmatp, evecfv(:, ist), 1, o(:, ist), 1))
            If (t1 .Gt. 0.d0) Then
               t1 = 1.d0 / Sqrt (t1)
               Do i = 1, nmatp
                  evecfv (i, ist) = t1 * evecfv (i, ist)
                  o (i, ist) = t1 * o (i, ist)
               End Do
            End If
         End Do
!$OMP END DO
!$OMP END PARALLEL
! end iteration loop
      End Do
      Deallocate (o)
      Call timesec (ts1)
!$OMP CRITICAL
      timefv = timefv + ts1 - ts0
!$OMP END CRITICAL
      Return
End Subroutine
