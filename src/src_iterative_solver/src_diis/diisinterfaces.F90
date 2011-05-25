
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Module diisinterfaces
!
      Implicit None
      Complex (8) zdotc
      Real (8) :: dlamch
      External zdotc, dlamch
      Real (8), External :: dznrm2
#ifdef USEDIISIF
      Interface
!
!
         Subroutine DIISseceqnfv (ik, ispn, apwalm, vgpc, evalfv, &
        & evecfv)
            Use modmain, Only: nstfv, vkl, ngk, igkig, nmat, vgkl, &
           & timemat, npmat, apwordmax, lmmaxapw, natmtot, nkpt, &
           & nmatmax, nspnfv, timefv, ngkmax
            Integer, Intent (In) :: ik
            Integer, Intent (In) :: ispn
            Real (8), Intent (In) :: vgpc (3, ngkmax)
            Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, &
           & lmmaxapw, natmtot)
            Real (8), Intent (Inout) :: evalfv (nstfv, nspnfv)
            Complex (8), Intent (Inout) :: evecfv (nmatmax, nstfv, &
           & nspnfv)
         End Subroutine
      End Interface
!
      Interface
!
!
         Subroutine seceqfvprecond (n, system, X, w, evalfv, evecfv)
            Use modmain, Only: nmatmax, nstfv
            Use modfvsystem
            Implicit None
            Type (evsystem) :: system
!
            Integer, Intent (In) :: n
!
            Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
            Real (8), Intent (Out) :: evalfv (nstfv), w (nmatmax)
            Complex (8), Intent (Out) :: X (nmatmax, nmatmax)
!
         End Subroutine seceqfvprecond
      End Interface
!
      Interface
!
!
         Subroutine prerotate_preconditioner (n, m, h, P)
            Use modmain, Only: nstfv, nmatmax
            Implicit None
            Integer, Intent (In) :: n, m
            Complex (8), Intent (In) :: h (n, n)
            Complex (8), Intent (Inout) :: P (nmatmax, nmatmax)
!
         End Subroutine prerotate_preconditioner
      End Interface
!
      Interface
!
!
         Subroutine precondspectrumupdate (n, m, hamilton, overlap, P, &
        & w)
            Use modmain, Only: nstfv, nmatmax
            Implicit None
            Integer, Intent (In) :: n, m
            Complex (8), Intent (In) :: hamilton (n, n), overlap (n, n)
            Complex (8), Intent (In) :: P (nmatmax, nmatmax)
            Real (8), Intent (Inout) :: w (nmatmax)
         End Subroutine precondspectrumupdate
      End Interface
      Interface
!
!
         Subroutine readprecond (ik, n, X, w)
            Use modmain
            Use modmpi
            Integer, Intent (In) :: n, ik
            Complex (8), Intent (Out) :: X (nmatmax, nmatmax)
            Real (8), Intent (Out) :: w (nmatmax)
         End Subroutine readprecond
      End Interface
      Interface
!
!
         Subroutine writeprecond (ik, n, X, w)
            Use modmain
            Use modmpi
            Integer, Intent (In) :: n, ik
            Complex (8), Intent (In) :: X (nmatmax, nmatmax)
            Real (8), Intent (In) :: w (nmatmax)
         End Subroutine writeprecond
      End Interface
      Interface
!
!
         Subroutine setuphsvect (n, m, system, evecfv, ldv, h, s)
!
            Use modmain, Only: nstfv, zone, zzero
            Use modfvsystem
            Implicit None
            Integer, Intent (In) :: n, m, ldv
            Type (evsystem) :: system
            Complex (8), Intent (In) :: evecfv (ldv, m)
            Complex (8), Intent (Out) :: h (n, m), s (n, m)
!
         End Subroutine setuphsvect
      End Interface
      Interface
!
!
         Subroutine rayleighqotient (n, m, evecfv, h, s, evalfv)
!
            Implicit None
            Integer, Intent (In) :: n, m
            Complex (8), Intent (In) :: h (n, m), s (n, m), evecfv (n, &
           & m)
            Real (8), Intent (Out) :: evalfv (m)
         End Subroutine rayleighqotient
      End Interface
      Interface
!
!
         Subroutine residualvectors (n, iunconverged, h, s, evalfv, r, &
        & rnorms)
            Use modmain, Only: nmatmax, nstfv
            Implicit None
            Integer, Intent (In) :: n, iunconverged
       !packed ut
            Complex (8), Intent (In) :: h (n, nstfv), s (n, nstfv)
            Complex (8), Intent (Out) :: r (n, nstfv)
            Real (8), Intent (In) :: evalfv (nstfv)
            Real (8), Intent (Out) :: rnorms (nstfv)
         End Subroutine residualvectors
      End Interface
      Interface
!
!
         Subroutine calcupdatevectors (n, iunconverged, P, w, r, &
        & evalfv, evecfv, phi)
            Use modmain, Only: nstfv, zzero, zone, nmatmax
!
            Implicit None
            Integer, Intent (In) :: n, iunconverged
            Complex (8), Intent (In) :: P (nmatmax, nmatmax)
            Complex (8), Intent (In) :: r (n, nstfv), evecfv (n, nstfv)
            Complex (8), Intent (Out) :: phi (n, nstfv)
            Real (8), Intent (In) :: w (nmatmax), evalfv (nstfv)
         End Subroutine calcupdatevectors
!
      End Interface
      Interface
!
!
         Subroutine diisupdate (idiis, icurrent, iunconverged, n, h, s, &
        & trialvec, evalfv, evecfv, infodiisupdate)
!
            Use modmain, Only: nstfv, zone, zzero
            Use sclcontroll, Only: diismax
            Implicit None
            Integer, Intent (In) :: idiis, icurrent, iunconverged, n
            Complex (8), Intent (In) :: h (n, nstfv, diismax)
            Complex (8), Intent (In) :: s (n, nstfv, diismax), trialvec &
           & (n, nstfv, diismax)
            Real (8), Intent (In) :: evalfv (nstfv, diismax)
            Complex (8), Intent (Out) :: evecfv (n, nstfv)
            Integer, Intent (Out) :: infodiisupdate
!
         End Subroutine diisupdate
      End Interface
      Interface
!
!
         Subroutine normalize (n, m, overlap, evecfv, ldv)
            Use modmain, Only: zone, zzero
!
            Implicit None
            Integer, Intent (In) :: n, m, ldv
            Complex (8), Intent (In) :: overlap (n, n)
            Complex (8), Intent (Out) :: evecfv (ldv, m)
         End Subroutine
      End Interface
      Interface
!
!
         Subroutine solvediis (m, Pmatrix, Qmatrix, c)
            Implicit None
            Integer, Intent (In) :: m
            Real (8), Intent (In) :: Pmatrix (m, m), Qmatrix (m, m)
            Real (8), Intent (Out) :: c (m)
         End Subroutine solvediis
      End Interface
      Interface
!
!
         Subroutine exactupdatevectors (n, iunconverged, hamilton, &
        & overlap, r, rhizvalue, eigenvector, trialvecs)
!calculate update equation with linsolver
!
!solvefor dA:  dA=(H-e*S)\R
!
! dA 	Update step to zero residual
! H 	Hamilton
! S 	Overlap
! e 	Rhitz Value
! R 	Residual
!
!trialvecs=eigenvector+dA
            Use modfvsystem
            Use modmain, Only: zone
            Integer, Intent (In) :: n, iunconverged
            Type (HermitianMatrix), Intent (In) :: hamilton, overlap
            Complex (8), Intent (In) :: r (n, iunconverged), &
           & eigenvector (n, iunconverged)
            Real (8), Intent (In) :: rhizvalue (iunconverged)
            Complex (8), Intent (Out) :: trialvecs (n, iunconverged)
!
         End Subroutine
      End Interface
#endif
End Module diisinterfaces
