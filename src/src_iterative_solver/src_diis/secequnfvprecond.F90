
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine seceqfvprecond (n, system, X, w, evalfv, evecfv)
      Use modmain, Only: nmatmax, nstfv
      Use modfvsystem
      Use diisinterfaces
      Implicit None
      Integer, Intent (In) :: n
      Type (evsystem) :: system
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
      Real (8), Intent (Out) :: evalfv (nstfv), w (nmatmax)
      Complex (8), Intent (Out) :: X (nmatmax, nmatmax)
!
  !local var
      Integer :: i
  !workarrays for lapaack
      Integer :: info, lwork
      Integer :: iwork (5*n)
      Integer :: ifail (n), mfound
      Real (8) :: v, abstol
      Real (8) :: rwork (7*n)
      Complex (8) :: work (2*n)
      lwork = 2 * n
      abstol = 2.d0 * dlamch ('S')
      Call zhegvx (1, 'V', 'A', 'U', n, system%hamilton%za, n, &
     & system%overlap%za, n, v, v, 1, nstfv, abstol, mfound, w, X, &
     & nmatmax, work, lwork, rwork, iwork, ifail, info)
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(seceqnfv): diagonalisation failed")')
         Write (*, '(" ZHPGVX returned INFO = ", I8)') info
         If (info .Gt. n) Then
            i = info - n
            Write (*, '(" The leading minor of the overlap matrix of or&
           &der ", I8)') i
            Write (*, '("  is not positive definite")')
            Write (*, '(" Order of overlap matrix : ", I8)') n
            Write (*,*)
         End If
         Stop
      End If
      Call dcopy (nstfv, w(1), 1, evalfv(1), 1)
      Call zcopy (nstfv*nmatmax, X, 1, evecfv(1, 1), 1)
!
!
End Subroutine seceqfvprecond
