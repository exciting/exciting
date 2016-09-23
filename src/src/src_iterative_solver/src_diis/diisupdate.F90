
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine diisupdate (idiis, icurrent, iunconverged, n, h, s, &
& trialvec, evalfv, evecfv, infodiisupdate)
      Use modmain, Only: nstfv, zone, zzero
      Use diisinterfaces
      Use sclcontroll, Only: recalculate_preconditioner
      Use sclcontroll, Only: diismax, maxdiisspace
      Implicit None
      Integer, Intent (In) :: idiis, icurrent, iunconverged, n
      Complex (8), Intent (In) :: h (n, nstfv, diismax)
      Complex (8), Intent (In) :: s (n, nstfv, diismax), trialvec (n, &
     & nstfv, diismax)
      Real (8), Intent (In) :: evalfv (nstfv, diismax)
      Complex (8), Intent (Out) :: evecfv (n, nstfv)
      Integer, Intent (Out) :: infodiisupdate
      Integer :: isubspace
      Logical :: lin
!
      Complex (8), Allocatable :: p (:, :)
      Real (8) :: nrm
      Integer :: i, j, ir, is
      Integer, Allocatable :: idx (:)
      Real (8), Allocatable :: Pmatrix (:, :), Qmatrix (:, :)
      Real (8), Allocatable :: c (:)
      Real (8) :: residnorm2
      Complex (8) :: z
      lin = .True.
!
      isubspace = Min (idiis, maxdiisspace)
      Allocate (p(n, isubspace))
      Allocate (Pmatrix(isubspace+1, isubspace+1), Qmatrix(isubspace+1, &
     & isubspace+1))
      Allocate (c(isubspace+1))
      Allocate (idx(isubspace))
      Pmatrix = 0.0
      Qmatrix = 0.0
      p = 0
      c = 0
      Do i = 1, iunconverged
     !calculate residuals
         Do j = 1, isubspace
            Call zcopy (n, h(1, i, j), 1, p(1, j), 1)
            z = cmplx (-evalfv(i, j), 0)
            Call zaxpy (n, z, s(1, i, j), 1, p(1, j), 1)
         End Do
         residnorm2 = dble (zdotc(n, p(1, icurrent), 1, p(1, icurrent), &
        & 1))
!
         Do ir = 1, isubspace
            Do is = 1, isubspace
               Pmatrix (is, ir) = dble (zdotc(n, p(1, is), 1, p(1, ir), &
              & 1)) / residnorm2
               If (dble(Pmatrix(is, ir)) .Lt. 1.e-4) Then
                  Write (889,*) "ir, is, p(1, i, ir)", ir, is, p (1, &
                 & is)
               End If
            End Do
         End Do
         If ( .Not. lin) Then
            Do ir = 1, isubspace
               Do is = 1, isubspace
                  Qmatrix (is, ir) = dble (zdotc(n, trialvec(1, i, is), &
                 & 1, s(1, i, ir), 1))
                  If (dble(Qmatrix(is, ir)) .Lt. 1.e-4) Then
                 ! write(*,*)"warning Qmatrix(is,ir)).lt.1.e-4 in diisupdate"
                 !   write(888,*)"ir,is,trialvec(1,i,is),s(1,i,ir)",ir,is,trialvec(:,i,is),"s\n\n",s(:,i,ir)
                  End If
!
               End Do
            End Do
         End If
    ! if(i==1 )write(*,*)"Pmatrix",Pmatrix
         If (lin) Then
            Call solvediislin (isubspace, Pmatrix, Qmatrix, c)
            Call sortidx (isubspace,-Abs(c), idx)
         Else
            Call solvediis (isubspace, Pmatrix, Qmatrix, c)
         End If
         If (recalculate_preconditioner .Eqv. .True.) Then
            infodiisupdate = 1
            Exit
         End If
     ! write(*,*) "c",c
!
!
         Call zcopy (n, zzero, 0, evecfv(:, i), 1)
         Do ir = 1, isubspace
            z = dcmplx (c(idx(ir)), 0.0)
            Call zaxpy (n, z, trialvec(1, i, idx(ir)), 1, evecfv(1, i), &
           & 1)
         End Do
      End Do
      Call zcopy (n*iunconverged, evecfv, 1, trialvec(1, 1, icurrent), &
     & 1)
      infodiisupdate = 0
End Subroutine diisupdate
