
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine solvediislin (m, Pmatrix, Qmatrix, c)
      Use diisinterfaces
      Use sclcontroll, Only: recalculate_preconditioner
      Implicit None
!
      Integer, Intent (In) :: m
!
      Real (8), Intent (Inout) :: Pmatrix (m+1, m+1), Qmatrix (m+1, &
     & m+1)
      Real (8), Intent (Out) :: c (m+1)
      Integer :: ipiv (m+1), info = 0, i, RANK, lwork = - 1
      Real (8), Allocatable :: WORK (:)
      Real (8) :: tmp (1)
      Pmatrix (m+1, m+1) = 0.0
      c (m+1) = - 1
      Do i = 1, m
         Pmatrix (i, m+1) = - 1.0
         Pmatrix (m+1, i) = - 1.0
         c (i) = 0
      End Do
      Call DGESV (m+1, 1, Pmatrix, m+1, ipiv, c, m+1, info)
 ! call DGELSY( m+1, m+1, 1, Pmatrix, m+1, c, m+1, IPIV,.10, RANK,&
 !           tmp, LWORK, INFO )
 !   LWORK=tmp(1)
 !   allocate(WORK(LWORK))
 !   INFO=1
 ! call DGELSY( m+1, m+1, 1, Pmatrix, m+1, c, m+1, IPIV,.10, RANK,&
        !        WORK, LWORK, INFO )
!
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(solvediis):  failed")')
         Write (*, '(" DGESV returned INFO = ", I8)') info
#ifdef DEBUG		
         Write (775,*) (Pmatrix)
         Write (776,*) (Qmatrix)
     ! stop
         recalculate_preconditioner = .True.
#endif
      End If
!
End Subroutine solvediislin
