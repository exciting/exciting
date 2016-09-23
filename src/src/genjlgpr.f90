!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genjlgpr (lmax, gpc, jlgpr)
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: gpc (ngvec)
      Real (8), Intent (Out) :: jlgpr (0:lmax, ngvec, nspecies)
! local variables
      Integer :: is, ig
      Real (8) :: t1
      Do is = 1, nspecies
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) SHARED(ngvec,gpc,rmt,lmax,jlgpr,is) PRIVATE(ig,t1)
!$OMP DO
#endif
         Do ig = 1, ngvec
            t1 = gpc (ig) * rmt (is)
            Call sbessel (lmax, t1, jlgpr(:, ig, is))
         End Do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL 
#endif

      End Do
      Return
End Subroutine
