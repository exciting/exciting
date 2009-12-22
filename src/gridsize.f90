!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gridsize
! !INTERFACE:
!
!
Subroutine gridsize
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Finds the ${\bf G}$-vector grid which completely contains the vectors with
!   $G<G_{\rm max}$ and is compatible with the FFT routine. The optimal sizes
!   are given by
!   $$ n_i=\frac{G_{\rm max}|{\bf a}_i|}{\pi}+1, $$
!   where ${\bf a}_i$ is the $i$th lattice vector.
!
! !REVISION HISTORY:
!   Created July 2003 (JKD)
!EOP
!BOC
      Implicit None
! find optimal grid size for potential and density
      ngrid (:) = Int (input%groundstate%gmaxvr*&
     & Sqrt(input%structure%crystal%basevect(1, :)**2+&
     & input%structure%crystal%basevect(2, :)**2+&
     & input%structure%crystal%basevect(3, :)**2)/pi) + 1
! find next largest FFT-compatible grid size
      Call nfftifc (ngrid(1))
      Call nfftifc (ngrid(2))
      Call nfftifc (ngrid(3))
      If ((ngrid(1) .Le. 0) .Or. (ngrid(2) .Le. 0) .Or. (ngrid(3) .Le. &
     & 0)) Then
         Write (*,*)
         Write (*, '("Error(gridsize): invalid ngrid : ", 3I8)') ngrid
         Write (*,*)
         Stop
      End If
! total number of points in grid
      ngrtot = ngrid (1) * ngrid (2) * ngrid (3)
! determine integer ranges for grid
      intgv (:, 1) = ngrid (:) / 2 - ngrid (:) + 1
      intgv (:, 2) = ngrid (:) / 2
      Return
End Subroutine
!EOC
