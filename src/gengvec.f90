!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gengvec
! !INTERFACE:
!
!
Subroutine gengvec
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates a set of ${\bf G}$-vectors used for the Fourier transform of the
!   charge density and potential and sorts them according to length. Integers
!   corresponding to the vectors in lattice coordinates are stored, as well as
!   the map from these integer coordinates to the ${\bf G}$-vector index. A map
!   from the ${\bf G}$-vector set to the standard FFT array structure is also
!   generated. Finally, the number of ${\bf G}$-vectors with magnitude less than
!   {\tt gmaxvr} is determined.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Increased number of G-vectors to ngrtot, July 2007 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ig, i1, i2, i3, j1, j2, j3, k
      Real (8) :: v (3), t1
! allocatable arrays
      Integer, Allocatable :: idx (:)
      Integer, Allocatable :: iar (:)
      Real (8), Allocatable :: rar (:)
! allocate local arrays
      Allocate (idx(ngrtot))
      Allocate (iar(ngrtot))
      Allocate (rar(ngrtot))
! allocate global G-vector arrays
      If (allocated(ivg)) deallocate (ivg)
      Allocate (ivg(3, ngrtot))
      If (allocated(ivgig)) deallocate (ivgig)
      Allocate (ivgig(intgv(1, 1) :intgv(1, 2), intgv(2, 1) :intgv(2, &
     & 2), intgv(3, 1) :intgv(3, 2)))
      If (allocated(igfft)) deallocate (igfft)
      Allocate (igfft(ngrtot))
      If (allocated(vgc)) deallocate (vgc)
      Allocate (vgc(3, ngrtot))
      If (allocated(gc)) deallocate (gc)
      Allocate (gc(ngrtot))
      ig = 0
      Do i1 = intgv (1, 1), intgv (1, 2)
         Do i2 = intgv (2, 1), intgv (2, 2)
            Do i3 = intgv (3, 1), intgv (3, 2)
               v (:) = dble (i1) * bvec (:, 1) + dble (i2) * bvec (:, &
              & 2) + dble (i3) * bvec (:, 3)
               t1 = v (1) ** 2 + v (2) ** 2 + v (3) ** 2
               ig = ig + 1
! map from G-vector to (i1,i2,i3) index
               ivg (1, ig) = i1
               ivg (2, ig) = i2
               ivg (3, ig) = i3
! length of each G-vector
               gc (ig) = Sqrt (t1)
            End Do
         End Do
      End Do
! sort by vector length
      Call sortidx (ngrtot, gc, idx)
! re-order arrays
      Do ig = 1, ngrtot
         rar (ig) = gc (ig)
      End Do
      Do ig = 1, ngrtot
         gc (ig) = rar (idx(ig))
      End Do
      Do k = 1, 3
         Do ig = 1, ngrtot
            iar (ig) = ivg (k, ig)
         End Do
         Do ig = 1, ngrtot
            ivg (k, ig) = iar (idx(ig))
         End Do
      End Do
      ivgig (:, :, :) = 0
      Do ig = 1, ngrtot
         i1 = ivg (1, ig)
         i2 = ivg (2, ig)
         i3 = ivg (3, ig)
! map from (i1,i2,i3) index to G-vector
         ivgig (i1, i2, i3) = ig
! assign G-vectors to global array
         vgc (:, ig) = dble (i1) * bvec (:, 1) + dble (i2) * bvec (:, &
        & 2) + dble (i3) * bvec (:, 3)
! Fourier transform index
         If (i1 .Ge. 0) Then
            j1 = i1
         Else
            j1 = ngrid (1) + i1
         End If
         If (i2 .Ge. 0) Then
            j2 = i2
         Else
            j2 = ngrid (2) + i2
         End If
         If (i3 .Ge. 0) Then
            j3 = i3
         Else
            j3 = ngrid (3) + i3
         End If
         igfft (ig) = j3 * ngrid (2) * ngrid (1) + j2 * ngrid (1) + j1 &
        & + 1
      End Do
! find the number of vectors with G < gmaxvr
      ngvec = 1
      Do ig = ngrtot, 1, - 1
         If (gc(ig) .Lt. input%groundstate%gmaxvr) Then
            ngvec = ig
            Go To 10
         End If
      End Do
10    Continue
      Deallocate (idx, iar, rar)
      Return
End Subroutine
!EOC
