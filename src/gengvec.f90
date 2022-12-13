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
      use sorting, only: sortidx
      use modmpi, only: mpiglobal
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
      End Do
! Fourier transform index
      call genigfft( ngrid, ivg, ngrtot, igfft)
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

!> generate the index map from a uniform grid of G-vectors
!> to the corresponding FFT grid
!>
!> For a given uniform grid of \(N_1 \times N_2 \times N_3\) \({\bf G}\)-vectors given by
!> \[ {\bf G}_n = {\bf G}_{(n_1,n_2,n_3)} = {\bf B} \cdot (n_1,n_2,n_3)^\top \]
!> with \(n_i = -\lfloor N_i/2 \rfloor, \dots, \lfloor N_i/2 \rfloor - 1\), 
!> for each \(n = (n_1,n_2,n_3)\) `igfft(n)` contains the index \(m\) of the point in the
!> corresponding FFT grid, given by
!> \[ m = n_3'\, N_1 N_2 + n_2'\, N_1 + n_3' + 1 \;, \]
!> with
!> \[ n_i' = \begin{cases} n_i & \text{ if } n_i \geq 0 \\ N_i + n_i & \text{ if } n_i < 0 \end{cases} \;. \]
subroutine genigfft( ngrid, ivg, ng, igfft)
  !> G-vector grid dimensions
  integer, intent(in) :: ngrid(3)
  !> integer components of G-vectors (lattice coordinates)
  integer, intent(in) :: ivg(3,*)
  !> number of G-vectors
  integer, intent(in) :: ng
  !> index map from G-vector index to FFT grid point index
  integer, intent(out) :: igfft(*)

  integer :: ig, i1, i2, i3, j1, j2, j3

  do ig = 1, ng
    i1 = ivg(1,ig)
    i2 = ivg(2,ig)
    i3 = ivg(3,ig)
    if( i1 >= 0) then
       j1 = i1
    else
       j1 = ngrid(1) + i1
    end if
    if( i2 >= 0) then
       j2 = i2
    else
       j2 = ngrid(2) + i2
    end if
    if( i3 >= 0) then
       j3 = i3
    else
       j3 = ngrid(3) + i3
    end if
    if( (j1 >= 0) .and. (j1 <= ngrid(1)) .and. &
        (j2 >= 0) .and. (j2 <= ngrid(2)) .and. &
        (j3 >= 0) .and. (j3 <= ngrid(3))) then
      igfft(ig) = j3*ngrid(2)*ngrid(1) + j2*ngrid(1) + j1 + 1
    else
      igfft(ig) = -1
    end if
  end do
end subroutine
