!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqkgir (iq, ik, igq)
      Use modmain
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik, igq
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqkgir'
      Integer :: ikq, ig, ig1, ig2, ig3, igk0, igk, iv (3), iv1 (3), &
     & iv3 (3), ivu (3)
      Integer, Allocatable :: aigk0 (:), aigk (:)
  ! grid index for k+q point
      ikq = ikmapikq (ik, iq)
      Allocate (aigk0(ngkmax0), aigk(ngkmax))
  ! positive umklapp G-vector
      ivu (:) = Nint (vkl0(:, ik)+vql(:, iq)-vkl(:, ikq))
  ! precalculate for speedup
      aigk0 (:) = igkig0 (:, 1, ik)
      aigk (:) = igkig (:, 1, ikq)
      ig3 = igqig (igq, iq)
      iv3 (:) = ivg (:, ig3)
      Do igk0 = 1, ngk0 (1, ik)
         ig1 = aigk0 (igk0)
         iv1 (:) = ivg (:, ig1) + iv3 (:)
         Do igk = 1, ngk (1, ikq)
            ig2 = aigk (igk)
        ! umklapp of k+q vector included
            iv (:) = iv1 (:) - (ivg(:, ig2)-ivu(:))
            ig = ivgig (iv(1), iv(2), iv(3))
            xihir (igk0, igk) = cfunig (ig)
         End Do
      End Do
      Deallocate (aigk0, aigk)
End Subroutine ematqkgir
