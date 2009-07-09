


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine hmllolon(hamilton, is, ia, ngp)
use modmain
use modinput
use modfvsystem
implicit none
! arguments
type(HermiteanMatrix)::hamilton
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp


! local variables
integer::ias, ilo1, ilo2, i, j, k
integer::l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
complex(8) zsum
ias=idxas(ia, is)
do ilo1=1, nlorb(is)
  l1=lorbl(ilo1, is)
  do m1=-l1, l1
    lm1=idxlm(l1, m1)
    i=ngp+idxlo(lm1, ilo1, ias)
    do ilo2=1, nlorb(is)
      l3=lorbl(ilo2, is)
      do m3=-l3, l3
	lm3=idxlm(l3, m3)
	j=ngp+idxlo(lm3, ilo2, ias)
	if (i.le.j) then
	  zsum=0.d0
	  do l2=0, input%groundstate%lmaxvr
	    if (mod(l1+l2+l3, 2).eq.0) then
	      do m2=-l2, l2
		lm2=idxlm(l2, m2)
		zsum=zsum+gntyry(lm1, lm2, lm3)*hlolo(ilo1, ilo2, lm2, ias)
	      end do
	    end if
	  end do

! calculate the matrix elements


       call Hermiteanmatrix_indexedupdate(hamilton, j, i, zsum)

	end if
      end do
    end do
  end do
end do
return
end subroutine
