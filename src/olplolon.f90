
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olplolon(is,ia,ngp)
use modmain
implicit none
! arguments

integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp


! local variables
complex(8)::zt
integer ias,ilo1,ilo2,l,m,lm,i,j,k
ias=idxas(ia,is)
do ilo1=1,nlorb(is)
  l=lorbl(ilo1,is)
  do ilo2=1,nlorb(is)
    if (lorbl(ilo2,is).eq.l) then
      do m=-l,l
        lm=idxlm(l,m)
        i=ngp+idxlo(lm,ilo1,ias)
        j=ngp+idxlo(lm,ilo2,ias)
        if (i.le.j) then
    
! calculate the matrix elements
          !  k=i+((j-1)*j)/2
          !  o(k)=o(k)+ololo(ilo1,ilo2,ias)
            zt= dcmplx(ololo(ilo1,ilo2,ias),0.0)
            call oupdate(zt,j,i) 
      
        end if
      end do
    end if
  end do
end do
return
end subroutine

