
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olplolo(tapp,is,ia,ngp,v,o)
use modmain
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: v(nmatmax)
complex(8), intent(inout) :: o(*)
! local variables
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
          if (tapp) then
! apply the overlap operator to v
            o(i)=o(i)+ololo(ilo1,ilo2,ias)*v(j)
            if (i.ne.j) o(j)=o(j)+ololo(ilo1,ilo2,ias)*v(i)
          else
! calculate the matrix elements
            k=i+((j-1)*j)/2
            o(k)=o(k)+ololo(ilo1,ilo2,ias)
          end if
        end if
      end do
    end if
  end do
end do
return
end subroutine

