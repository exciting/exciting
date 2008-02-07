
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpalon(is,ia,ngp,apwalm)
use modmain
implicit none
! arguments

integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)


! local variables
integer ias,ilo,io,l,m,lm,i,j,k
complex(8) zsum
ias=idxas(ia,is)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    j=ngp+idxlo(lm,ilo,ias)
   
! calculate the matrix elements
      !k=((j-1)*j)/2
      do i=1,ngp
       ! k=k+1
        zsum=0.d0
        do io=1,apword(l,is)
          zsum=zsum+conjg(apwalm(i,io,lm,ias))*oalo(io,ilo,ias)
        end do
        !o(k)=o(k)+zsum
        call oupdate(zsum,j,i) 
        
      end do
   
  end do
end do
return
end subroutine

