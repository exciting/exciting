
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findsubgrps(n,sld,lspl,mt,sg)
  implicit none
  ! arguments
  integer, intent(in) :: n
  integer, intent(in) :: sld(n)
  integer, intent(in) :: lspl(n)
  integer, intent(in) :: mt(n,n)
  integer, intent(out) :: sg(n,*)
  ! local variables
  integer :: i,j,l,done(n),map(n)

  ! find generators and favour those with positive determinant
  j=0
  do i=1,n
     l=lspl(i)
     if (sld(l).eq.1) then
        j=j+1
        map(j)=i
     end if
  end do
  do i=1,n
     l=lspl(i)
     if (sld(l).eq.-1) then
        j=j+1
        map(j)=i
     end if
  end do
  
  do i=1,n
     write(*,*) i,map(i),sld(lspl(i))
  end do

end subroutine findsubgrps
