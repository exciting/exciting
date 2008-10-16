
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findsymgenr(nc,mt,ngen,nsgen,gen,orbgen)
  use modsym
  implicit none
  ! arguments
  integer, intent(in) :: nc
  integer, intent(in) :: mt(nc,nc)
  integer, intent(out) :: ngen
  integer, intent(out) :: nsgen(nc)
  integer, intent(out) :: gen(nc)
  integer, intent(out) :: orbgen(nc,nc)
  ! local variables
  integer :: i,is,j,js,norb(nc),orb(nc,nc),idx(nc),done(nc)
  ! set up orbits of group elements
  norb(:)=0
  norb(1)=1
  do i=1,nc
     orb(i,1)=i
     do j=2,nc
        orb(i,j)=mt(orb(i,j-1),i)
	if (orb(i,j).eq.1) then
           norb(i)=j-1
           exit
        end if
     end do
  end do
  ! sort orbits according to their number of elements
  call sortidx(nc,dble(norb),idx)
  ! add largest generator to set
  ngen=1
  gen(:)=0
  gen(ngen)=orb(idx(nc),1)
  nsgen(:)=0
  nsgen(1)=norb(idx(nc))
  orbgen(:,:)=0
  orbgen(1,:)=orb(idx(nc),:nsgen(1))
  do i=nc-1,1,-1
     is=idx(i)
     do j=nc,i+1,-1
       js=idx(j)
        ! discard orbit if one of its elements is equal to previous generators
        if (any(orb(js,:norb(js)).eq.orb(is,1))) goto 10
        ! discard trivial generator
        if (orb(is,1).eq.1) goto 10
     end do
     ! add new generator to set
     ngen=ngen+1
     nsgen(ngen)=norb(is)
     gen(ngen)=orb(is,1)
     orbgen(ngen,:nsgen(ngen))=orb(is,:nsgen(ngen))
     10 continue
  end do
  ! check if orbits cover the symmetry group
  done(:)=0
  do i=1,ngen
     done(orbgen(i,:nsgen(i)))=done(orbgen(i,:nsgen(i)))+1
  end do
  if (any(done.eq.0)) then
     write(*,*)
     write(*,'("Error(findsymgenr): Generators do not cover the symmetry &
       &group")')
     write(*,*)
!!!!!!!!     stop
  end if
end subroutine findsymgenr
