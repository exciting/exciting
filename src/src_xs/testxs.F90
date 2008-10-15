
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine testxs
use indices
!  use modmain
!  use modxs
  implicit none
!  logical :: tabel
!  integer :: scmut(maxsymcrys,maxsymcrys),i,j

!call test_indices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$  call init0
!!$  call init1
!!$  call gensymcmut(epslat,maxsymcrys,nsymcrys,symlat,lsplsymc,vtlsymc,scmut, &
!!$     tabel)
!!$
!!$  do i=1,nsymcrys
!!$     do j=1,nsymcrys
!!$        write(*,*) i,j,scmut(i,j)
!!$     end do
!!$  end do
!!$  do i=1,nsymcrys
!!$     write(8899,'(48i3.2)') scmut(i,:nsymcrys)
!!$  end do
!!$  do i=1,nsymcrys
!!$     write(8999,'(2i5.2)') i,scmut(i,i)
!!$  end do
!!$
!!$  call findsubgrps(nsymcrys,scmut)
!!$
!!$  write (*,*) 'abelian ? ',tabel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$  integer, parameter :: n=50
!!$  complex(8) :: m(n,n),m2(n,n),m3(n,n), z
!!$  real(8) :: rnd1(n,n),rnd2(n,n),cpu0,cpu1,cpu2,cpumm,cpuzg
!!$  integer :: j,k,l,i
!!$
!!$  cpumm=0.d0
!!$  cpuzg=0.d0
!!$
!!$  call random_number(rnd1)
!!$  call random_number(rnd2)
!!$  m2(:,:)=rnd1(:,:)+(0.d0,1.d0)*rnd2(:,:)
!!$
!!$  call random_number(rnd1)
!!$  call random_number(rnd2)
!!$  m3(:,:)=rnd1(:,:)+(0.d0,1.d0)*rnd2(:,:)
!!$
!!$
!!$  m(:,:)=0.d0
!!$  z=(2.2d0,3.4d0)
!!$  do j=1,1
!!$     call cpu_time(cpu0)
!!$
!!$     do i=1,10
!!$        do l=1,n
!!$           do k=1,n
!!$              m(l,k)=z*l+k
!!$           end do
!!$        end do
!!$     end do
!!$     call cpu_time(cpu1)
!!$     cpumm=cpumm+cpu1-cpu0
!!$     write(*,*) 'elapsed time (2 loops)',cpu1-cpu0
!!$
!!$     do i=1,10
!!$        forall(l=1:n,k=1:n)
!!$           m(l,k)=z*l+k
!!$        end forall
!!$     end do
!!$     call cpu_time(cpu2)
!!$     cpuzg=cpuzg+cpu2-cpu1
!!$     write(*,*) 'elapsed time  (forall statement)',cpu2-cpu1
!!$     write(*,*)
!!$  end do
!!$
!!$  write(*,*)

end subroutine testxs
