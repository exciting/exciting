
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public 
! License. See the file COPYING for license details.

module summations
  implicit none
contains

  subroutine doublesummation_simple_cz(z,za,zm,zb,alpha,beta,tbal)
    implicit none
    complex(8), intent(inout) :: z(:,:)
    complex(8), intent(in) :: za(:,:),zm(:,:),zb(:,:)
    complex(8), intent(in) :: alpha,beta
    logical, intent(in) :: tbal
    ! local variables
    integer :: n1,n2,m1,m2,r1,r2,s1,s2
    complex(8), parameter :: zzero=cmplx(0.d0,0.d0,8), zone=cmplx(1.d0,0.d0,8)
    complex(8), allocatable :: zw(:,:)
    r1=size(z,1);  r2=size(z,2)
    m1=size(za,1); n1=size(za,2)
    m2=size(zb,1); n2=size(zb,2)
    s1=size(zm,1); s2=size(zm,2)
    ! check output-array consistency
    if ((r1.ne.n1).or.(r2.ne.n2)) then
       write(*,*)
       write(*,'("Error(doublesummation_simple_cz): inconsistent output&
            & array size: (r1,r2,n1,n2)",4i8)') r1,r2,n1,n2
       write(*,*)
       stop
    end if
    ! multiplication consistency
    if ((s1.ne.m1).or.(s2.ne.m2)) then
       write(*,*)
       write(*,'("Error(doublesummation_simple_cz): inconsistent multiplication&
            & array size: (s1,s2,m1,m2)",4i8)') s1,s2,m1,m2
       write(*,*)
       stop
    end if
    ! choose the order of the two matrix multiplications in a way such that the
    ! number of operations is balanced for botz multiplications
    if ((n1.ge.n2).and.tbal) then
       allocate(zw(m1,n2))
       ! calculate double summation :
       ! Z(n1,n2) = alpha * sum{m1,m2} x
       !          x  conjg(A(m1,n1)) [M(m1,m2) B(m2,n2)] + beta*Z(n1,n2)
       ! as two matrix multiplications
       ! (i)  W(m1,n2) = M(m1,m2) . B(m2,n2)
       ! (ii) Z(n1,n2)= alpha * conjg(A'(n1,m1)) . W(m1,n2) + beta*Z(n1,n2)
       ! ad (i)
       call zgemm('n','n', m1, n2, m2, zone, zm, &
            m1, zb, m2, zzero, zw, m1 )
       ! ad (ii)
       call zgemm('c','n', n1, n2, m1, alpha, za, &
            m1, zw, m1, beta, z, n1 )
    else
       allocate(zw(n1,m2))
       ! calculate double summation :
       ! Z(n1,n2) = alpha * sum{m1,m2} x
       !          x  [conjg(A(m1,n1)) M(m1,m2)] B(m2,n2) + beta*Z(n1,n2)
       ! as two matrix multiplications
       ! (i)  W(n1,m2) = conjg(A'(n1,m1)) . M(m1,m2)
       ! (ii) Z(n1,n2)= alpha * W(n1,m2) .  B(m2,n2) + beta*Z(n1,n2)
       ! ad (i)
       call zgemm('c','n', n1, m2, m1, zone, za, &
            m1, zm, m1, zzero, zw, n1 )
       ! ad (ii)
       call zgemm('n','n', n1, n2, m2, alpha, zw, &
            n1, zb, m2, beta, z, n1 )
    end if
    deallocate(zw)
  end subroutine doublesummation_simple_cz

end module summations

!///////////////////////////////////////////////////////////////////////////////

!!$program test
!!$  use summations
!!$  implicit none
!!$
!!$  integer, parameter :: n1=200,n2=40,m1=2000,m2=2000
!!$  complex(8) :: z2(n1,n2),z(n1,n2),za(m1,n1),za2(n1,m1),zb(m2,n2),zm(m1,m2)
!!$  complex(8) :: zw(m1,n2),zw2(n1,m2),zi,alpha,beta
!!$  real(8) :: rz(n1,n2),ra(m1,n1),rb(m2,n2),rm(m1,m2)
!!$  real(8) :: cpu0,cpu1,cpuzgemm,cpumatmul
!!$
!!$  zi=cmplx(0.d0,1.d0,8)
!!$  alpha=1.d0
!!$  beta=0.d0
!!$
!!$  za=0.d0
!!$  zb=0.d0
!!$  zm=0.d0
!!$  call random_number(ra)
!!$  call random_number(rb)
!!$  call random_number(rm)
!!$  za=za+ra/dble(m1)
!!$  zb=zb+rb/dble(m2)
!!$  zm=zm+rm
!!$  call random_number(ra)
!!$  call random_number(rb)
!!$  call random_number(rm)
!!$  za=za+zi*ra/dble(m1)
!!$  zb=zb+zi*rb/dble(m2)
!!$  zm=zm+zi*rm
!!$
!!$  za2=transpose(za)
!!$
!!$  call cpu_time(cpu0)
!!$  call doublesummation_simple_cz(z,za,zm,zb,alpha,beta)
!!$  call cpu_time(cpu1)
!!$  cpuzgemm=cpu1-cpu0
!!$if (n1.ge.n2) then
!!$   zw=matmul(zm,zb)
!!$   z2=matmul(conjg(transpose(za)),zw)
!!$else
!!$   zw2=matmul(conjg(transpose(za)),zm)
!!$   z2=matmul(zw2,zb)
!!$end if
!!$  call cpu_time(cpu0)
!!$  cpumatmul=cpu0-cpu1
!!$  write(111,*) z
!!$  write(222,*) z2
!!$  write(333,*) z-z2
!!$
!!$write(*,'(a,f12.3)') 'CPU-time (zgemm) :',cpuzgemm
!!$write(*,'(a,f12.3)') 'CPU-time (matmul):',cpumatmul
!!$write(*,*) 'max deviation:',maxval(abs(z-z2))
!!$
!!$end program test
