
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
    ! number of operations is balanced for both multiplications
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

module blaswrappers
  implicit none
contains
  subroutine zgemm_wrap(z,ta,za,tb,zb,alpha,beta)
    implicit none
    ! arguments
    complex(8), intent(inout) :: z(:,:)
    complex(8), intent(in) :: za(:,:),zb(:,:)
    character(1), intent(in) :: ta,tb
    complex(8), intent(in) :: alpha,beta
    ! local variables
    integer :: nz1,nz2,na1,na2,nb1,nb2
    nz1=size(z,1); nz2=size(z,2)
    if ((ta.eq.'n').and.(tb.eq.'n')) then
       na1=size(za,1); na2=size(za,2)
       nb1=size(zb,1); nb2=size(zb,2)
       ! z(nz1,nz2)~alpha*za(na1,na2)*zb(nb1,nb2)
       if ((na1.gt.nz1).or.(na2.ne.nb1).or.(nb2.gt.nz2)) then
          write(*,*)
          write(*,'("Error(zgemm_wrap): inconsistent array dimensions")')
          write(*,'(" Z = alpha A B + beta Z")')
          write(*,'(" Z : ",2i8)') nz1,nz2
          write(*,'(" A : ",2i8)') na1,nb2
          write(*,'(" B : ",2i8)') nb1,nb2
          write(*,*)
          stop
       end if
       ! call to BLAS routine
       call zgemm(ta,tb,na1,nb2,na2,alpha,za,na1,zb,nb1,beta,z,nz1)
    else if ((ta.eq.'n').and.(tb.ne.'n')) then
       na1=size(za,1); na2=size(za,2)
       nb1=size(zb,1); nb2=size(zb,2)
       ! z(nz1,nz2)~alpha*za(na1,na2)* j(zb(nb2,nb1))
       if ((na1.gt.nz1).or.(na2.ne.nb2).or.(nb1.gt.nz2)) then
          write(*,*)
          write(*,'("Error(zgemm_wrap): inconsistent array dimensions")')
          write(*,'(" Z = alpha A j(B) + beta Z")')
          write(*,'(" Z    : ",2i8)') nz1,nz2
          write(*,'(" A    : ",2i8)') na1,nb2
          write(*,'(" j(B) : ",2i8)') nb2,nb1
          write(*,*)
          stop
       end if
       ! call to BLAS routine
       call zgemm(ta,tb,na1,nb1,na2,alpha,za,na1,zb,nb1,beta,z,nz1)
    else if ((ta.ne.'n').and.(tb.eq.'n')) then
       na1=size(za,1); na2=size(za,2)
       nb1=size(zb,1); nb2=size(zb,2)
       ! z(nz1,nz2)~alpha*j(za(na2,na1))*zb(nb1,nb2)
       if ((na2.gt.nz1).or.(na1.ne.nb1).or.(nb2.gt.nz2)) then
          write(*,*)
          write(*,'("Error(zgemm_wrap): inconsistent array dimensions")')
          write(*,'(" Z = alpha A B + beta Z")')
          write(*,'(" Z    : ",2i8)') nz1,nz2
          write(*,'(" j(A) : ",2i8)') na2,nb1
          write(*,'(" B    : ",2i8)') nb1,nb2
          write(*,*)
          stop
       end if
       ! call to BLAS routine
       call zgemm(ta,tb,na2,nb2,na1,alpha,za,na1,zb,nb1,beta,z,nz1)
    else if ((ta.ne.'n').and.(tb.ne.'n')) then
       na1=size(za,1); na2=size(za,2)
       nb1=size(zb,1); nb2=size(zb,2)
       ! z(nz1,nz2)~alpha*j(za(na2,na1)) * j(zb(nb2,nb1))
       if ((na2.gt.nz1).or.(na1.ne.nb2).or.(nb1.gt.nz2)) then
          write(*,*)
          write(*,'("Error(zgemm_wrap): inconsistent array dimensions")')
          write(*,'(" Z = alpha A j(B) + beta Z")')
          write(*,'(" Z    : ",2i8)') nz1,nz2
          write(*,'(" j(A) : ",2i8)') na2,nb1
          write(*,'(" j(B) : ",2i8)') nb2,nb1
          write(*,*)
          stop
       end if
       ! call to BLAS routine
       call zgemm(ta,tb,na2,nb1,na1,alpha,za,na1,zb,nb1,beta,z,nz1)
    end if
  end subroutine zgemm_wrap
end module blaswrappers
