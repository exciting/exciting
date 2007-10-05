
! Copyright (C) 2006 S. Sagmeister.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_pade
  implicit none
contains

  subroutine pade(m,z,n,iw,ih,h)
    !
    ! Implementation of a Pade approximant using Thiele's method.
    ! Expressions taken from K. Lee, Phys. Rev. B 54 R8286 (1996)
    !
    ! Created Dec. 2006 (Sagmeister)
    !
    use m_ctdfrac
    implicit none
    ! arguments
    integer, intent(in) :: m
    complex(8), intent(in) :: z(m)
    integer, intent(in) :: n
    complex(8), intent(in) :: iw(n), ih(n)
    complex(8), intent(out) :: h(m)
    ! local variables
    character(*), parameter :: thisnam = 'pade'
    complex(8), allocatable :: acoef(:), bcoef(:), a(:), aa(:), bb(:), c(:,:)
    complex(8) :: zz
    integer :: j,l,k

    ! require order higher than two
    if (n < 2) then
       write(*,*) 'Error('//thisnam//'): approximant order too small (< 2)'
       call terminate
    end if

    ! allocate
    allocate(acoef(n),bcoef(0:n),a(n),aa(0:n),bb(0:n),c(n,n))

    ! coefficients for numerator and denominator
    a(1) = ih(1)
    c(1,:)=ih(:)
    do j=2,n
       do l=2,n
          c(j,l)=(a(j-1)-c(j-1,l))/((iw(l)-iw(j-1))*c(j-1,l))
       end do
       a(j)=c(j,j)
    end do

    ! calculate values at given frequencies
    do k=1,m
       zz = z(k)
       bcoef(:) = (1.d0,0.d0)
       bcoef(0) = (0.d0,0.d0)
       acoef(1) = a(1)
       do j=2,n
          acoef(j) = a(j)*(zz-iw(j-1))
       end do
       ! continued fraction evaluation of the Pade approximant
       call ctdfrac(n,acoef,bcoef,h(k))
    end do

    deallocate(acoef,bcoef,a,aa,bb,c)

  end subroutine

end module

