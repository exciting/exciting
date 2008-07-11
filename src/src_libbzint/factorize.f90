!BOP
!
! !ROUTINE: factorize
!
! !INTERFACE:
subroutine factorize(n,x,k,div)

! !DESCRIPTION: 
!   This subroutine factorizes the real coordinates of a vector x,
!   the output is an integer vector k, such that k(i)/div=x(i)
!
! !REVISION HISTORY:
!   Created July 2008 by Sagmeister
!EOP
!BOC
  implicit none
  ! arguments
  integer(4), intent(in) :: n
  real(8), intent(in) :: x(n)
  integer(4), intent(out) :: div
  integer(4), intent(out) :: k(n)
  ! local variables
  integer, parameter :: maxint=10**4
  real(8), parameter :: eps=1.d-5
  real(8) :: dx
  do div=1,maxint
     k(:)=nint(dble(div)*x(:))
     dx=maxval(abs(dble(k)/dble(div)-x))
     if (dx.lt.eps) exit
  end do
  if (dx.ge.eps) then
     write(*,*)
     write(*,'("Error(factorize): factorization failed")')
     write(*,'(" maximum integer :",i12)') maxint
     write(*,'(" tolerance       :",g18.10)') eps
     write(*,'(" deviation       :",g18.10)') dx
     write(*,*)
     stop
  end if
  if (dx.gt.1.d-12) then
     write(*,*)
     write(*,'("Warning(factorize): small deviation in factorization")')
     write(*,'(" maximum deviation :",g18.10)') dx
     write(*,*)
  end if
end subroutine factorize
!EOC
