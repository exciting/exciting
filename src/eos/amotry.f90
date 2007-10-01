real(8) function amotry(p,y,psum,nparam,ihi,fac)
implicit none
! arguments
integer, intent(in) :: ihi
integer, intent(in) :: nparam
real(8), intent(in) :: fac
real(8), intent(inout) :: p(nparam,nparam+1)
real(8), intent(inout) :: psum(nparam)
real(8), intent(inout) :: y(nparam+1)
! local variables
integer j
real(8) fac1,fac2,ytry
! external functions
real(8) fopt
external fopt
! automatic arrays
real(8) ptry(nparam)
fac1=(1.d0-fac)/dble(nparam)
fac2=fac1-fac
do j=1,nparam
  ptry(j)=psum(j)*fac1-p(j,ihi)*fac2
end do
ytry=fopt(ptry)
if (ytry.lt.y(ihi)) then
  y(ihi)=ytry
  do j=1,nparam
    psum(j)=psum(j)-p(j,ihi)+ptry(j)
    p(j,ihi)=ptry(j)
  end do
end if
amotry=ytry
return
end function
