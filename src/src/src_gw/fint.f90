!BOP
! !ROUTINE: fint
! !INTERFACE:
subroutine fint(np,n,x,f,g)
! !INPUT/OUTPUT PARAMETERS:
!   np : order of fitting polynomial (in,integer)
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   f  : function array (in,real(n))
!   g  : integrated function (out,real(n))
! !DESCRIPTION:
!   Calculates the integrals $g(x_i)$ of a function $f$ defined on a set of
!   points $x_i$ as follows
!   $$ g(x_i)=\int_0^{x_i} f(x)\,dx. $$
!   This is performed by piecewise fitting the function to polynomials of order
!   $n_p-1$ and performing the integrations analytically.
!
! !REVISION HISTORY:
!   Created May 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: np
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(in) :: f(n)
real(8), intent(out) :: g(n)
! local variables
integer i,i0,npo2
! automatic arrays
real(8) c(np)
! external functions
real(8) polynom
external polynom
if (n.lt.2) then
  write(*,*)
  write(*,'("Error(fint): n < 2 : ",I8)') n
  write(*,*)
  stop
end if
if (np.lt.2) then
  write(*,*)
  write(*,'("Error(fint): np < 2 : ",I8)') np
  write(*,*)
  stop
end if
if (n.lt.np) then
  write(*,*)
  write(*,'("Error(fint): n < np : ",2I8)') n,np
  write(*,*)
  stop
end if
npo2=np/2
g(1)=0.d0
do i=2,n
  if (i.le.npo2) then
    i0=1
  else if (i.gt.n-npo2) then
    i0=n-np+1
  else
    i0=i-npo2
  end if
  g(i)=polynom(-1,np,x(i0),f(i0),c,x(i))+g(i0)
end do
return
end subroutine
!EOC
