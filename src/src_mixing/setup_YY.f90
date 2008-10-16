subroutine setup_YY(n,noldstepsmax,S,Y,YY)
implicit none
integer,intent(in)::n,noldstepsmax
real(8),intent(in)::S(n,noldstepsmax),Y(n,noldstepsmax)
real(8),intent(out)::YY(noldstepsmax,noldstepsmax)
YY=0
end subroutine
