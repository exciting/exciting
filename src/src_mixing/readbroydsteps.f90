subroutine readbroydsteps_and_init_SY(noldsteps,noldstepsmax,n,S,Y)
implicit none
integer, intent(in)::noldsteps,noldstepsmax,n
real(8),INTENT(OUT)::S(n,noldstepsmax),Y(n,noldstepsmax)
integer i
open(23,file="BROYDEN.OUT",ACTION= "READ",FORM='UNFORMATTED')

Do i=1,noldsteps
read(23)s(:,i),y(:,i)
end do
close(23)

!setup Y ...


end subroutine
