subroutine appent_current_to_broyden_file(n,potential,residual)
implicit none
integer ,intent(in)::n
real(8),intent(in)::potential(n),residual(n)
open(23,file="BROYDEN.OUT",ACCESS="APPEND",FORM='UNFORMATTED')
write(23)potential,residual
close(23)
end subroutine
