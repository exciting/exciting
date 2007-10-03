subroutine gndstate_solvepotential(init,mu,nu,beta,f,n,dv)
  use modmain
  use modmpi
  integer,intent(in):: n
  logical,intent(inout)::init
  real(8),intent(inout)::nu(n)
  real(8),intent(inout)::mu(n),beta(n),f(n)
  real(8),intent(inout)::dv

end subroutine gndstate_solvepotential
