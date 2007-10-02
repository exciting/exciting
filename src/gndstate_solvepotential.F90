subroutine gndstate_solvepotential(init,mu,nu,beta,f,n,dv)
  use modmain
  use modmpi
  integer,intent(in):: n
  logical,intent(inout)::init
  real(8),intent(inout)::nu(n)
  real(8),intent(inout)::mu(n),beta(n),f(n)
  real(8),intent(inout)::dv

  ! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(.false.,lradstp,ndmag,magmt,magir)
! convert the density from a coarse to a fine radial mesh
  call rfmtctof(rhomt)
  ! convert the magnetisation from a coarse to a fine radial mesh
  do idm=1,ndmag
     call rfmtctof(magmt(1,1,1,idm))
  end do
  ! add the core density to the total density
  call addrhocr
  ! calculate the charges
  call charge
  ! calculate the moments
  if (spinpol) call moment
  ! normalise the density
  call rhonorm
  ! compute the effective potential
  call poteff
  ! pack interstitial and muffin-tin effective potential and field into one array
  call packeff(.true.,n,nu)
  ! mix in the old potential and field with the new
  if(rank.eq.0) 	call mixer(init,beta0,betamax,n,nu,mu,beta,f,dv)
#ifdef MPI
  call  MPI_BCAST(nu(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  call  MPI_BCAST(mu(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  call  MPI_BCAST(f(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  call  MPI_BCAST(beta(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif
  ! unpack potential and field
  call packeff(.false.,n,nu)
  ! add the fixed spin moment effective field
  if (fixspin) call fsmfield
  ! Fourier transform effective potential to G-space
  call genveffig
end subroutine gndstate_solvepotential
