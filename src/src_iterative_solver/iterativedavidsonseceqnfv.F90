subroutine  iterativedavidsonseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain
  use sclcontroll
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   ispn   : first-variational spin index (in,integer)
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   vgpc   : G+k-vectors in Cartesian coordinates
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  ! This routine will perform several Bock Davidson iterations following the sceme:
  ! 1. for each of m bands:
  !   a. calculate Residual
  !$$
  !\ket{\mathbf{R}\left(\ket{\mathbf{A}^{ap}},E^{ap}\right)}=(\mathbf{H}-E^{ap}\mathbf{S})\ket{ \mathbf{A}^{ap}}
  !$$
  !   b. calculate $\delta \mathbf{A}$
  ! 2. solve Projected system in evecsv+$\delta \mathbf{A}$ subspace
  !EOP
  !BOC
  implicit none
  ! arguments
  integer, 	intent(in) 		:: ik
  integer, 	intent(in) 		:: ispn
  real(8),    intent(in)      :: vgpc(3,ngkmax)
  complex(8), intent(in) 		:: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8), 	intent(inout) 	:: evalfv(nstfv,nspnfv)
  complex(8), intent(inout) 	:: evecfv(nmatmax,nstfv,nspnfv)

  ! local variables


  integer 	::is,ia,i,n,np,ievec,iv
  complex(8)	::evecp(2*nstfv,nstfv),r(nmat(ik,ispn))
  real(8)  	::vl,vu,abstol,evalp(nstfv)
  real(8) 	::cpu0,cpu1
  real(8) 	::eps,rnorm
  complex(8) 	:: h(npmat(ik,ispn)),hprojected(nstfv*2*(nstfv*2+1)/2)
  complex(8) 	:: o(npmat(ik,ispn)),oprojected(nstfv*2*(nstfv*2+1)/2)
  complex(8) 	:: hminuses(npmat(ik,ispn)),da(nmat(ik,ispn),nstfv)
  complex(8)X(nmatmax,nmatmax)
  real(8)::w(nmatmax)

  if ((ik.lt.1).or.(ik.gt.nkpt)) then
     write(*,*)
     write(*,'("Error(seceqnfv): k-point out of range : ",I8)') ik
     write(*,*)
     stop
  end if
  n=nmat(ik,ispn)
  np=npmat(ik,ispn)

  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!
  call cpu_time(cpu0)
  call hamiltonandoverlapsetup(np,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,h,o)
  call cpu_time(cpu1)
#ifdef DEBUG
  write(112,*)"h",h
  write(113,*)"o",o
#endif
  !$OMP CRITICAL
  timemat=timemat+cpu1-cpu0
  !$OMP END CRITICAL
  !update eigenvectors with iteration
  call cpu_time(cpu0)
  if(calculate_preconditioner()) then
     call seceqfvprecond(ik,n,h,o,X,evalfv,evecfv)
  else
     call readprecond(ik,n,X,w)
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call getevalfv(vkl(1,ik),evalfv)

     do i=1,3
        do ievec=1,nstfv
           call residualvector(n,np,h,o,evecfv(:,ievec,ispn),r(:),rnorm)
           if  (rnorm.lt.reps) exit
           call calcupdatevector(n,X,r(:),HminuseS(:),da(:,ievec) ) 
           call diisupdate(n,da,o,evecfv)
        end do

        do ievec=1,nstfv
           !calculate new eigenvalues from rayreigh quotient
           call rayleighqotient(n,evecfv(:,ievec,ispn)&
                ,h,o,evalfv(ievec,ispn))
        end do
     end do
	 call prerotate_preconditioner(n,h,o,evecfv,X)
     call cpu_time(cpu1)
  endif
  timefv=timefv+cpu1-cpu0
  return
end subroutine iterativedavidsonseceqnfv
!EOC
