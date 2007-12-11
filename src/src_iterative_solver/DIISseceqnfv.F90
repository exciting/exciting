subroutine  DIISseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain, only: nstfv,vkl,ngk,igkig,nmat,vgkl,timemat,npmat&
       ,apwordmax,lmmaxapw,natmtot,nkpt,nmatmax,nspnfv,timefv,ngkmax
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
  real(8),    intent(in)    :: vgpc(3,ngkmax)
  complex(8), intent(in) 	:: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8), 	intent(inout) 	:: evalfv(nstfv,nspnfv)
  complex(8), intent(inout) :: evecfv(nmatmax,nstfv,nspnfv)

  ! local variables


  integer 	::is,ia,idiis,n,np,ievec,iv
  real(8)  	::vl,vu,abstol
  real(8) 	::cpu0,cpu1
  real(8) 	::eps,rnorm
  complex(8) 	:: hamilton(nmat(ik,ispn),nmat(ik,ispn)),hprojected(nstfv*2*(nstfv*2+1)/2)
  complex(8) 	:: overlap(nmat(ik,ispn),nmat(ik,ispn)),oprojected(nstfv*2*(nstfv*2+1)/2)
  complex(8)::P(nmatmax,nmatmax), h(nmat(ik,ispn),nstfv,diismax) ,&
       s(nmat(ik,ispn),nstfv,diismax),&
       r(nmat(ik,ispn),nstfv),subspacevectors(nmat(ik,ispn),nstfv,diismax)
  real(8)::w(nmatmax),rnorms(nstfv)
  integer evecmap(nstfv),  iunconverged
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
  call hamiltonandoverlapsetupnotpacked(n,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,hamilton,overlap)
  
  call cpu_time(cpu1)

  !$OMP CRITICAL
  timemat=timemat+cpu1-cpu0
  !$OMP END CRITICAL
  !update eigenvectors with iteration
  call cpu_time(cpu0)
  if(calculate_preconditioner()) then
     call seceqfvprecond(ik,n,hamilton,overlap,P,evalfv,evecfv)
  else
     iunconverged=nstfv	
     call readprecond(ik,n,P,w)    	
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call getevalfv(vkl(1,ik),evalfv)
     call prerotate_preconditioner(n,hamilton,overlap,P,w)
     do idiis=1,diismax
        call setuphsvect(n,hamilton,overlap,evecfv,h,s)
        call rayleighqotient(n,iunconverged,evecfv(:,:,ispn)&
             ,hamilton,overlap,evalfv(:,ispn))
        call residualvectors(n,iunconverged,h(:,:,idiis),s(:,:,idiis),evalfv(:,ispn),r,rnorms)
        if  (allconverged(nstfv,rnorms)) exit	
        call remove_converged(evecmap(nstfv),iunconverged,r,h,s,subspacevectors)

        call calcupdatevectors(n,idiis,iunconverged,P,h,s,subspacevectors) 
        call diisupdate(idiis,iunconverged,n,h(:,:,idiis),s(:,:,idiis), subspacevectors(:,:,idiis),evecfv)
        do ievec=1,nstfv
           !calculate new eigenvalues from rayreigh quotient
           call rayleighqotient(n,iunconverged,evecfv(:,ievec,ispn)&
                ,hamilton,overlap,evalfv(ievec,ispn))
        end do
     end do

     call cpu_time(cpu1)
  endif
  timefv=timefv+cpu1-cpu0
  return
end subroutine DIISseceqnfv
!EOC
