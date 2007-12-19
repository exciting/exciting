subroutine  DIISseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain, only: nstfv,vkl,ngk,igkig,nmat,vgkl,timemat,npmat&
       ,apwordmax,lmmaxapw,natmtot,nkpt,nmatmax,nspnfv,timefv,ngkmax,zzero,zone
  use sclcontroll
  use diisinterfaces
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
  ! argumentstrialvec
  integer, 	intent(in) 		:: ik
  integer, 	intent(in) 		:: ispn
  real(8),    intent(in)    :: vgpc(3,ngkmax)
  complex(8), intent(in) 	:: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8), 	intent(inout) 	:: evalfv(nstfv,nspnfv)
  complex(8), intent(inout) :: evecfv(nmatmax,nstfv,nspnfv)

  ! local variables


  integer 	::is,ia,idiis,n,np,ievec,i
  real(8)  	::vl,vu,abstol
  real(8) 	::cpu0,cpu1
  real(8) 	::eps,rnorm
  complex(8):: hamilton(nmat(ik,ispn),nmat(ik,ispn))
  complex(8):: overlap(nmat(ik,ispn),nmat(ik,ispn))
  complex(8):: P(nmatmax,nmatmax)
  complex(8):: h(nmat(ik,ispn),nstfv,diismax) 
  complex(8):: s(nmat(ik,ispn),nstfv,diismax)
  complex(8):: r(nmat(ik,ispn),nstfv)
  complex(8):: trialvecs(nmat(ik,ispn),nstfv,diismax)
  real(8)::w(nmatmax),rnorms(nstfv)
  complex(8)::z
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
     call seceqfvprecond(n,hamilton,overlap,P,w,evalfv(:,ispn),evecfv(:,:,ispn))
     call writeprecond(ik,n,P,w)
  else
     iunconverged=nstfv	
     call readprecond(ik,n,P,w)    	
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call getevalfv(vkl(1,ik),evalfv)

     if( doprerotate_preconditioner()) then
        !      call prerotate_preconditioner(n,2*nstfv,hamilton,evecfv,P)
        !     call precondspectrumupdate(n,2*nstfv,hamilton,overlap,P,w)
     endif
     do idiis=1,diismax
        write(*,*)"diisiter", idiis
        !h(:,:,diis) holds matrix with current aproximate 
        !vectors multiplied with hamilton
        !o: same for overlap*evecfv
        call setuphsvect(n,iunconverged,hamilton,overlap,evecfv(:,:,ispn),nmatmax,&
             h(:,:,idiis),s(:,:,idiis))
        call rayleighqotient(n,iunconverged,evecfv(:,:,ispn)&
             , h(:,:,idiis),s(:,:,idiis),evalfv(:,ispn))
        write (777,*)w(:)
        write (778,*)evalfv(:,ispn)
        call residualvectors(n,iunconverged,h(:,:,idiis),s(:,:,idiis)&
             ,evalfv(:,ispn),r,rnorms)
        write(*,*)"rnorms",rnorms
        if  (allconverged(nstfv,rnorms)) exit	
        ! call remove_converged(evecmap(nstfv),iunconverged,r,h,s,trialvecs)
      
        call calcupdatevectors(n,iunconverged,P,w,r,evalfv,&
             evecfv(:,:,ispn),trialvecs(:,:,idiis)) 
        call setuphsvect(n,iunconverged,hamilton,overlap,trialvecs(:,:,ispn),n,&
             h(:,:,idiis),s(:,:,idiis))
            
        if(idiis.gt.1)then
           call diisupdate(idiis,iunconverged,n,h,s, trialvecs&
                ,evalfv(:,ispn),evecfv(:,:,ispn))
        endif
  		stop
end do

call cpu_time(cpu1)
endif
timefv=timefv+cpu1-cpu0
return
end subroutine DIISseceqnfv
!EOC
