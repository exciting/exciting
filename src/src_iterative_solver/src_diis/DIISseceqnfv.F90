subroutine  DIISseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain, only: nstfv,vkl,ngk,igkig,nmat,vgkl,timemat,npmat&
       ,apwordmax,lmmaxapw,natmtot,nkpt,nmatmax,nspnfv,timefv,ngkmax,zzero,zone
  use sclcontroll
  use diisinterfaces
  use modfvsystem 
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

  type(evsystem)::system
  logical::packed,jacdav
  integer 	::is,ia,idiis,n,np,ievec,i,info,flag
  real(8)  	::vl,vu,abstol
  real(8) 	::cpu0,cpu1
  real(8) 	::eps,rnorm
  complex(8),allocatable:: P(:,:)
  complex(8),allocatable::         h(:,:,:) 
  complex(8),allocatable::         s(:,:,:)
  complex(8),allocatable::         r(:,:)
  complex(8),allocatable:: trialvecs(:,:,:)
  complex(8),allocatable:: eigenvector(:,:)
  real(8),allocatable:: eigenvalue(:)

  real(8)::w(nmatmax),rnorms(nstfv)
  complex(8)::z
  integer iunconverged,evecmap(nstfv)

  if ((ik.lt.1).or.(ik.gt.nkpt)) then
     write(*,*)
     write(*,'("Error(seceqnfv): k-point out of range : ",I8)') ik
     write(*,*)
     stop
  end if
  n=nmat(ik,ispn)
  np=npmat(ik,ispn)
  allocate( P(nmatmax,nmatmax))
  allocate(h(nmat(ik,ispn),nstfv,diismax)) 
  allocate( s(nmat(ik,ispn),nstfv,diismax))
  allocate(    r(nmat(ik,ispn),nstfv))
  allocate(  trialvecs(nmat(ik,ispn),nstfv,diismax))
  allocate(eigenvector(nmat(ik,ispn),nstfv))
  allocate(eigenvalue(nstfv))

  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!
  call cpu_time(cpu0)
  packed=.false.
jacdav=.false.
  call newsystem(system,packed,n)
  call hamiltonandoverlapsetup(system,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc)

  call cpu_time(cpu1)

  !$OMP CRITICAL
  timemat=timemat+cpu1-cpu0
  !$OMP END CRITICAL
  !update eigenvectors with iteration
  recalculate_preconditioner=.false.
  call cpu_time(cpu0)
  if(calculate_preconditioner()) then
     P=0
     w=0
     call seceqfvprecond(n,system,P,w,evalfv(:,ispn),evecfv(:,:,ispn))
     call writeprecond(ik,n,P,w)
  else

     iunconverged=nstfv	
     call readprecond(ik,n,P,w)    	
     !    write(*,*)"readeigenvalues",w
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call getevalfv(vkl(1,ik),evalfv)

      call zlarnv(2, iseed, n*nstfv, eigenvector)
      eigenvector=cmplx(dble(eigenvector),0.)
      call zscal(n*nstfv,dcmplx(1e-3/n/nstfv,0.),eigenvector,1)
     do i=1,nstfv
        call zcopy(n ,evecfv(1,i,ispn),1,eigenvector(1,i),1)
   !     call zaxpy(n ,zone,evecfv(1,i,ispn),1,eigenvector(1,i),1)
        eigenvalue(i)=evalfv(i,ispn)
        evecmap(i)=i
     end do

     if( doprerotate_preconditioner()) then

        !write(777,*)P
        !call prerotate_preconditioner(n,2*nstfv,hamilton,P)
        !call normalize(n,2*nstfv,overlap,P,nmatmax)	
        !call precondspectrumupdate(n,2*nstfv,hamilton,overlap,P,w)
        !write(778,*)P     
        !stop 
     endif
    	 if(jacdav)   call jacdavblock(n, iunconverged, system, n, & 
          eigenvector, h(:,:,idiis), s(:,:,idiis), eigenvalue, &
          trialvecs(:,:,idiis), h(:,:,idiis), 0)
     do idiis=1,diismax
        write(*,*)"diisiter", idiis
        !h(:,:,diis) holds matrix with current aproximate 
        !vectors multiplied with hamilton
        !o: same for overlap*evecfv
        call setuphsvect(n,iunconverged,system,eigenvector,n,&
             h(:,:,idiis),s(:,:,idiis))
        call rayleighqotient(n,iunconverged,eigenvector&
             , h(:,:,idiis),s(:,:,idiis),eigenvalue)
        call residualvectors(n,iunconverged,h(:,:,idiis),s(:,:,idiis)&
             ,eigenvalue,r,rnorms)
             


        do i=1,nstfv
           if(evecmap(i).ne.0) then
              call zcopy (n,eigenvector(1,evecmap(i)), 1,evecfv(1,i,ispn),1)
          
              evalfv(i,ispn)=eigenvalue(evecmap(i))
               write(203,*)idiis,ik,i,rnorms(evecmap(i))
           endif
        end do
       
          
        if  (allconverged(iunconverged,rnorms).or. idiis.eq.(diismax-1)) exit	
        call remove_converged(evecmap,iunconverged,&
             rnorms,n,r,h,s,eigenvector,eigenvalue,trialvecs)
        if (rnorms(idamax(iunconverged,rnorms,1)).gt.1e-1.and.(idiis.gt.1)) then
           recalculate_preconditioner=.true.
           write(*,*)"recalculate preconditioner"
           exit
        endif
  
      write(*,*)"norm,beforecalc" , dznrm2( n, eigenvector, 1)
        if(.not.jacdav)then
           call calcupdatevectors(n,iunconverged,P,w,r,eigenvalue,&
                eigenvector,trialvecs(:,:,idiis))  
   write(*,*)"norm,afterecalc" , dznrm2( n, eigenvector, 1)
        else
           call jacdavblock(n, iunconverged, system, n, & 
                eigenvector, h(:,:,idiis), s(:,:,idiis), eigenvalue, &
                trialvecs(:,:,idiis), h(:,:,idiis), 1) 
           call zaxpy(n*iunconverged,zone,trialvecs(1,1,idiis),1,eigenvector(1,1),1)
             call zcopy(n*iunconverged,trialvecs(1,1,idiis),1,eigenvector(1,1),1)
        endif
        
        call setuphsvect(n,iunconverged,system,eigenvector,n,&
             h(:,:,idiis),s(:,:,idiis)) 
        if(idiis.gt.1)then
           call diisupdate(idiis,iunconverged,n,h,s, trialvecs&
                ,eigenvalue,eigenvector,info)
                 write(*,*)"norm,afterdiis" , dznrm2( n, eigenvector, 1)
                 if(dznrm2( n, eigenvector, 1).lt.0.1e-2) stop
                  
        endif
     end do
   call normalize(n,nstfv,system%overlap%za,evecfv(:,:,ispn),nmatmax)	
     if ( recalculate_preconditioner .or. (idiis .gt. diismax-1)) then 
        call seceqfvprecond(n,system,P,w,evalfv(:,ispn),evecfv(:,:,ispn))
        call writeprecond(ik,n,P,w)
        write(*,*)"recalculate preconditioner"
     endif
     call cpu_time(cpu1)
	 if(jacdav)     call jacdavblock(n, iunconverged, system, n, & 
          eigenvector, h(:,:,idiis), s(:,:,idiis), eigenvalue, &
          trialvecs(:,:,idiis), h(:,:,idiis), -1) 

  endif

  deallocate(P,h,s, r,trialvecs,eigenvector, eigenvalue)



  call deleteystem(system)
  timefv=timefv+cpu1-cpu0

  return
end subroutine DIISseceqnfv
!EOC
