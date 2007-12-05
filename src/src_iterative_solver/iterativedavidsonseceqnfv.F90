subroutine  iterativedavidsonseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)

 !USES:
use modmain
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
real(8) 	::residualeps=1e-6
logical 	::blockdavidsonconverged
real(8) 	::eps,rnorm
complex(8) 	:: h(npmat(ik,ispn)),hprojected(nstfv*2*(nstfv*2+1)/2)
complex(8) 	:: o(npmat(ik,ispn)),oprojected(nstfv*2*(nstfv*2+1)/2)
complex(8) 	:: hminuses(npmat(ik,ispn)),da(nmat(ik,ispn),nstfv)


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

call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
call getevalfv(vkl(1,ik),evalfv)


#ifdef DEBUG
write(114,*)"evecfv" ,evecfv
#endif

do i=1,3
	do ievec=1,nstfv	
	!do ievec=1,1	
	
   	! blas call means : HminuseS(:)=h(:)-evalfv(ievec,ispn)*o(:)
   		call zcopy(np,h,1,hminuses,1)
   		call zaxpy(np,dcmplx(-evalfv(ievec,ispn),0),o,1,hminuses,1)
#ifdef DEBUG
write(115,*)"hminuses",hminuses
#endif
      	call residualvector(n,np,hminuses(:),evecfv(:,ievec,ispn),nmatmax,r(:),rnorm)
      	if (rnorm.lt.residualeps)then
      		 blockdavidsonconverged=.true.
      	end if

      	call calcupdatevector(n,r(:),HminuseS(:),da(:,ievec) ) 
	end do
#ifdef DEBUG
	write(331,*)"r",r
	write(332,*)"HminuseS",HminuseS
	write(333,*)"da",da
#endif
	call setupprojectedhamilton(n,nstfv,h,o,nmatmax,evecfv(:,:,ispn),evalfv(:,ispn),da(:,:),hprojected(:),oprojected(:))
#ifdef DEBUG
	write(334,*)"hprojected",hprojected
	write(335,*)"oprojected",oprojected
#endif
	call projectedsecequn(nstfv,hprojected(:),oprojected(:),evecp(:,:),evalp(:))
#ifdef DEBUG
	write(336,*)"evalp",evalp
	write(337,*)"evalfv",evalfv(:,ispn)
#endif
write(338,*)"evecfv",evecfv(:,:,ispn)
	call updateevecfv(n,nstfv,da(:,:),nmatmax,evecfv(:,:,ispn),evalfv(:,ispn),evecp(:,:),evalp(:))
do ievec=1,nstfv
!calculate new eigenvalues from rayreigh quotient
call rayleighqotient(n,evecfv(:,ievec,ispn),h,o,evalfv(ievec,ispn))
end do

end do

call cpu_time(cpu1)
!$OMP CRITICAL
timefv=timefv+cpu1-cpu0
!$OMP END CRITICAL 

call putevecfv(ik,evecfv)
call putevalfv(ik,evalfv)

return
end subroutine
!EOC
