!Driver for jakobi Davidson librarie from Gerard Sleijpen http://www.math.ruu.nl/people/sleijpen/
subroutine  jdseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)
use modmain, only: nstfv,vkl,ngk,igkig,nmat,vgkl,timemat,npmat&
       ,apwordmax,lmmaxapw,natmtot,nkpt,nmatmax,nspnfv,timefv,ngkmax,zzero,zone
  use sclcontroll
  use jacobidavidsoncommon
  implicit none
  ! argumentstrialvec
  integer, 	intent(in) 		:: ik
  integer, 	intent(in) 		:: ispn
  real(8),    intent(in)    :: vgpc(3,ngkmax)
  complex(8), intent(in) 	:: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8), 	intent(inout) 	:: evalfv(nstfv,nspnfv)
  complex(8), intent(inout) :: evecfv(nmatmax,nstfv,nspnfv)
  integer n,np
  !interface vars for JDQZ
  complex(8)::alpha(4*nstfv),beta(4*nstfv),eivec(nmat(ik,ispn),nstfv)
  real(8),parameter:: eps=1e-9,lock=1e-9
   integer::jmax,jmin,ldeg,lwork,i,mparam
  complex(8),allocatable::zwork(:,:)
   jmax=4*nstfv 
   jmin=2*nstfv
   mparam=30
   ldeg=2
  ! lwork=10 + 6*ldeg + 5*jmax + 3*nstfv
  lwork=4+mparam+ 5*jmax + 3*nstfv
  np=npmat(ik,ispn)
  n=nmat(ik,ispn)
  
allocate(zwork(n,lwork),hamilton(np),overlap(np))

 call hamiltonandoverlapsetup(np,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,hamilton,overlap)

 
if (iscl.gt.1) then
call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
do i=1,nstfv
call zcopy(n,evecfv(1,i,ispn),1,eivec(1,i),1)
end do
 !call getevalfv(vkl(1,ik),evalfv)
!else
!init evalfv
endif

call initprecon(n,lowesteval,hamilton,overlap)

call JDQZ( alpha, beta, eivec, .true., n, lowesteval, eps,& 
	nstfv,jmax , jmin, 1,mparam, ldeg, 100, 1000,&
	lock, -1, 1, zwork, lwork )

 call normalizep(n,nstfv,overlap,eivec,n)	
do i=1,nstfv
call zcopy(n,eivec(1,i),1,evecfv(1,i,ispn),1)
evalfv(i,ispn)=dble (alpha(i)/beta(i))
end do
call deallocprecon()
deallocate(zwork,hamilton,overlap)
end subroutine