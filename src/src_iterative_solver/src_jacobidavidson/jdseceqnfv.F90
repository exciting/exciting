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
  complex(8)::alpha,beta,eivec(nmat(ik,ispn),nstfv)
  real(8),parameter:: eps=1e-9
   integer::jmax,jmin,ldeg,lwork
  complex(8),allocatable::zwork(:)
   jmax=4*nstfv
   jmin=2*nstfv
   ldeg=2
   lwork=10 + 6*ldeg + 5*jmax + 3*nstfv
  np=npmat(ik,ispn)
  n=nmat(ik,ispn)
  
allocate(hamilton(np),overlap(np))
 call hamiltonandoverlapsetup(np,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,hamilton,overlap)

 
!if (not firste) then
! call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
 !call getevalfv(vkl(1,ik),evalfv)
!else
!init evalfv
!endif

call initprecon()
allocate (zwork(n*lwork))
call JDQZ( alpha, beta, eivec, .true., n, lowesteval, eps,& 
	nstfv,jmax , jmin, 2, 30, ldeg, 100, 1000,&
	1.0e-9, 1, 1, zwork, lwork )
deallocate (zwork)

deallocate(hamilton,overlap)

end subroutine