!Driver for jakobi Davidson librarie from Gerard Sleijpen http://www.math.ruu.nl/people/sleijpen/
subroutine  jdseceqnfv(ik,ispn,apwalm,vgpc,evalfv,evecfv)
  use modmain, only: nstfv,vkl,ngk,igkig,nmat,vgkl,timemat,npmat&
       ,apwordmax,lmmaxapw,natmtot,nkpt,nmatmax,nspnfv,timefv,ngkmax,zzero,zone&
       ,scrpath,filext
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
  integer n,np,info
  real:: cpu0,cpu1
  integer::recl
  character(256):: outfilenamestring,krange,scrpathtmp
  !interface vars for JDQZ
  complex(8)::alpha(4*nstfv),beta(4*nstfv),eivec(nmat(ik,ispn),nstfv),target
  real(8),parameter:: eps=1e-8,lock=1e-8
  integer::jmax,jmin,ldeg,lwork,i,mparam
  complex(8),allocatable::zwork(:,:)
  jmax=3*nstfv 
  jmin=2*nstfv
  mparam=30
  ldeg=2 
  ! lwork=10 + 6*ldeg + 5*jmax + 3*nstfv
  lwork=4+mparam+ 5*jmax + 3*nstfv 
  np=npmat(ik,ispn)
  n=nmat(ik,ispn)

  allocate(zwork(n,lwork),hamilton(np),overlap(np))

  call hamiltonandoverlapsetup(np,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,hamilton,overlap)
  call cpu_time(cpu0)
  write(krange,*)ik
  write(krange,*)trim(krange),ispn
  scrpathtmp=scrpath
  outfilenamestring=trim(scrpathtmp)//"work"//trim(krange)//trim(filext)   
  inquire(iolength=recl)zwork
  if (iscl.gt.1) then
  info=-1
     open(80,file=outfilenamestring,action='READ', &
          form='UNFORMATTED',access='DIRECT',recl=recl)
     read(80,rec=1)zwork
     close(80)

     !add random
     else
     info=0
  endif
  target=dcmplx(lowesteval,0)
  call initprecon(n,target,hamilton,overlap)

  call JDQZ( alpha, beta, eivec, .true., n, target, eps,& 
       nstfv,jmax , jmin, 1,mparam, ldeg, 7, 1000,&
       lock, -1, 3, zwork, lwork ,info)

  open(80,file=outfilenamestring,action='WRITE', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  write(80,rec=1)zwork
  close(80)
  
  call normalizep(n,nstfv,overlap,eivec,n)	
  
  do i=1,nstfv
     call zcopy(n,eivec(1,i),1,evecfv(1,i,ispn),1)
     evalfv(i,ispn)=dble (alpha(i)/beta(i))
  end do
  
  call deallocprecon()
  deallocate(zwork,hamilton,overlap)
  call cpu_time(cpu1)

  timefv=timefv+cpu1-cpu0
end subroutine jdseceqnfv
