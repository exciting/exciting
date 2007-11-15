subroutine iterativearpacksecequn(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain
  use modmpi
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   ispn   : first-variational spin index (in,integer)
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   vgpc   : G+k-vectors in Cartesian coordinates
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  ! This routine will perform several ARPACK iterations 

  !BOC
  implicit none
  ! arguments
  integer,intent(in)       :: ik
  integer,intent(in)       :: ispn
  real(8),intent(in)       :: vgpc(3,ngkmax)
  complex(8),intent(in)    :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8),intent(inout)    :: evalfv(nstfv,nspnfv)
  complex(8),intent(inout) :: evecfv(nmatmax,nstfv,nspnfv)

  ! local variables
  complex(8) 	:: h(npmat(ik,ispn))
  complex(8) 	:: o(npmat(ik,ispn))
  integer ::n
  real:: cpu0,cpu1
  Complex(8)::                 zero, one
  parameter         (zero = (0.0D+0, 0.0D+0) ,one = (1.0D+0, 0.0D+0) )
  !ARPACK Interface vars

  integer:: ido, nev, ncv, lworkl, info,infoznaupd, info2, j,i
  integer :: nevmax, ncvmax,nmax
  integer:: nconv, maxitr, ishfts, mode, ldv
  integer::iparam(11), ipntr(14)
  complex(8),allocatable::workd(:),resid(:),v(:,:),workev(:),workl(:),d(:)
  real(8),allocatable:: rwork(:),rd(:)
  integer,allocatable::idx(:)
  complex(8)::sigma
  character:: bmat*1, which*2
  real(8):: tol
  logical::rvec
  logical::          select(nmat(ik,ispn))
  !ZHPTR interface vars
  integer ::IPIV(nmat(ik,ispn))
  !parameters
  nev=nstfv +2
  ncv=2*nev+2
  nevmax=nev
  ncvmax= ncv
  nmax=nmatmax
  n=nmat(ik,ispn)
  ldv=n
  lworkl =3*ncvmax*ncvmax+5*ncvmax 
  allocate(workd(3*nmax))
  allocate(resid(nmax))
  allocate(v(ldv,ncvmax))
  allocate(workev(2*ncvmax))
  allocate(workl(lworkl))
  allocate(d(ncvmax))
  allocate(rwork(ncvmax))
  allocate(rd(ncvmax),idx(ncvmax))
  bmat  = 'G'
  which = 'LM'

  sigma=dcmplx(lowesteval,0)

  tol    = 1e-8
  ido    = 0
  info   = 0
  ishfts = 1
  maxitr = 2000
  mode   = 3
  iparam(1) = ishfts
  iparam(3) = maxitr  
  iparam(7) = mode 
  if(iscl.ne.1) then
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call getevalfv(vkl(1,ik),evalfv)
  endif
  call hamiltonandoverlapsetup(npmat(ik,ispn),ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,h,o)
  !calculate LU decomposition to be used in the reverse communication loop
#ifdef DEBUG
  write (333,*)"h",h,"o",o
#endif

  call cpu_time(cpu0)
  call zaxpy(npmat(ik,ispn),-sigma,o,1,h,1)
  call zhptrf('U', n, h, IPIV, info )
  if (info.ne.0)then
     write(*,*)"error in iterativearpacksecequn zhptrf ",info
     stop
  endif
  !################################################
  !# reverse comunication loop of arpack library: #
  !################################################
  do i=1,maxitr 
     call znaupd  &
          ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,  &
          ipntr, workd, workl, lworkl, rwork, infoznaupd)
#ifdef DEBUG
     ! write(*,*) "ido",ido
#endif

     if (ido .eq. -1 .or. ido .eq. 1) then

	call zhpmv("U",n,dcmplx(1.0,0.0),o,workd(ipntr(1)), 1,&
             dcmplx(0,0),workd(ipntr(2)), 1)
        call zhptrs('U', N, 1, h, IPIV, workd(ipntr(2)), n, INFO )
        if (info.ne.0)then
           write(*,*)"error in iterativearpacksecequn zhptrs ",info
           stop
        endif
     else if(ido .eq.1) then
        call zcopy (n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
        call zhptrs('U', N, 1, h, IPIV, workd(ipntr(2)), n, INFO )
     else if (ido .eq. 2) then
 	call zhpmv("U",n,dcmplx(1.0,0.0),o,workd(ipntr(1)), 1,&
             dcmplx(0,0),workd(ipntr(2)), 1)
     else 
        exit
     endif
  end do
  if ( infoznaupd .ne. 0 ) then
     print *, ' ' 
     print *, ' Error with znaupd, info = ', infoznaupd
     print *, ' Check the documentation of znaupd'
     print *, ' '
     stop
  else
     rvec = .true.

     call zneupd  (rvec,'A',select,d,v,n,sigma,&
          workev,bmat,n,which,nev,tol,resid,ncv,v,&
          n, iparam, ipntr, workd, workl, lworkl, rwork,&
          info2 )
     !         (rvec , howmny, select, d     ,
     !    &                   z    , ldz   , sigma , workev,
     !    &                   bmat , n     , which , nev   ,
     !    &                   tol  , resid , ncv   , v     ,
     !    &                   ldv  , iparam, ipntr , workd ,
     !    &                   workl, lworkl, rwork , info  )
     if ( info2 .ne. 0 ) then
        print *, ' ' 
        print *, ' Error with zneupd, info = ', info2
        print *, ' Check the documentation of zneupd'
        print *, ' '	 
	write(*,*)"eval",d(1:nev)
        write(*,*)"iter",i	
        stop
     endif

  endif
  call cpu_time(cpu1)
  timefv=timefv+cpu1-cpu0
#ifdef DEBUG
  write(*,*)"eval",d(1:nev)	
  write(*,*)"iterations",i
#endif
  if(rank.eq.0)write(60,*)"ARPACK iterations ", i
  rd=real(d)
  call sortidx (nevmax,rd(:),idx(:))
  do j=1,nstfv
     evecfv(:,j,ispn)=v(:,idx(j))
     evalfv(j,ispn)=rd(idx(j))
  end do
  evecfv(:,1:nstfv,ispn)=v(:,1:nstfv)
  evalfv(1:nstfv,ispn)=d(1:nstfv)
  call putevecfv(ik,evecfv)
  call putevalfv(ik,evalfv)
  deallocate(workd,resid,v,workev,workl,d)

  deallocate(rwork)
  return
end subroutine iterativearpacksecequn
!EOC
