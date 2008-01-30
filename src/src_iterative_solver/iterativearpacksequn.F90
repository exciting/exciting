subroutine iterativearpacksecequn(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain
  use modmpi
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
  ! This routine will perform several ARPACK iterations 

  !BOC
  implicit none
#ifdef DEBUG
  !include declarations for timing output of ARPACK
#include "./debugf90.h"
#endif
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
  real:: cpu0,cpu1,cpu2
  Complex(8)::                 zero, one
  parameter         (zero = (0.0D+0, 0.0D+0) ,one = (1.0D+0, 0.0D+0) )
  !IO vars
  integer::koffset,recl
  character(256):: outfilenamestring,filetag
  external outfilenamestring

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
  logical:: select(nmat(ik,ispn))
  !ZHPTR interface vars
  integer ::IPIV(nmat(ik,ispn))
#ifdef DEBUG
  ndigit = -3
  logfil = 6
  mngets = 1
  mnaitr = 1
  mnapps = 1
  mnaupd = 1
  mnaup2 = 1
  mneigh = 1
  mneupd = 1
  open (logfil,file="ARPACK.OUT",action="WRITE") 
#endif
  !##################	
  !ARPACK parameters
  !##################
  nev=nstfv
  ncv=2*nev
  ncv=min(2*nev,maxncv)
  ncv=max(ncv,nev+2)
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
  resid(:)=0.0
  tol    = 1.e-8
  ido    = 0
  info   = 0
  ishfts = 1
  maxitr = 2000
  mode= 3
  iparam(1) = ishfts
  iparam(3) = maxitr  
  iparam(7) = mode 
  !################################
  !open file with previous residual
  !################################
  inquire(iolength=recl)resid  
  koffset=ik-firstk(procofk(ik))+1      
  infoznaupd=0
  
  !##################
  !setup hamiltonian#
  !##################
  call hamiltonandoverlapsetup(npmat(ik,ispn),ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,h,o)
  call cpu_time(cpu0)
  !#######################################################################
  !calculate LU decomposition to be used in the reverse communication loop
  !#######################################################################

  call zaxpy(npmat(ik,ispn),-sigma,o,1,h,1)
  !call zhptrf('U', n, h, IPIV, info )
 call ZPPTRF( 'U', n, h, INFO )
  if (info.ne.0)then
     write(*,*)"error in iterativearpacksecequn zhptrf ",info
     stop
  endif
  call cpu_time(cpu1)
  !################################################
  !# reverse comunication loop of arpack library: #
  !################################################
  do i=1,maxitr 
     call znaupd  &
          ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,  &
          ipntr, workd, workl, lworkl, rwork, infoznaupd)


     if (ido .eq. -1 .or. ido .eq. 1) then

	call zhpmv("U",n,dcmplx(1.0,0.0),o,workd(ipntr(1)), 1,&
             dcmplx(0,0),workd(ipntr(2)), 1)

        !call zhptrs('U', N, 1, h, IPIV, workd(ipntr(2)), n, INFO )
        call ZPPTRS( 'U', N, 1, h,workd(ipntr(2)) , n, INFO )
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
  !###############
  ! errorhandling
  !###############
  if ( infoznaupd .ne. 0 ) then
     print *, ' ' 
     print *, ' Error with znaupd, info = ', infoznaupd
     print *, ' Check the documentation of znaupd'
     print *, ' '
     stop
  else
    
 
     !########################
     !post processing of evec
     !########################
     rvec = .true.
     select=.true.
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
  call cpu_time(cpu2)
  timefv=timefv+cpu2-cpu0
#ifdef DEBUG
  close(logfil)
#endif
  if(rank.eq.0)write(60,*)"ARPACK iter", i
  if(rank.eq.0)write(60,*)"matrixsize",n&
       ,"time LU/iter",cpu1-cpu0,"/",cpu2-cpu1

  !##########################
  !sort and copy eigenvectors
  !##########################
  rd=real(d)
  call sortidx (nstfv,rd(:),idx(:))
  do j=1,nstfv
     evecfv(:,j,ispn)=v(:,idx(j))
     evalfv(j,ispn)=rd(idx(j))
  end do
  deallocate(workd,resid,v,workev,workl,d)
  deallocate(rwork,rd,idx)
  return
end subroutine iterativearpacksecequn
!EOC
