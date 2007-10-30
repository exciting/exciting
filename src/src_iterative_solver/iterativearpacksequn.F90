subroutine iterativearpacksecequn(ik,ispn,apwalm,vgpc,evalfv,evecfv)

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
  Complex(8)::                 zero, one
  parameter         (zero = (0.0D+0, 0.0D+0) ,one = (1.0D+0, 0.0D+0) )
  !ARPACK Interface vars

  integer:: ido, nev, ncv, lworkl, info, ierr, j,i
  integer:: nconv, maxitr, ishfts, mode, ldv
  integer::iparam(11), ipntr(14)
  complex(8)::workd(3*nmat(ik,ispn)),resid(nmat(ik,ispn)),v(nmat(ik,ispn),nstfv)
  complex(8)::d(nmat(ik,ispn))
  complex(8)::workl(3*9*nstfv*nstfv+5*3*nstfv)
  complex(8):: sigma, workev(9*nstfv)
  character         bmat*1, which*2
  real(8):: rwork(nmat(ik,ispn)), tol
  logical::rvec
  logical::          select(nmat(ik,ispn))
  !ZHPTR interface vars
  integer ::IPIV(nmat(ik,ispn))
  !parameters
  nev=nstfv
  ncv=3*nstfv
  n=npmat(ik,ispn)
  bmat  = 'G'
  which = 'LM'
  sigma = zero
  lworkl =3*9*nstfv*nstfv+5*3*nstfv 
  tol    = 0.0 
  ido    = 0
  info   = 0
  ishfts = 1
  maxitr = 300
  mode   = 2
  iparam(1) = ishfts
  iparam(3) = maxitr  
  iparam(7) = mode 
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevalfv(vkl(1,ik),evalfv)
  call hamiltonandoverlapsetup(npmat(ik,ispn),ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,h,o)
  !calculate Lu decomposition to be used in the reverse communication loop
  call zhptrf('U', N, o, IPIV, info )

  !################################################
  !# reverse comunication loop of arpack library: #
  !################################################
  do 
     call znaupd(ido,bmat,n,which,nstfv,tol,resid,ncv,v,nmat(ik,ispn),\
     iparam,ipntr,workd,workl,lworkl,rwork,info)
     if (ido .eq. -1 .or. ido .eq. 1) then

	call zhpmv("U",n,dcmplx(1.0,0.0),h,workd(ipntr(1)), 1, dcmplx(0,0),workd(ipntr(2)), 1)

        call zhptrs('U', N, 1, o, IPIV, workd(ipntr(2)), n, INFO )


        cycle
     else if (ido .eq. 2) then
 	call zhpmv("U",n,dcmplx(1.0,0.0),o,workd(ipntr(1)), 1, dcmplx(0,0),workd(ipntr(2)), 1)
        cycle
     else 
        exit
     endif
  end do
  if ( ierr .ne. 0 ) then
     print *, ' ' 
     print *, ' Error with znaupd, info = ', ierr
     print *, ' Check the documentation of znaupd'
     print *, ' '

  else
     rvec = .true.

     call zneupd  (rvec,'A',select,d,v,ldv,sigma,&
          workev,bmat,n,which,nev,tol,resid,ncv,v,&
          ldv, iparam, ipntr, workd, workl, lworkl, rwork,&
          ierr )

  endif
  evecfv(:,:,ispn)=v(:,:)
  evalfv(:,ispn)=d(:)
  call putevecfv(ik,evecfv)
  call putevalfv(ik,evalfv)

  return
end subroutine iterativearpacksecequn
!EOC
