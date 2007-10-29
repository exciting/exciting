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
integer, 	intent(in) 		:: ik
integer, 	intent(in) 		:: ispn
real(8),    intent(in)      :: vgpc(3,ngkmax)
complex(8), intent(in) 		:: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), 	intent(inout) 	:: evalfv(nstfv,nspnfv)
complex(8), intent(inout) 	:: evecfv(nmatmax,nstfv,nspnfv)

! local variables
complex(8) 	:: h(npmat(ik,ispn))
complex(8) 	:: o(npmat(ik,ispn))
!ARPACK Interface vars

integer:: ido, n, nev, ncv, lworkl, info, ierr, j,nconv, maxitr, ishfts, mode
integer::iparam(11), ipntr(14), ipiv(nmat(ik,ispn))
complex(8)::workd(3*nmat(ik,ispn)),resid(nmat(ik,ispn)),v(nmat(ik,ispn),nstfv)
 complex(8)::workl(3*9*nstfv*nstfv+5*3*nstfv)
 character         bmat*1, which*2
 real(8):: rwork(nmat(ik,ispn)), tol
call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
call getevalfv(vkl(1,ik),evalfv)
call hamiltonandoverlapsetup(npmat(ik,ispn),ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc,h,o)
do
         call znaupd( ido, bmat, n, which, nev, tol, resid, ncv,\
v, nmat(ik,ispn), iparam, ipntr, workd, workl, lworkl,\
rwork, info )
     if (ido .eq. -1 .or. ido .eq. 1) then
     !   call matvecA (n, workd(ipntr(1)), temp_array)
     !   call solveM  (n, temp_array, workd(ipntr(2)))
        cycle
    else if (ido .eq. 2) then
      !  call matvecM (n, workd(ipntr(1)), workd(ipntr(2)))
        cycle
    else 
    exit
    endif
end do
  
call putevecfv(ik,evecfv)
call putevalfv(ik,evalfv)

return
end subroutine
!EOC
