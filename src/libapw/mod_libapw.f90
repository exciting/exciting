module mod_libapw

integer nrfmtmax
integer ordrfmtmax
integer, allocatable :: ordrfmt(:,:)
integer, allocatable :: nrfmt(:)
real(8), allocatable :: hmltrad(:,:,:,:)
real(8), allocatable :: beffrad(:,:,:,:,:)
real(8), allocatable :: ovlprad(:,:,:,:)
real(8), allocatable :: beffir(:,:)
real(8), allocatable :: densmt(:,:,:,:,:,:)
real(8), allocatable :: densir(:,:,:)
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: hrfmt(:,:,:)
integer, allocatable :: lrf(:,:)

contains

subroutine libapw_init 
use modmain
implicit none
integer is,l,io,ilo,ik,ikloc

#ifdef _LIBAPW_
if (allocated(nrfmt)) deallocate(nrfmt)
allocate(nrfmt(nspecies))
if (allocated(ordrfmt)) deallocate(ordrfmt)
allocate(ordrfmt(0:lmaxapw,nspecies))
nrfmt=0
ordrfmt=0
do is=1,nspecies
  do l=0,lmaxapw
    do io=1,apword(l,is)
      nrfmt(is)=nrfmt(is)+1
      ordrfmt(l,is)=ordrfmt(l,is)+1
    enddo
  enddo
  do ilo=1,nlorb(is)
    nrfmt(is)=nrfmt(is)+1
    ordrfmt(lorbl(ilo,is),is)=ordrfmt(lorbl(ilo,is),is)+1
  enddo
enddo
nrfmtmax=maxval(nrfmt)
ordrfmtmax=maxval(ordrfmt)
if (allocated(hmltrad)) deallocate(hmltrad)
allocate(hmltrad(lmmaxvr,nrfmtmax,nrfmtmax,natmtot))
hmltrad=0.d0
if (allocated(beffrad)) deallocate(beffrad)
allocate(beffrad(lmmaxvr,nrfmtmax,nrfmtmax,natmtot,ndmag))
hmltrad=0.d0
if (allocated(ovlprad)) deallocate(ovlprad)
allocate(ovlprad(0:lmaxapw,ordrfmtmax,ordrfmtmax,natmtot))
ovlprad=0.d0
if (allocated(beffir)) deallocate(beffir)
allocate(beffir(ngrtot,ndmag))
beffir=0.d0
if (allocated(densmt)) deallocate(densmt)
allocate(densmt(nrfmtmax,nrfmtmax,lmmaxvr,natmtot,nspinor,nspinor))
densmt=0.d0
if (allocated(densir)) deallocate(densir)
allocate(densir(ngrtot,nspinor,nspinor))
densir=0.d0
if (allocated(rfmt)) deallocate(rfmt)
allocate(rfmt(nrmtmax,nrfmtmax,natmcls))
if (allocated(hrfmt)) deallocate(hrfmt)
allocate(hrfmt(nrmtmax,nrfmtmax,natmcls))
if (allocated(lrf)) deallocate(lrf)
allocate(lrf(nrfmtmax,natmcls))

call lapw_load_global(natmtot,nspecies,lmaxvr,lmaxapw,apwordmax,nrmtmax,&
  &ngkmax,ngvec,ngrtot,nlomax,ias2is,intgv,ivg,ivgig,ngrid,igfft,cfunir, &
  &cfunig,gntyry,nstfv,nstsv,nmatmax,nrfmtmax,ordrfmtmax,evaltol,spinpol,&
  &ndmag,omega,natmcls,ic2ias,natoms_in_class)
do is=1,nspecies
  call lapw_load_species(is,nlorb(is),lorbl(1,is),apword(0,is),rmt(is),nrmt(is))
enddo
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  call lapw_load_kpoint(ngk(1,ik),igkig(1,1,ikloc),vgkc(1,1,1,ikloc),wkpt(ik))
enddo
call lapw_init
#endif
return
end subroutine

#ifdef _LIBAPW_
subroutine libapw_seceqn_init
use modmain
implicit none
integer ir,is,ia,ic,ias,l1,l2,io1,io2,ilo1,ilo2,i1,i2,nr,i
integer l1tmp(0:lmaxapw),l2tmp,lm
real(8) cb,t1
real(8), allocatable :: bmt(:,:,:)
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
!
! collect radial functions
do ic=1,natmcls
  ias=ic2ias(ic)
  is=ias2is(ias)
  i1=0
  do l1=0,lmaxapw
    do io1=1,apword(l1,is)
      i1=i1+1
      rfmt(:,i1,ic)=apwfr(:,1,io1,l1,ias)
      hrfmt(:,i1,ic)=apwfr(:,2,io1,l1,ias)
      lrf(i1,ic)=l1
    enddo
  enddo
  do ilo1=1,nlorb(is)
    i1=i1+1
    rfmt(:,i1,ic)=lofr(:,1,ilo1,ias)
    hrfmt(:,i1,ic)=lofr(:,2,ilo1,ias)
    lrf(i1,ic)=lorbl(ilo1,is)
  enddo
enddo
allocate(bmt(lmmaxvr,nrmtmax,ndmag))
cb=gfacte/(4.d0*solsc)
do ias=1,natmtot
  ic=ias2ic(ias)
  ia=ias2ia(ias)
  is=ias2is(ias)
  nr=nrmt(is)
  do ir=1,nr
    r2(ir)=spr(ir,is)**2
  end do
! hamiltonian radial integrals
  do i1=1,nrfmt(is)
    do i2=1,nrfmt(is)
      if (lrf(i1,ic).eq.lrf(i2,ic)) then
        do ir=1,nr
          fr(ir)=rfmt(ir,i1,ic)*hrfmt(ir,i2,ic)*r2(ir)
        enddo
        call fderiv(-1,nr,spr(:,is),fr,gr,cf)
        hmltrad(1,i1,i2,ias)=gr(nr)/y00
      else
        hmltrad(1,i1,i2,ias)=0.d0
      endif
      if (i1.ge.i2) then
        do lm=2,lmmaxvr
          do ir=1,nr
            fr(ir)=rfmt(ir,i1,ic)*rfmt(ir,i2,ic)*r2(ir)*veffmt(lm,ir,ias)
          enddo
          call fderiv(-1,nr,spr(:,is),fr,gr,cf)
          hmltrad(lm,i1,i2,ias)=gr(nr)
          hmltrad(lm,i2,i1,ias)=gr(nr)
        enddo
      endif
    enddo
  enddo 
! overlap integrals
  do l1=0,lmaxapw
    do io1=1,apword(l1,is)
      ovlprad(l1,io1,io1,ias)=1.d0
    enddo
  enddo
  l1tmp=0
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    l1tmp(l1)=l1tmp(l1)+1
    do io1=1,apword(l1,is)
      ovlprad(l1,io1,apword(l1,is)+l1tmp(l1),ias)=oalo(io1,ilo1,ias)
      ovlprad(l1,apword(l1,is)+l1tmp(l1),io1,ias)=oalo(io1,ilo1,ias)
    enddo
    l2tmp=0
    do ilo2=1,nlorb(is)
      if (l1.eq.lorbl(ilo2,is)) then
        l2tmp=l2tmp+1
        ovlprad(l1,apword(l1,is)+l1tmp(l1),apword(l1,is)+l2tmp,ias)=ololo(ilo1,ilo2,ias)
        ovlprad(l1,apword(l1,is)+l2tmp,apword(l1,is)+l1tmp(l1),ias)=ololo(ilo2,ilo1,ias)
      endif
    enddo
  enddo
! generate radial integals for magnetic field
  if (spinpol) then
    bmt=0.d0
    ! z
    bmt(:,:,1)=bxcmt(:,:,ias,ndmag)
    t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
    bmt(1,:,1)=bmt(1,:,1)+t1/y00
    if (ndmag.eq.3) then
      ! x
      bmt(:,:,2)=bxcmt(:,:,ias,1)
      t1=cb*(bfcmt(1,ia,is)+bfieldc(1))
      bmt(1,:,2)=bmt(1,:,2)+t1/y00
      ! y
      bmt(:,:,3)=bxcmt(:,:,ias,2)
      t1=cb*(bfcmt(2,ia,is)+bfieldc(2))
      bmt(1,:,3)=bmt(1,:,3)+t1/y00
    endif
    do i=1,ndmag
      do i1=1,nrfmt(is)
        do i2=1,nrfmt(is)
          if (i1.ge.i2) then
            do lm=1,lmmaxvr
              do ir=1,nr
                fr(ir)=rfmt(ir,i1,ic)*rfmt(ir,i2,ic)*r2(ir)*bmt(lm,ir,i)
              enddo
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              beffrad(lm,i1,i2,ias,i)=gr(nr)
              beffrad(lm,i2,i1,ias,i)=gr(nr)
            enddo !lm
          endif
        enddo
      enddo
    enddo !i
  endif ! spinpol
enddo !ias
deallocate(bmt)
if (spinpol) then
  ! z
  do ir=1,ngrtot
    beffir(ir,1)=bxcir(ir,ndmag)+cb*bfieldc(3)
  enddo
  if (ndmag.eq.3) then
    ! x
    do ir=1,ngrtot
      beffir(ir,2)=bxcir(ir,1)+cb*bfieldc(1)
    enddo
    ! y
    do ir=1,ngrtot
      beffir(ir,3)=bxcir(ir,2)+cb*bfieldc(2)
    enddo
  endif
endif
call lapw_seceqn_init(hmltrad,ovlprad,beffrad,apwfr,apwdfr,beffir,veffig) 
densmt=0.d0
densir=0.d0
return
end subroutine

subroutine libapw_rhomag
use modmain
implicit none
integer n,ias,ic,is,ir,ispn1,ispn2,lm3,j1,j2
real(8), allocatable :: fr(:,:,:)

n=nrfmtmax*nrfmtmax*lmmaxvr*natmtot*nspinor*nspinor
call mpi_grid_reduce(densmt(1,1,1,1,1,1),n,dims=(/dim_k/))
call mpi_grid_reduce(densir(1,1,1),ngrtot*nspinor*nspinor,dims=(/dim_k/))

do j1=1,nrfmtmax
  densmt(j1,j1,:,:,:,:)=0.5*densmt(j1,j1,:,:,:,:)
enddo
call timer_start(t_rho_mag_conv)
allocate(fr(nrmtmax,nspinor,nspinor))
do ias=1,natmtot
  ic=ias2ic(ias)
  is=ias2is(ias)
  do lm3=1,lmmaxvr
    fr=0.d0
    do ispn1=1,nspinor
      do ispn2=1,nspinor
        if ((ndmag.eq.1.and.(ispn1.eq.ispn2)).or.ndmag.ne.1) then
          do j2=1,nrfmtmax
            do j1=1,j2
              fr(:,ispn1,ispn2)=fr(:,ispn1,ispn2)+2*densmt(j1,j2,lm3,ias,ispn1,ispn2)*&
                &rfmt(:,j1,ic)*rfmt(:,j2,ic)
            enddo
          enddo
        endif
      enddo 
    enddo !l1
    if (spinpol) then
      rhomt(lm3,:,ias)=fr(:,1,1)+fr(:,2,2)
      magmt(lm3,:,ias,ndmag)=fr(:,1,1)-fr(:,2,2)
      if (ndmag.eq.3) then
        magmt(lm3,:,ias,1)=fr(:,1,2)
        magmt(lm3,:,ias,2)=fr(:,2,1)
      endif
    else
      rhomt(lm3,:,ias)=fr(:,1,1)
    endif
  enddo !lm3
enddo !ias
deallocate(fr)
if (spinpol) then
  rhoir(:)=densir(:,1,1)+densir(:,2,2)
  magir(:,ndmag)=densir(:,1,1)-densir(:,2,2)
  if (ndmag.eq.3) then
    magir(:,1)=densir(:,1,2)
    magir(:,2)=densir(:,2,1)
  endif
else
  rhoir(:)=densir(:,1,1)
endif
call timer_stop(t_rho_mag_conv)
call timer_start(t_rho_mag_sym)
! symmetrise the density
call symrf(1,rhomt,rhoir)
! symmetrise the magnetisation
if (spinpol) call symrvf(1,magmt,magir)
call timer_stop(t_rho_mag_sym)
! add the core density to the total density
call addrhocr
! calculate the charges
call charge
! calculate the moments
if (spinpol) call moment
! normalise the density
call rhonorm

return
end subroutine



#endif

end module
