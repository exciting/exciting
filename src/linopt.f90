
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linopt
use modmain
implicit none
! local variables
integer ik,nsk(3),iw,jw
integer isym,lspl,iop
integer i,i1,i2,n,recl
real(8), parameter :: eps=1.d-8
real(8) wd(2),wplas,t1,t2
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: fw(:)
real(8), allocatable :: g(:)
real(8), allocatable :: cf(:,:)
real(8), allocatable :: e(:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: sc(:,:,:)
real(8), allocatable :: d(:)
real(8), allocatable :: eps1(:)
real(8), allocatable :: eps2(:)
real(8), allocatable :: sigma1(:)
real(8), allocatable :: sigma2(:)
real(8), allocatable :: delta(:,:)
real(8), allocatable :: pmatint(:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: pmat(:,:,:)
if ((usegdft).and.(xctype.lt.0)) then
  write(*,*)
  write(*,'("Error(linopt): generalised DFT cannot be used with exact &
   &exchange")')
  write(*,*)
  stop
end if
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
! allocate local arrays
allocate(w(nwdos))
allocate(fw(nwdos))
allocate(g(nwdos))
allocate(cf(3,nwdos))
n=nstsv*nstsv
allocate(e(n,nkpt))
allocate(f(n,nkpt))
allocate(sc(3,3,nsymcrys))
allocate(d(nsymcrys))
allocate(eps1(nwdos),eps2(nwdos))
allocate(sigma1(nwdos),sigma2(nwdos))
! allocate first-variational eigenvector array
allocate(evecfv(nmatmax,nstfv))
! allocate second-variational eigenvector array
allocate(evecsv(nstsv,nstsv))
! allocate the momentum matrix elements array
allocate(pmat(3,nstsv,nstsv))
! allocate intraband matrix elements
allocate(pmatint(nstsv,nkpt))
! set up for generalised DFT correction if required
allocate(delta(nstsv,nstsv))
if (usegdft) then
! initialisation required for generalised DFT
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  call readstate
  call poteff
  call linengy
  call genapwfr
  call genlofr
end if
! pre-calculation for symmetrisation
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  sc(:,:,isym)=dble(symlat(:,:,lspl))
  call r3mtm(sc(1,1,isym),binv,sc(1,1,isym))
  call r3mm(bvec,sc(1,1,isym),sc(1,1,isym))
  d(isym)=sc(1,1,isym)*sc(2,2,isym)-sc(1,2,isym)*sc(2,1,isym)
end do
! energy interval should start from zero
wdos(1)=0.d0
! generate energy grid
t1=(wdos(2)-wdos(1))/dble(nwdos)
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do
! find the record length for momentum matrix element file
inquire(iolength=recl) pmat
open(50,file='PMAT.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! loop over number of desired optical components
do iop=1,noptcomp
  i1=optcomp(1,iop)
  i2=optcomp(2,iop)
! open files for writting
  write(fname,'("EPSILON_",2I1,".OUT")') i1,i2
  open(60,file=trim(fname),action='WRITE',form='FORMATTED')
  write(fname,'("SIGMA_",2I1,".OUT")') i1,i2
  open(61,file=trim(fname),action='WRITE',form='FORMATTED')
  e(:,:)=0.d0
  f(:,:)=0.d0
  do ik=1,nkpt
! compute generalised DFT correction if required
    if (usegdft) then
      call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
      call getevecsv(vkl(1,ik),evecsv)
      call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
      call gdft(ik,apwalm,evecfv,evecsv,delta)
    end if
! read matrix elements from direct-access file
    read(50,rec=ik) pmat
    call linoptk(ik,i1,i2,sc,d,delta,pmat,e(1,ik),f(1,ik),pmatint(1,ik))
  end do
  t1=-4.d0*(pi**2)/omega
  f(:,:)=t1*f(:,:)
! calculate imaginary part of the interband dielectric function
  call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,n,n,e,f,eps2)
  do iw=1,nwdos
    t1=w(iw)
    if (abs(t1).gt.eps) then
      eps2(iw)=eps2(iw)/(t1**2)
    else
      eps2(iw)=0.d0
    end if
  end do
! calculate the intraband Drude-like contribution and plasma frequency
  if (intraband) then
    write(fname,'("PLASMA_",2I1,".OUT")') i1,i2
    open(62,file=trim(fname),action='WRITE',form='FORMATTED')
    wd(1)=efermi-0.001d0
    wd(2)=efermi+0.001d0
    call brzint(nsmdos,ngridk,nsk,ikmap,3,wd,nstsv,nstsv,evalsv,pmatint,g)
    wplas=sqrt(g(2)*8.d0*pi/omega)
    do iw=1,nwdos
      if (abs(w(iw)).gt.eps) then
! add the intraband contribution to the imaginary part of the tensor
        eps2(iw)=eps2(iw)+swidth*(wplas**2)/(w(iw)*(w(iw)**2+swidth**2))
      end if
    end do
! write plasma frequency to file
    write(62,'(G18.10," : plasma frequency")') wplas
    close(62)
  end if
! calculate real part of the dielectric function
  if (i1.eq.i2) then
    t1=1.d0
  else
    t1=0.d0
  end if
! Kramers-Kronig transform to find real part of dielectric tensor
  do iw=1,nwdos
    do jw=1,nwdos
      t2=w(jw)**2-w(iw)**2
      if (abs(t2).gt.eps) then
        fw(jw)=w(jw)*eps2(jw)/t2
      else
        fw(jw)=0.d0
      end if
    end do
    call fderiv(-1,nwdos,w,fw,g,cf)
    eps1(iw)=t1+(2.d0/pi)*g(nwdos)
  end do
! write dielectric function to a file
  do iw=1,nwdos
    write(60,'(2G18.10)') w(iw),eps1(iw)
  end do
  write(60,'("     ")')
  do iw=1,nwdos
    write(60,'(2G18.10)') w(iw),eps2(iw)
  end do
! calculate optical conductivity
  sigma1(:)=(eps2(:))*w(:)/(4.d0*pi)
  sigma2(:)=-(eps1(:)-t1)*w(:)/(4.d0*pi)
! write the interband optical conductivity to a file
  do iw=1,nwdos
    write(61,'(2G18.10)') w(iw),sigma1(iw)
  end do
  write(61,'("     ")')
  do iw=1,nwdos
    write(61,'(2G18.10)') w(iw),sigma2(iw)
  end do
  close(60)
  close(61)
! end loop over number of components
end do
close(50)
write(*,*)
write(*,'("Info(linopt):")')
write(*,'(" dielectric tensor written to EPSILON_ij.OUT")')
write(*,*)
write(*,'(" optical conductivity written to SIGMA_ij.OUT")')
if (intraband) then
  write(*,*)
  write(*,'(" plasma frequency written to PLASMA_ij.OUT")')
end if
write(*,*)
write(*,'(" for all requested components i,j")')
deallocate(w,fw,g,cf,e,f,sc,d)
deallocate(eps1,eps2)
deallocate(sigma1,sigma2)
deallocate(evecfv,evecsv,pmat,pmatint)
if (usegdft) deallocate(delta,apwalm)
return
end subroutine

