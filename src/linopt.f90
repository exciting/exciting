
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linopt
use modmain
implicit none
! local variables
integer ik,ist1,ist2
integer isym,recl,d,i,m,n
integer i1,i2,iw,jw,nsk(3)
real(8), parameter :: eps=1.d-8
real(8) e12,t1,t2
real(8) sum,wint(2)
complex(8) zt1(3)
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: fw(:)
real(8), allocatable :: g(:)
real(8), allocatable :: cf(:,:)
real(8), allocatable :: e(:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: a(:,:,:)
real(8), allocatable :: s(:)
real(8), allocatable :: eps1(:,:)
real(8), allocatable :: eps2(:,:)
real(8), allocatable :: sigma1(:,:)
real(8), allocatable :: sigma2(:,:)
real(8), allocatable :: delta(:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: kerr(:)
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
allocate(a(3,3,nsymcrys))
allocate(s(nsymcrys))
d=1
if (spinorb) d=2
allocate(eps1(nwdos,d),eps2(nwdos,d))
allocate(sigma1(nwdos,d),sigma2(nwdos,d))
allocate(kerr(nwdos))
! allocate first-variational eigenvector array
allocate(evecfv(nmatmax,nstfv))
! allocate second-variational eigenvector array
allocate(evecsv(nstsv,nstsv))
! allocate the momentum matrix elements array
allocate(pmat(3,nstsv,nstsv))
! set up for generalised DFT correction if required
if (usegdft) then
  allocate(delta(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  call readstate
  call poteff
  call linengy
  call genapwfr
  call genlofr
end if
! precalculate for symmetrisation
do isym=1,nsymcrys
  a(:,:,isym)=dble(symcrys(:,:,isym))
  call r3mtm(a(1,1,isym),binv,a(1,1,isym))
  call r3mm(bvec,a(1,1,isym),a(1,1,isym))
  s(isym)=a(1,1,isym)*a(2,2,isym)-a(1,2,isym)*a(2,1,isym)
end do
! dielectric tensor components
i1=optcomp(1)
i2=optcomp(2)
if (spinorb) then
  i1=1
  i2=1
end if
d=1
10 continue
! find the record length
inquire(iolength=recl) pmat
open(50,file='PMAT.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
e(:,:)=0.d0
f(:,:)=0.d0
do ik=1,nkpt
! get the eigenvalues/vectors and occupancies from file
  call getevalsv(vkl(1,ik),evalsv(1,ik))
  call getoccsv(vkl(1,ik),occsv(1,ik))
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! compute generalised DFT correction if required
  if (usegdft) then
    call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
    call gdft(ik,apwalm,evecfv,evecsv,delta)
  end if
! read matrix elements from direct-access file
  read(50,rec=ik) pmat
  m=0
  do ist1=1,nstsv
    do ist2=1,nstsv
      m=m+1
! symmetrise the matrix elements
      sum=0.d0
      do isym=1,nsymcrys
        zt1(:)=0.d0
        if (i1.eq.i2) then
          do i=1,3
            zt1(i1)=zt1(i1)+a(i,i1,isym)*pmat(i,ist1,ist2)
          end do
        else
          do i=1,3
            zt1(i1)=zt1(i1)+a(i,i1,isym)*pmat(i,ist1,ist2)
            zt1(i2)=zt1(i2)+a(i,i2,isym)*pmat(i,ist1,ist2)
          end do
        end if
! calculate the desired dielectric tensor components
        if (i1.eq.i2) then
          sum=sum+dble(zt1(i1)*conjg(zt1(i1)))
        else
          if (s(isym).lt.0.d0) then
            sum=sum+aimag(zt1(i2)*conjg(zt1(i1)))
          else
            sum=sum+aimag(zt1(i1)*conjg(zt1(i2)))
          endif
        end if
 ! end of symmetrisation
      end do
      sum=sum/dble(nsymcrys)
      e12=evalsv(ist1,ik)-evalsv(ist2,ik)
! scissors correction
      if (evalsv(ist1,ik).gt.efermi) e12=e12+scissor
      if (evalsv(ist2,ik).gt.efermi) e12=e12-scissor
! generalised DFT correction
      if (usegdft) e12=e12+delta(ist1,ist2)
      e(m,ik)=e12
      f(m,ik)=(occsv(ist1,ik)-occsv(ist2,ik))*sum
    end do
  end do
end do
close(50)
t1=-4.d0*(pi**2)/omega
f(:,:)=t1*f(:,:)
! energy interval
wint(1)=0.d0
wint(2)=wintdos(2)
! generate energy grid
t1=(wint(2)-wint(1))/dble(nwdos)
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wint(1)
end do
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do
! integrate over the Brillouin zone
call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wint,n,n,e,f,eps2(1,d))
do iw=1,nwdos
  t1=w(iw)-scissor
  if (abs(t1).gt.eps) then
    eps2(iw,d)=eps2(iw,d)/(t1**2)
  end if
end do
! calculate the real part of dielectric function
if (i1.eq.i2) then
  t1=1.d0
else
  t1=0.d0
end if
do iw=1,nwdos
  do jw=1,nwdos
    t2=w(jw)**2-w(iw)**2
    if (abs(t2).gt.eps) then
      fw(jw)=w(jw)*eps2(jw,d)/t2
    else
      fw(jw)=0.d0
    end if
  end do
  call fderiv(-1,nwdos,w,fw,g,cf)
  eps1(iw,d)=t1+(2.d0/pi)*g(nwdos)
end do
sigma1(:,d)=eps2(:,d)*w(:)/(4.d0*pi)
sigma2(:,d)=-(eps1(:,d)-t1)*w(:)/(4.d0*pi)
! calculate the second component for MOKE
if ((spinorb).and.(d.eq.1)) then
  d=2
  i1=1
  i2=2
  goto 10
end if
! calculate the MOKE
if (spinorb) then
  do iw=2,nwdos
    zt1(1)=sigma1(iw,2)+zi*sigma2(iw,2)
    zt1(2)=(sigma1(iw,1)+zi*sigma2(iw,1))*sqrt(eps1(iw,1)+zi*eps2(iw,1))
    if (abs(zt1(2)).gt.0.d0) then
      kerr(iw)=-zt1(1)/zt1(2)
    else
      kerr(iw)=0.d0
    end if
  end do
end if
! write the dielectric tensor to file
open(50,file='EPSILON.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),eps1(iw,1)
end do
write(50,'("     ")')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),eps2(iw,1)
end do
! write the second component for MOKE
if (spinorb) then
  write(50,'("     ")')
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw),eps1(iw,2)
  end do
  write(50,'("     ")')
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw),eps2(iw,2)
  end do
end if
! write the optical conductivity
open(51,file='SIGMA.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  write(51,'(2G18.10)') w(iw),sigma1(iw,1)
end do
write(51,'("     ")')
do iw=1,nwdos
  write(51,'(2G18.10)') w(iw),sigma2(iw,1)
end do
! write the second component for MOKE
if (spinorb) then
  write(51,'("     ")')
  do iw=1,nwdos
    write(51,'(2G18.10)') w(iw),sigma1(iw,2)
  end do
  write(51,'("     ")')
  do iw=1,nwdos
    write(51,'(2G18.10)') w(iw),sigma2(iw,2)
  end do
end if
close(50)
close(51)
if (spinorb) then
  open(50,file='KERR.OUT',action='WRITE',form='FORMATTED')
  do iw=2,nwdos
    write(50,'(2G18.10)') w(iw),dble(kerr(iw))
  end do
  write(50,'("     ")')
  do iw=2,nwdos
    write(50,'(2G18.10)') w(iw),aimag(kerr(iw))
  end do
  close(50)
end if
write(*,*)
write(*,'("Info(linopt):")')
write(*,'(" element (",I1,",",I1,") of the dielectric tensor")') i1,i2
write(*,'(" written to EPSILON.OUT")')
write(*,*)
write(*,'(" element (",I1,",",I1,") of the optical conductivity")') i1,i2
write(*,'(" written to SIGMA.OUT")')
if (spinorb) then
  write(*,*)
  write(*,'(" Kerr angle in radians written to KERR.OUT")')
end if
write(*,*)
deallocate(w,fw,g,cf,e,f,a,s)
deallocate(eps1,eps2,sigma1,sigma2)
deallocate(evecfv,evecsv,kerr,pmat)
if (usegdft) deallocate(delta,apwalm)
return
end subroutine

