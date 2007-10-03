
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnsv(ik,apwalm,evalfv,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer is,ia,ias,ist1,ist2,ir,irc,ispn
integer nsc,i,j,k,lm,igk,ifg,lwork,info
! fine structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
! electron g factor
real(8), parameter :: ge=2.0023193043718d0
real(8), parameter :: ga4=ge*alpha/4.d0
real(8), parameter :: a24=alpha**2/4.d0
real(8) rm,t1
real(8) cpu0,cpu1
! automatic arrays
complex(8) zftp1(lmmaxvr),zftp2(lmmaxvr)
complex(8) zlflm(lmmaxvr,3)
! allocatable arrays
real(8), allocatable :: bmt(:,:,:)
real(8), allocatable :: bir(:,:)
real(8), allocatable :: vr(:)
real(8), allocatable :: drv(:)
real(8), allocatable :: cf(:,:)
real(8), allocatable :: sor(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: wfmt1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: zfft1(:)
complex(8), allocatable :: zfft2(:)
complex(8), allocatable :: zv(:,:)
complex(8), allocatable :: work(:)
! external functions
complex(8) zdotc,zfmtinp
external zdotc,zfmtinp
! spin-unpolarised case
if (.not.spinpol) then
  do i=1,nstsv
    evalsv(i,ik)=evalfv(i)
    spnchr(1,i,ik)=1.d0
  end do
  evecsv(:,:)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
! number of spin combinations after application of Hamiltonian
if ((ndmag.eq.3).or.(spinorb)) then
  nsc=4
else
  nsc=2
end if
call cpu_time(cpu0)
allocate(bmt(lmmaxvr,nrcmtmax,3))
allocate(bir(ngrtot,3))
allocate(vr(nrmtmax))
allocate(drv(nrmtmax))
allocate(cf(3,nrmtmax))
allocate(sor(nrcmtmax))
allocate(rwork(3*nstsv))
allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv))
allocate(wfmt2(lmmaxvr,nrcmtmax,nsc))
allocate(zfft1(ngrtot))
allocate(zfft2(ngrtot))
allocate(zv(ngkmax,nsc))
lwork=2*nstsv
allocate(work(lwork))
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! exchange-correlation magnetic field in spherical coordinates
    if (ndmag.eq.3) then
! non-collinear
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        do i=1,3
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw, &
           bxcmt(1,ir,ias,i),1,0.d0,bmt(1,irc,i),1)
        end do
      end do
    else
! collinear
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        bmt(:,irc,1:2)=0.d0
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw, &
         bxcmt(1,ir,ias,1),1,0.d0,bmt(1,irc,3),1)
      end do
    end if
! external muffin-tin magnetic field
    do irc=1,nrcmt(is)
      do i=1,3
        bmt(:,irc,i)=bmt(:,irc,i)+ga4*(bfcmt(i,ia,is)+bfieldc(i))
      end do
    end do
! spin-orbit radial function
    if (spinorb) then
! radial derivative of the spherical part of the potential
      vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias)*y00
      call fderiv(1,nrmt(is),spr(1,is),vr,drv,cf)
! spin-orbit coupling prefactor
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        rm=1.d0-0.5d0*(alpha**2)*vr(ir)
        sor(irc)=a24*drv(ir)/(spr(ir,is)*rm**2)
      end do
    end if
! compute the first-variational wavefunctions
    do ist1=1,nstfv
      call wavefmt(lradstp,lmaxvr,is,ia,ngk(ik,1),apwalm,evecfv(1,ist1), &
       lmmaxvr,wfmt1(1,1,ist1))
    end do
    do ist2=1,nstfv
      do irc=1,nrcmt(is)
! convert wavefunction to spherical coordinates
        call zgemv('N',lmmaxvr,lmmaxvr,zone,zbshtapw,lmmaxapw, &
         wfmt1(1,irc,ist2),1,zzero,zftp1,1)
! apply effective magnetic field and convert to spherical harmonics
        zftp2(:)=zftp1(:)*bmt(:,irc,3)
        call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
         wfmt2(1,irc,1),1)
        wfmt2(:,irc,2)=-wfmt2(:,irc,1)
        if (nsc.eq.4) then
          zftp2(:)=zftp1(:)*cmplx(bmt(:,irc,1),bmt(:,irc,2),8)
          call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
           wfmt2(1,irc,3),1)
          zftp2(:)=zftp1(:)*cmplx(bmt(:,irc,1),-bmt(:,irc,2),8)
          call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
           wfmt2(1,irc,4),1)
        end if
! apply spin-orbit coupling if required
        if (spinorb) then
          call lopzflm(lmaxvr,wfmt1(1,irc,ist2),lmmaxvr,zlflm)
          t1=sor(irc)
          do lm=1,lmmaxvr
            wfmt2(lm,irc,1)=wfmt2(lm,irc,1)+t1*zlflm(lm,3)
            wfmt2(lm,irc,2)=wfmt2(lm,irc,2)-t1*zlflm(lm,3)
            wfmt2(lm,irc,3)=wfmt2(lm,irc,3)+t1*(zlflm(lm,1)+zi*zlflm(lm,2))
            wfmt2(lm,irc,4)=wfmt2(lm,irc,4)+t1*(zlflm(lm,1)-zi*zlflm(lm,2))
          end do
        end if
      end do
! second-variational Hamiltonian matrix
      do ist1=1,nstfv
        do k=1,nsc
          if (k.eq.1) then
            i=ist1
            j=ist2
          else if (k.eq.2) then
            i=ist1+nstfv
            j=ist2+nstfv
          else if (k.eq.3) then
            i=ist1+nstfv
            j=ist2
          else
            i=ist1
            j=ist2+nstfv
          end if
          if (i.le.j) then
            evecsv(i,j)=evecsv(i,j)+zfmtinp(lmaxmat,nrcmt(is),rcmt(1,is), &
             lmmaxvr,wfmt1(1,1,ist1),wfmt2(1,1,k))
          end if
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
if (ndmag.eq.3) then
! non-collinear
  do ir=1,ngrtot
    bir(ir,:)=(bxcir(ir,:)+ga4*bfieldc(:))*cfunir(ir)
  end do
else
! collinear
  do ir=1,ngrtot
    bir(ir,1:2)=0.d0
    bir(ir,3)=(bxcir(ir,1)+ga4*bfieldc(3))*cfunir(ir)
  end do
end if
do ist2=1,nstfv
  zfft1(:)=0.d0
  do igk=1,ngk(ik,1)
    ifg=igfft(igkig(igk,ik,1))
    zfft1(ifg)=evecfv(igk,ist2)
  end do
! Fourier transform wavefunction to real-space
  call zfftifc(3,ngrid,1,zfft1)
  do k=1,nsc
    if (k.eq.1) then
      zfft2(:)=zfft1(:)*bir(:,3)
    else if (k.eq.2) then
      zfft2(:)=-zfft1(:)*bir(:,3)
    else if (k.eq.3) then
      zfft2(:)=zfft1(:)*(bir(:,1)+zi*bir(:,2))
    else
      zfft2(:)=zfft1(:)*(bir(:,1)-zi*bir(:,2))
    end if
    call zfftifc(3,ngrid,-1,zfft2)
    do igk=1,ngk(ik,1)
      ifg=igfft(igkig(igk,ik,1))
      zv(igk,k)=zfft2(ifg)
    end do
  end do
  do ist1=1,nstfv
    do k=1,nsc
      if (k.eq.1) then
        i=ist1
        j=ist2
      else if (k.eq.2) then
        i=ist1+nstfv
        j=ist2+nstfv
      else if (k.eq.3) then
        i=ist1+nstfv
        j=ist2
      else
        i=ist1
        j=ist2+nstfv
      end if
      if (i.le.j) then
        evecsv(i,j)=evecsv(i,j)+zdotc(ngk(ik,1),evecfv(1,ist1),1,zv(1,k),1)
      end if
    end do
  end do
end do
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist1=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist1)
  end do
end do
! diagonalise second-variational Hamiltonian
if ((ndmag.eq.3).or.(spinorb)) then
! non-collinear: full diagonalisation
  call zheev('V','U',nstsv,evecsv,nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 10
else
! collinear: block diagonalise H
  call zheev('V','U',nstfv,evecsv,nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 10
  i=nstfv+1
  call zheev('V','U',nstfv,evecsv(i,i),nstsv,evalsv(i,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 10
  do i=1,nstfv
    do j=1,nstfv
      evecsv(i,j+nstfv)=0.d0
      evecsv(i+nstfv,j)=0.d0
    end do
  end do
end if
deallocate(bmt,bir,vr,drv,cf,sor,rwork)
deallocate(wfmt1,wfmt2,zfft1,zfft2,zv,work)
call cpu_time(cpu1)
!$OMP CRITICAL
timesv=timesv+cpu1-cpu0
!$OMP END CRITICAL
return
10 continue
write(*,*)
write(*,'("Error(seceqnsv): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ik
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
stop
end subroutine
