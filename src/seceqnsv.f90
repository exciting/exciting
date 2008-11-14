
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl
! F. Bultmark, F. Cricchio and L. Nordstrom.
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
integer ispn,jspn,ia,is,ias
integer ist,jst,i,j,k,l,lm,nm
integer ir,irc,igk,ifg
integer nsc,lwork,info
! fine structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
! electron g factor
real(8), parameter :: ge=2.0023193043718d0
real(8), parameter :: ga4=ge*alpha/4.d0
real(8), parameter :: a24=alpha**2/4.d0
real(8) rm,t1
real(8) ts0,ts1
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
if ((.not.spinpol).and.(ldapu.eq.0)) then
  do i=1,nstsv
    evalsv(i,ik)=evalfv(i)
  end do
  evecsv(:,:)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
! number of spin combinations after application of Hamiltonian
if (spinpol) then
  if ((ncmag).or.(spinorb)) then
    nsc=3
  else
    nsc=2
  end if
else
  nsc=1
end if
call timesec(ts0)
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
    if (spinpol) then
! exchange-correlation magnetic field in spherical coordinates
      if (ncmag) then
! non-collinear
        irc=0
        do ir=1,nrmt(is),lradstp
          irc=irc+1
          do i=1,3
            call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
             bxcmt(:,ir,ias,i),1,0.d0,bmt(:,irc,i),1)
          end do
        end do
      else
! collinear
        irc=0
        do ir=1,nrmt(is),lradstp
          irc=irc+1
          bmt(:,irc,1:2)=0.d0
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
           bxcmt(:,ir,ias,1),1,0.d0,bmt(:,irc,3),1)
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
        call fderiv(1,nrmt(is),spr(:,is),vr,drv,cf)
! spin-orbit coupling prefactor
        irc=0
        do ir=1,nrmt(is),lradstp
          irc=irc+1
          rm=1.d0-0.5d0*(alpha**2)*vr(ir)
          sor(irc)=a24*drv(ir)/(spr(ir,is)*rm**2)
        end do
      end if
    end if
! compute the first-variational wavefunctions
    do ist=1,nstfv
      call wavefmt(lradstp,lmaxvr,is,ia,ngk(1,ik),apwalm,evecfv(:,ist), &
       lmmaxvr,wfmt1(:,:,ist))
    end do
! begin loop over states
    do jst=1,nstfv
      if (spinpol) then
        do irc=1,nrcmt(is)
! convert wavefunction to spherical coordinates
          call zgemv('N',lmmaxvr,lmmaxvr,zone,zbshtvr,lmmaxvr, &
           wfmt1(:,irc,jst),1,zzero,zftp1,1)
! apply effective magnetic field and convert to spherical harmonics
          zftp2(:)=zftp1(:)*bmt(:,irc,3)
          call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
           wfmt2(:,irc,1),1)
          wfmt2(:,irc,2)=-wfmt2(:,irc,1)
          if (nsc.eq.3) then
            zftp2(:)=zftp1(:)*cmplx(bmt(:,irc,1),-bmt(:,irc,2),8)
            call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
             wfmt2(:,irc,3),1)
          end if
! apply spin-orbit coupling if required
          if (spinorb) then
            call lopzflm(lmaxvr,wfmt1(:,irc,jst),lmmaxvr,zlflm)
            t1=sor(irc)
            do lm=1,lmmaxvr
              wfmt2(lm,irc,1)=wfmt2(lm,irc,1)+t1*zlflm(lm,3)
              wfmt2(lm,irc,2)=wfmt2(lm,irc,2)-t1*zlflm(lm,3)
              wfmt2(lm,irc,3)=wfmt2(lm,irc,3)+t1*(zlflm(lm,1)-zi*zlflm(lm,2))
            end do
          end if
        end do
      else
        wfmt2(:,:,:)=0.d0
      end if
! apply LDA+U potential if required
      if ((ldapu.ne.0).and.(llu(is).ge.0)) then
        l=llu(is)
        nm=2*l+1
        lm=idxlm(l,-l)
        do k=1,nsc
          if (k.eq.1) then
            ispn=1
            jspn=1
          else if (k.eq.2) then
            ispn=2
            jspn=2
          else 
            ispn=1
            jspn=2
          end if
          call zgemm('N','N',nm,nrcmt(is),nm,zone,vmatlu(lm,lm,ispn,jspn,ias), &
           lmmaxlu,wfmt1(lm,1,jst),lmmaxvr,zone,wfmt2(lm,1,k),lmmaxvr)
        end do 
      end if
! second-variational Hamiltonian matrix
      do ist=1,nstfv
        do k=1,nsc
          if (k.eq.1) then
            i=ist
            j=jst
          else if (k.eq.2) then
            i=ist+nstfv
            j=jst+nstfv
          else
            i=ist
            j=jst+nstfv
          end if
          if (i.le.j) then
            evecsv(i,j)=evecsv(i,j)+zfmtinp(.true.,lmaxmat,nrcmt(is), &
             rcmt(:,is),lmmaxvr,wfmt1(:,:,ist),wfmt2(:,:,k))
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
if (spinpol) then
  if (ncmag) then
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
  do jst=1,nstfv
    zfft1(:)=0.d0
    do igk=1,ngk(1,ik)
      ifg=igfft(igkig(igk,1,ik))
      zfft1(ifg)=evecfv(igk,jst)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngrid,1,zfft1)
! multiply with magnetic field and transform to G-space
    zfft2(:)=zfft1(:)*bir(:,3)
    call zfftifc(3,ngrid,-1,zfft2)
    do igk=1,ngk(1,ik)
      ifg=igfft(igkig(igk,1,ik))
      zv(igk,1)=zfft2(ifg)
      zv(igk,2)=-zfft2(ifg)
    end do
    if (nsc.eq.3) then
      zfft2(:)=zfft1(:)*cmplx(bir(:,1),-bir(:,2),8)
      call zfftifc(3,ngrid,-1,zfft2)
      do igk=1,ngk(1,ik)
        ifg=igfft(igkig(igk,1,ik))
        zv(igk,3)=zfft2(ifg)
      end do
    end if
! add to Hamiltonian matrix
    do ist=1,nstfv
      do k=1,nsc
        if (k.eq.1) then
          i=ist
          j=jst
        else if (k.eq.2) then
          i=ist+nstfv
          j=jst+nstfv
        else
          i=ist
          j=jst+nstfv
        end if
        if (i.le.j) then
          evecsv(i,j)=evecsv(i,j)+zdotc(ngk(1,ik),evecfv(:,ist),1,zv(:,k),1)
        end if
      end do
    end do
  end do
end if
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist)
  end do
end do
! diagonalise second-variational Hamiltonian
if (ndmag.eq.1) then
! collinear: block diagonalise H
  call zheev('V','U',nstfv,evecsv,nstsv,evalsv(:,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  i=nstfv+1
  call zheev('V','U',nstfv,evecsv(i,i),nstsv,evalsv(i,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  do i=1,nstfv
    do j=1,nstfv
      evecsv(i,j+nstfv)=0.d0
      evecsv(i+nstfv,j)=0.d0
    end do
  end do
else
! non-collinear or spin-unpolarised: full diagonalisation
  call zheev('V','U',nstsv,evecsv,nstsv,evalsv(:,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
end if
deallocate(bmt,bir,vr,drv,cf,sor,rwork)
deallocate(wfmt1,wfmt2,zfft1,zfft2,zv,work)
call timesec(ts1)
!$OMP CRITICAL
timesv=timesv+ts1-ts0
!$OMP END CRITICAL
return
20 continue
write(*,*)
write(*,'("Error(seceqnsv): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ik
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
stop
end subroutine
