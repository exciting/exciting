
! Copyright (C) 2006 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnss(ik,apwalm,evalfv,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
real(8), intent(in) :: evalfv(nstfv,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn,jspn,is,ia,ias
integer ist,jst,i,j,k,l,lm,nm
integer ir,irc,igk,ifg
integer lwork,info
! fine structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
! electron g factor
real(8), parameter :: ge=2.0023193043718d0
real(8), parameter :: ga4=ge*alpha/4.d0
real(8) ts0,ts1
! automatic arrays
complex(8) zftp1(lmmaxvr,nspnfv),zftp2(lmmaxvr)
! allocatable arrays
real(8), allocatable :: bmt(:,:,:)
real(8), allocatable :: bir(:,:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: zfft1(:,:)
complex(8), allocatable :: zfft2(:)
complex(8), allocatable :: zv(:,:)
complex(8), allocatable :: work(:)
! external functions
complex(8) zdotc
complex(8) zfmtinp
external zdotc,zfmtinp
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(seceqnss): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
call timesec(ts0)
allocate(bmt(lmmaxvr,nrcmtmax,3))
allocate(bir(ngrtot,3))
allocate(rwork(3*nstsv))
allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax,3))
allocate(zfft1(ngrtot,nspnfv))
allocate(zfft2(ngrtot))
allocate(zv(ngkmax,3))
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
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      do i=1,3
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,bxcmt(:,ir,ias,i), &
         1,0.d0,bmt(:,irc,i),1)
      end do
    end do
! external muffin-tin magnetic field
    do irc=1,nrcmt(is)
      do i=1,3
        bmt(:,irc,i)=bmt(:,irc,i)+ga4*(bfcmt(i,ia,is)+bfieldc(i))
      end do
    end do
! compute the first-variational wavefunctions
    do ispn=1,nspnfv
      do ist=1,nstfv
        call wavefmt(lradstp,lmaxvr,is,ia,ngk(ispn,ik),apwalm(:,:,:,:,ispn), &
         evecfv(:,ist,ispn),lmmaxvr,wfmt1(:,:,ist,ispn))
      end do
    end do
    do jst=1,nstfv
      do irc=1,nrcmt(is)
! convert wavefunctions to spherical coordinates
        do ispn=1,nspnfv
          call zgemv('N',lmmaxvr,lmmaxvr,zone,zbshtvr,lmmaxvr, &
           wfmt1(:,irc,jst,ispn),1,zzero,zftp1(:,ispn),1)
        end do
! apply effective magnetic field and convert to spherical harmonics
        zftp2(:)=zftp1(:,1)*bmt(:,irc,3)
        call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
         wfmt2(:,irc,1),1)
        zftp2(:)=-zftp1(:,2)*bmt(:,irc,3)
        call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
         wfmt2(:,irc,2),1)
        zftp2(:)=zftp1(:,2)*cmplx(bmt(:,irc,1),-bmt(:,irc,2),8)
        call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zftp2,1,zzero, &
         wfmt2(:,irc,3),1)
      end do
! apply LDA+U potential if required
      if ((ldapu.ne.0).and.(llu(is).ge.0)) then
        l=llu(is)
        nm=2*l+1
        lm=idxlm(l,-l)
        do k=1,3
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
           lmmaxlu,wfmt1(lm,1,jst,jspn),lmmaxvr,zone,wfmt2(lm,1,k),lmmaxvr)
        end do
      end if
! second-variational Hamiltonian matrix
      do ist=1,nstfv
        do k=1,3
          if (k.eq.1) then
            ispn=1
            i=ist
            j=jst
          else if (k.eq.2) then
            ispn=2
            i=ist+nstfv
            j=jst+nstfv
          else
            ispn=1
            i=ist
            j=jst+nstfv
          end if
          if (i.le.j) then
            evecsv(i,j)=evecsv(i,j)+zfmtinp(.true.,lmaxmat,nrcmt(is), &
             rcmt(:,is),lmmaxvr,wfmt1(:,:,ist,ispn),wfmt2(:,:,k))
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
do ir=1,ngrtot
  bir(ir,:)=(bxcir(ir,:)+ga4*bfieldc(:))*cfunir(ir)
end do
do jst=1,nstfv
  do ispn=1,nspnfv
    zfft1(:,ispn)=0.d0
    do igk=1,ngk(ispn,ik)
      ifg=igfft(igkig(igk,ispn,ik))
      zfft1(ifg,ispn)=evecfv(igk,jst,ispn)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngrid,1,zfft1(:,ispn))
  end do
! multiply with magnetic field and transform to G-space
  do k=1,3
    if (k.eq.1) then
      ispn=1
      zfft2(:)=zfft1(:,1)*bir(:,3)
    else if (k.eq.2) then
      ispn=2
      zfft2(:)=-zfft1(:,2)*bir(:,3)
    else
      ispn=1
      zfft2(:)=zfft1(:,2)*cmplx(bir(:,1),-bir(:,2),8)
    end if
    call zfftifc(3,ngrid,-1,zfft2)
    do igk=1,ngk(ispn,ik)
      ifg=igfft(igkig(igk,ispn,ik))
      zv(igk,k)=zfft2(ifg)
    end do
  end do
  do ist=1,nstfv
    do k=1,3
      if (k.eq.1) then
        ispn=1
        i=ist
        j=jst
      else if (k.eq.2) then
        ispn=2
        i=ist+nstfv
        j=jst+nstfv
      else
        ispn=1
        i=ist
        j=jst+nstfv
      end if
      if (i.le.j) then
        evecsv(i,j)=evecsv(i,j) &
         +zdotc(ngk(ispn,ik),evecfv(:,ist,ispn),1,zv(:,k),1)
      end if
    end do
  end do
end do
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist,ispn)
  end do
end do
! diagonalise the Hamiltonian
call zheev('V','U',nstsv,evecsv,nstsv,evalsv(:,ik),work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnss): diagonalisation of the second-variational &
   &Hamiltonian failed")')
  write(*,'(" for k-point ",I8)') ik
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(bmt,bir,rwork)
deallocate(wfmt1,wfmt2,zfft1,zfft2,work)
call timesec(ts1)
!$OMP CRITICAL
timesv=timesv+ts1-ts0
!$OMP END CRITICAL
return
end subroutine

