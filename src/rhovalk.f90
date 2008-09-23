
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhovalk
! !INTERFACE:
subroutine rhovalk(ik,evecfv,evecsv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the partial valence charge density from the eigenvectors at
!   $k$-point {\tt ik}. In the muffin-tin region, the wavefunction is obtained
!   in terms of its $(l,m)$-components from both the APW and local-orbital
!   functions. Using a backward spherical harmonic transform (SHT), the
!   wavefunction is converted to real-space and the density obtained from its
!   modulus squared. This density is then transformed with a forward SHT and
!   accumulated in the global variable {\tt rhomt}. A similar proccess is used
!   for the intersitial density in which the wavefunction in real-space is
!   obtained from a Fourier transform of the sum of APW functions. The
!   interstitial density is added to the global array {\tt rhoir}. See routines
!   {\tt wavefmt}, {\tt genshtmat} and {\tt seceqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer nsd,ispn,jspn,is,ia,ias,ist
integer ir,irc,itp,igk,ifg,i,j,n
real(8) wo,t1,t2,t3
real(8) ts0,ts1
complex(8) zt1,zt2,zt3
! allocatable arrays
logical, allocatable :: done(:,:)
real(8), allocatable :: rflm(:,:)
real(8), allocatable :: rfmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfmt2(:,:,:,:)
complex(8), allocatable :: wfmt3(:,:,:)
complex(8), allocatable :: zfft(:,:)
call timesec(ts0)
if (spinpol) then
  if (ncmag) then
    nsd=4
  else
    nsd=2
  end if
else
  nsd=1
end if
allocate(done(nstfv,nspnfv))
allocate(rflm(lmmaxvr,nsd))
allocate(rfmt(lmmaxvr,nrcmtmax,nsd))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(wfmt1(lmmaxvr,nrcmtmax))
if (tevecsv) allocate(wfmt2(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt3(lmmaxvr,nrcmtmax,nspinor))
allocate(zfft(ngrtot,nspinor))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
!----------------------------!
!     muffin-tin density     !
!----------------------------!
do is=1,nspecies
  n=lmmaxvr*nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    done(:,:)=.false.
    rfmt(:,:,:)=0.d0
    do j=1,nstsv
      wo=wkpt(ik)*occsv(j,ik)
      if (abs(wo).gt.epsocc) then
        if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
          wfmt3(:,:,:)=0.d0
          i=0
          do ispn=1,nspinor
            if (spinsprl) then
              jspn=ispn
            else
              jspn=1
            end if
            do ist=1,nstfv
              i=i+1
              zt1=evecsv(i,j)
              if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                if (.not.done(ist,jspn)) then
                  call wavefmt(lradstp,lmaxvr,is,ia,ngk(jspn,ik), &
                   apwalm(:,:,:,:,jspn),evecfv(:,ist,jspn),lmmaxvr,wfmt1)
! convert from spherical harmonics to spherical coordinates
                  call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr, &
                   lmmaxvr,wfmt1,lmmaxvr,zzero,wfmt2(:,:,ist,jspn),lmmaxvr)
                  done(ist,jspn)=.true.
                end if
! add to spinor wavefunction
                call zaxpy(n,zt1,wfmt2(:,:,ist,jspn),1,wfmt3(:,:,ispn),1)
              end if
            end do
          end do
        else
! spin-unpolarised wavefunction
          call wavefmt(lradstp,lmaxvr,is,ia,ngk(1,ik),apwalm,evecfv(:,j,1), &
           lmmaxvr,wfmt1)
! convert from spherical harmonics to spherical coordinates
          call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr,lmmaxvr, &
           wfmt1,lmmaxvr,zzero,wfmt3,lmmaxvr)
        end if
! add to the spin density matrix
        if (spinpol) then
! spin-polarised
          do irc=1,nrcmt(is)
            do itp=1,lmmaxvr
              zt1=wfmt3(itp,irc,1)
              zt2=wfmt3(itp,irc,2)
              zt3=zt1*conjg(zt2)
              rfmt(itp,irc,1)=rfmt(itp,irc,1)+wo*(dble(zt1)**2+aimag(zt1)**2)
              rfmt(itp,irc,2)=rfmt(itp,irc,2)+wo*(dble(zt2)**2+aimag(zt2)**2)
              if (ncmag) then
                rfmt(itp,irc,3)=rfmt(itp,irc,3)+wo*dble(zt3)
                rfmt(itp,irc,4)=rfmt(itp,irc,4)+wo*aimag(zt3)
              end if
            end do
          end do
        else
! spin-unpolarised
          do irc=1,nrcmt(is)
            do itp=1,lmmaxvr
              zt1=wfmt3(itp,irc,1)
              rfmt(itp,irc,1)=rfmt(itp,irc,1)+wo*(dble(zt1)**2+aimag(zt1)**2)
            end do
          end do
        end if
      end if
    end do
! convert to spherical harmonics and add to rhomt and magmt
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      do i=1,nsd
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt(:,irc,i),1, &
         0.d0,rflm(:,i),1)
      end do
!$OMP CRITICAL
      if (spinpol) then
! spin-polarised
        if (ncmag) then
          magmt(:,ir,ias,1)=magmt(:,ir,ias,1)+2.d0*rflm(:,3)
          magmt(:,ir,ias,2)=magmt(:,ir,ias,2)-2.d0*rflm(:,4)
          magmt(:,ir,ias,3)=magmt(:,ir,ias,3)+rflm(:,1)-rflm(:,2)
        else
          magmt(:,ir,ias,1)=magmt(:,ir,ias,1)+rflm(:,1)-rflm(:,2)
        end if
        rhomt(:,ir,ias)=rhomt(:,ir,ias)+rflm(:,1)+rflm(:,2)
      else
! spin-unpolarised
        rhomt(:,ir,ias)=rhomt(:,ir,ias)+rflm(:,1)
      end if
!$OMP END CRITICAL
    end do
  end do
end do
!------------------------------!
!     interstitial density     !
!------------------------------!
do j=1,nstsv
  wo=wkpt(ik)*occsv(j,ik)
  if (abs(wo).gt.epsocc) then
    t1=wo/omega
    zfft(:,:)=0.d0
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        if (spinsprl) then
          jspn=ispn
        else
          jspn=1
        end if
        do ist=1,nstfv
          i=i+1
          zt1=evecsv(i,j)
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
            do igk=1,ngk(jspn,ik)
              ifg=igfft(igkig(igk,jspn,ik))
              zfft(ifg,ispn)=zfft(ifg,ispn)+zt1*evecfv(igk,ist,jspn)
            end do
          end if
        end do
      end do
    else
! spin-unpolarised wavefunction
      do igk=1,ngk(1,ik)
        ifg=igfft(igkig(igk,1,ik))
        zfft(ifg,1)=evecfv(igk,j,1)
      end do
    end if
! Fourier transform wavefunction to real-space
    do ispn=1,nspinor
      call zfftifc(3,ngrid,1,zfft(:,ispn))
    end do
!$OMP CRITICAL
    if (spinpol) then
! spin-polarised
      do ir=1,ngrtot
        zt1=zfft(ir,1)
        zt2=zfft(ir,2)
        zt3=zt1*conjg(zt2)
        t2=dble(zt1)**2+aimag(zt1)**2
        t3=dble(zt2)**2+aimag(zt2)**2
        rhoir(ir)=rhoir(ir)+t1*(t2+t3)
        if (ncmag) then
          magir(ir,1)=magir(ir,1)+2.d0*t1*dble(zt3)
          magir(ir,2)=magir(ir,2)-2.d0*t1*aimag(zt3)
          magir(ir,3)=magir(ir,3)+t1*(t2-t3)
        else
          magir(ir,1)=magir(ir,1)+t1*(t2-t3)
        end if
      end do
    else
! spin-unpolarised
      do ir=1,ngrtot
        zt1=zfft(ir,1)
        rhoir(ir)=rhoir(ir)+t1*(dble(zt1)**2+aimag(zt1)**2)
      end do
    end if
!$OMP END CRITICAL
  end if
end do
deallocate(done,rflm,rfmt,apwalm,wfmt1,wfmt3,zfft)
if (tevecsv) deallocate(wfmt2)
call timesec(ts1)
!$OMP CRITICAL
timerho=timerho+ts1-ts0
!$OMP END CRITICAL
return
end subroutine
!EOC
