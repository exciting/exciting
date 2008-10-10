
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwf2val(ik,evecfv,evecsv,gw2fmt,gw2fir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
real(8), intent(inout) :: gw2fmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(inout) :: gw2fir(ngrtot)
! local variables
integer ispn,ist,is,ia,ias
integer ir,itp,igk,ifg,i,j,n
real(8) wo,t1
complex(8) zt1
! automatic arrays
logical done(nstfv)
complex(8) zftp(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: gzfmt(:,:,:)
complex(8), allocatable :: zfft1(:,:)
complex(8), allocatable :: zfft2(:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt1(lmmaxvr,nrmtmax,nstfv))
allocate(wfmt2(lmmaxvr,nrmtmax,nspinor))
allocate(gzfmt(lmmaxvr,nrmtmax,3))
allocate(zfft1(ngrtot,nspinor))
allocate(zfft2(ngrtot))
! find the matching coefficients
call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
  n=lmmaxvr*nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    done(:)=.false.
    do j=1,nstsv
      wo=wkpt(ik)*occsv(j,ik)
      if (abs(wo).gt.epsocc) then
        if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
          wfmt2(:,:,:)=0.d0
          i=0
          do ispn=1,nspinor
            do ist=1,nstfv
              i=i+1
              zt1=evecsv(i,j)
              if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                if (.not.done(ist)) then
                  call wavefmt(1,lmaxvr,is,ia,ngk(1,ik),apwalm,evecfv(:,ist), &
                   lmmaxvr,wfmt1(:,:,ist))
                  done(ist)=.true.
                end if
! add to spinor wavefunction
                call zaxpy(n,zt1,wfmt1(:,:,ist),1,wfmt2(:,:,ispn),1)
              end if
            end do
          end do
        else
! spin-unpolarised wavefunction
          call wavefmt(1,lmaxvr,is,ia,ngk(1,ik),apwalm,evecfv(:,j),lmmaxvr, &
           wfmt2)
        end if
! compute the gradient of the wavefunction
        do ispn=1,nspinor
          call gradzfmt(lmaxvr,nrmt(is),spr(:,is),lmmaxvr,nrmtmax, &
           wfmt2(:,:,ispn),gzfmt)
! convert gradient from spherical harmonics to spherical coordinates
          do i=1,3
            do ir=1,nrmt(is)
              call zgemv('N',lmmaxvr,lmmaxvr,zone,zbshtvr,lmmaxvr, &
               gzfmt(:,ir,i),1,zzero,zftp,1)
              do itp=1,lmmaxvr
                t1=wo*(dble(zftp(itp))**2+aimag(zftp(itp))**2)
                gw2fmt(itp,ir,ias)=gw2fmt(itp,ir,ias)+t1
              end do
            end do
          end do
        end do
      end if
    end do
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
do j=1,nstsv
  wo=wkpt(ik)*occsv(j,ik)
  if (abs(wo).gt.epsocc) then
    t1=wo/omega
    zfft1(:,:)=0.d0
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        do ist=1,nstfv
          i=i+1
          zt1=evecsv(i,j)
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
            do igk=1,ngk(1,ik)
              ifg=igfft(igkig(igk,1,ik))
              zfft1(ifg,ispn)=zfft1(ifg,ispn)+zt1*evecfv(igk,ist)
            end do
          end if
        end do
      end do
    else
! spin-unpolarised wavefunction
      do igk=1,ngk(1,ik)
        ifg=igfft(igkig(igk,1,ik))
        zfft1(ifg,1)=evecfv(igk,j)
      end do
    end if
! compute gradient of wavefunction
    do ispn=1,nspinor
      do i=1,3
        zfft2(:)=0.d0
        do igk=1,ngk(1,ik)
          ifg=igfft(igkig(igk,1,ik))
          zfft2(ifg)=zi*vgkc(i,igk,1,ik)*zfft1(ifg,ispn)
        end do
! Fourier transform gradient to real-space
        call zfftifc(3,ngrid,1,zfft2)
        do ir=1,ngrtot
          gw2fir(ir)=gw2fir(ir)+t1*(dble(zfft2(ir))**2+aimag(zfft2(ir))**2)
        end do
      end do
    end do
  end if
end do
deallocate(apwalm,wfmt1,wfmt2,gzfmt,zfft1,zfft2)
return
end subroutine
