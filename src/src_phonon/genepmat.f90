
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genepmat(iq,vpl,dveffmt,dveffir,epmat)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: dveffmt(lmmaxapw,nrcmtmax,natmtot,3*natmtot)
complex(8), intent(in) :: dveffir(ngrtot,3*natmtot)
complex(8), intent(out) :: epmat(nstsv,nstsv,3*natmtot)
! local variables
integer is,ia,ias
integer ngp,ngpq,igp,ifg
integer nrc,irc,iv(3)
integer ist,jst,ispn
integer i,j,k,l,m,n
integer i1,i2,i3,ir
real(8) vpc(3),vpql(3),vpqc(3)
real(8) v1(3),v2(3),v3(3),t1
complex(8) zt1
! allocatable arrays
integer, allocatable :: igpig(:)
integer, allocatable :: igpqig(:)
real(8), allocatable :: vgpl(:,:)
real(8), allocatable :: vgpc(:,:)
real(8), allocatable :: gpc(:)
real(8), allocatable :: tpgpc(:,:)
real(8), allocatable :: vgpql(:,:)
real(8), allocatable :: vgpqc(:,:)
real(8), allocatable :: gpqc(:)
real(8), allocatable :: tpgpqc(:,:)
complex(8), allocatable :: sfacgp(:,:)
complex(8), allocatable :: sfacgpq(:,:)
complex(8), allocatable :: apwalm1(:,:,:,:)
complex(8), allocatable :: apwalm2(:,:,:,:)
complex(8), allocatable :: evecfv1(:,:)
complex(8), allocatable :: evecfv2(:,:)
complex(8), allocatable :: evecsv1(:,:)
complex(8), allocatable :: evecsv2(:,:)
complex(8), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wfmt3(:,:)
complex(8), allocatable :: zfir1(:)
complex(8), allocatable :: zfir2(:)
complex(8), allocatable :: zfir3(:)
complex(8), allocatable :: zv(:)
complex(8), allocatable :: epm(:,:,:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
n=3*natmtot
! allocate local arrays
allocate(igpig(ngkmax))
allocate(igpqig(ngkmax))
allocate(vgpl(3,ngkmax))
allocate(vgpc(3,ngkmax))
allocate(gpc(ngkmax))
allocate(tpgpc(2,ngkmax))
allocate(vgpql(3,ngkmax))
allocate(vgpqc(3,ngkmax))
allocate(gpqc(ngkmax))
allocate(tpgpqc(2,ngkmax))
allocate(sfacgp(ngkmax,natmtot))
allocate(sfacgpq(ngkmax,natmtot))
allocate(apwalm1(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(apwalm2(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv1(nmatmax,nstfv))
allocate(evecfv2(nmatmax,nstfv))
if (tevecsv) then
  allocate(evecsv1(nstsv,nstsv))
  allocate(evecsv2(nstsv,nstsv))
end if
allocate(wfmt1(lmmaxapw,nrcmtmax))
allocate(wfmt2(lmmaxapw,nrcmtmax,nstfv))
allocate(wfmt3(lmmaxapw,nrcmtmax))
allocate(zfir1(ngrtot))
allocate(zfir2(ngrtot))
allocate(zfir3(ngrtot))
allocate(zv(ngkmax))
allocate(epm(nstfv,nstfv,n))
! p-vector in Cartesian coordinates
call r3mv(bvec,vpl,vpc)
! generate the G+p vectors
call gengpvec(vpl,vpc,ngp,igpig,vgpl,vgpc,gpc,tpgpc)
! generate the structure factors
call gensfacgp(ngp,vgpc,ngkmax,sfacgp)
! find the matching coefficients for k-point p
call match(ngp,gpc,tpgpc,sfacgp,apwalm1)
! get the eigenvectors for k-point p
call getevecfv(vpl,vgpl,evecfv1)
! p+q-vector in lattice coordinates
vpql(:)=vpl(:)+vql(:,iq)
! map vector components to [0,1) interval
call r3frac(epslat,vpql,iv)
! p+q-vector in Cartesian coordinates
call r3mv(bvec,vpql,vpqc)
! generate the G+p+q-vectors
call gengpvec(vpql,vpqc,ngpq,igpqig,vgpql,vgpqc,gpqc,tpgpqc)
! generate the structure factors
call gensfacgp(ngpq,vgpqc,ngkmax,sfacgpq)
! find the matching coefficients for k-point p+q
call match(ngpq,gpqc,tpgpqc,sfacgpq,apwalm2)
! get the eigenvectors for k-point p+q
call getevecfv(vpql,vgpql,evecfv2)
! set the first-variational matrix element array to zero
epm(:,:,:)=0.d0
!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,nstfv
! calculate the wavefunction for k-point p+q
      call wavefmt(lradstp,lmaxapw,is,ia,ngpq,apwalm2,evecfv2(:,ist),lmmaxapw, &
       wfmt1)
! convert from spherical harmonics to spherical coordinates
      call zgemm('N','N',lmmaxapw,nrc,lmmaxapw,zone,zbshtapw,lmmaxapw,wfmt1, &
       lmmaxapw,zzero,wfmt2(:,:,ist),lmmaxapw)
    end do
    do jst=1,nstfv
! calculate the wavefunction for k-point p
      call wavefmt(lradstp,lmaxapw,is,ia,ngp,apwalm1,evecfv1(:,jst),lmmaxapw, &
       wfmt1)
! convert from spherical harmonics to spherical coordinates
      call zgemm('N','N',lmmaxapw,nrc,lmmaxapw,zone,zbshtapw,lmmaxapw,wfmt1, &
       lmmaxapw,zzero,wfmt3,lmmaxapw)
! loop over phonon branches
      do i=1,n
! multiply the wavefunction by the change in effective potential
        do irc=1,nrc
          wfmt1(:,irc)=wfmt3(:,irc)*dveffmt(:,irc,ias,i)
        end do
! add to the first-variational matrix elements
        do ist=1,nstfv
          epm(ist,jst,i)=epm(ist,jst,i)+zfmtinp(.false.,lmaxapw,nrc, &
           rcmt(:,is),lmmaxapw,wfmt2(:,:,ist),wfmt1)
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
! store G=q+p-p', where p' is the p+q-vector mapped to [0,1)
v1(:)=vqc(:,iq)+vpc(:)-vpqc(:)
! compute exp(i(q+p-p').r) for each r-vector on the grid
ir=0
do i3=0,ngrid(3)-1
  v2(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    v2(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      v2(1)=dble(i1)/dble(ngrid(1))
      ir=ir+1
      call r3mv(avec,v2,v3)
      t1=v1(1)*v3(1)+v1(2)*v3(2)+v1(3)*v3(3)
      zfir1(ir)=cmplx(cos(t1),sin(t1),8)
    end do
  end do
end do
! compute interstitial wavefunctions for k-point p
do jst=1,nstfv
  zfir2(:)=0.d0
  do igp=1,ngp
    ifg=igfft(igpig(igp))
    zfir2(ifg)=evecfv1(igp,jst)
  end do
! Fourier transform wavefunction to real-space
  call zfftifc(3,ngrid,1,zfir2)
! multiply with the phase factor
  zfir2(:)=zfir2(:)*zfir1(:)
! loop over phonon branches
  do i=1,n
! multiply the wavefunction with the change in effective potential
    zfir3(:)=zfir2(:)*dveffir(:,i)
! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,zfir3)
! store as wavefunction with G+p+q index
    do igp=1,ngpq
      ifg=igfft(igpqig(igp))
      zv(igp)=zfir3(ifg)
    end do
! add to the first-variational matrix elements
    do ist=1,nstfv
      epm(ist,jst,i)=epm(ist,jst,i)+zdotc(ngpq,evecfv2(:,ist),1,zv,1)
    end do
  end do
end do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
if (tevecsv) then
! get the second-variational eigenvectors
  call getevecsv(vpl,evecsv1)
  call getevecsv(vpql,evecsv2)
  epmat(:,:,:)=0.d0
  do i=1,nstsv
    do j=1,nstsv
      k=0
      do ispn=1,nspinor
        do ist=1,nstfv
          k=k+1
          l=(ispn-1)*nstfv
          do jst=1,nstfv
            l=l+1
            zt1=conjg(evecsv2(k,i))*evecsv1(l,j)
            do m=1,n
              epmat(i,j,m)=epmat(i,j,m)+epm(ist,jst,m)*zt1
            end do
          end do
        end do
      end do
    end do
  end do
else
  epmat(:,:,:)=epm(:,:,:)
end if
deallocate(igpig,igpqig,vgpl,vgpc,gpc,tpgpc,vgpql,vgpqc,gpqc,tpgpqc)
deallocate(sfacgp,sfacgpq,apwalm1,apwalm2,evecfv1,evecfv2)
if (tevecsv) deallocate(evecsv1,evecsv2)
deallocate(wfmt1,wfmt2,wfmt3,zfir1,zfir2,zfir3,zv,epm)
return
end subroutine

