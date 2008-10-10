
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpiqr(ik,emat)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: emat(nstsv,nstsv)
! local variables
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3,ist,jst
integer is,ia,i1,i2,i3,iv(3)
integer ngkq,igk,ifg,ir,irc
integer i,j,k,l,ispn
real(8) vecqc(3),qc,tp(2)
real(8) vkql(3),vkqc(3),x,t1
real(8) v1(3),v2(3),v3(3)
complex(8) zsum,zt1,zt2,zt3
! automatic arrays
complex(8) ylm(lmmaxvr),zl(0:lmaxvr)
complex(8) zflm(lmmaxvr)
! allocatable arrays
integer, allocatable :: igkqig(:)
real(8), allocatable :: gnt(:,:,:)
real(8), allocatable :: jlqr(:,:)
real(8), allocatable :: vgkql(:,:)
real(8), allocatable :: vgkqc(:,:)
real(8), allocatable :: gkqc(:)
real(8), allocatable :: tpgkqc(:,:)
complex(8), allocatable :: sfacgkq(:,:)
complex(8), allocatable :: apwalm1(:,:,:,:)
complex(8), allocatable :: apwalm2(:,:,:,:)
complex(8), allocatable :: evecfv1(:,:)
complex(8), allocatable :: evecfv2(:,:)
complex(8), allocatable :: evecsv1(:,:)
complex(8), allocatable :: evecsv2(:,:)
complex(8), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wfmt3(:,:)
complex(8), allocatable :: wfir(:)
complex(8), allocatable :: zfir1(:)
complex(8), allocatable :: zfir2(:)
complex(8), allocatable :: em(:,:)
! external functions
real(8) gaunt
complex(8) zfmtinp,zdotc
external gaunt,zfmtinp,zdotc
! check if q-vector is zero
t1=vecql(1)**2+vecql(2)**2+vecql(3)**2
if (t1.lt.epslat) then
  emat(:,:)=0.d0
  do i=1,nstsv
    emat(i,i)=1.d0
  end do
  return
end if
! check q-vector is commensurate with k-point grid
v1(:)=dble(ngridk(:))*vecql(:)
v2(:)=abs(v1(:)-nint(v1(:)))
if ((v2(1).gt.epslat).or.(v2(2).gt.epslat).or.(v2(3).gt.epslat)) then
  write(*,*)
  write(*,'("Error(genexpiqr): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! allocate local arrays
allocate(igkqig(ngkmax))
allocate(gnt(lmmaxvr,lmmaxvr,lmmaxvr))
allocate(jlqr(0:lmaxvr,nrcmtmax))
allocate(vgkql(3,ngkmax))
allocate(vgkqc(3,ngkmax))
allocate(gkqc(ngkmax))
allocate(tpgkqc(2,ngkmax))
allocate(sfacgkq(ngkmax,natmtot))
allocate(apwalm1(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(apwalm2(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv1(nmatmax,nstfv))
allocate(evecfv2(nmatmax,nstfv))
if (tevecsv) then
  allocate(evecsv1(nstsv,nstsv))
  allocate(evecsv2(nstsv,nstsv))
end if
allocate(wfmt1(lmmaxvr,nrcmtmax))
allocate(wfmt2(lmmaxvr,nrcmtmax,nstfv))
allocate(wfmt3(lmmaxvr,nrcmtmax))
allocate(wfir(ngkmax))
allocate(zfir1(ngrtot),zfir2(ngrtot))
allocate(em(nstfv,nstfv))
! compute the Gaunt coefficients
do l1=0,lmaxvr
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    do l2=0,lmaxvr
      do m2=-l2,l2
        lm2=idxlm(l2,m2)
        do l3=0,lmaxvr
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
            gnt(lm1,lm2,lm3)=gaunt(l1,l2,l3,m1,m2,m3)
          end do
        end do
      end do
    end do
  end do
end do
! q-vector in Cartesian coordinates
call r3mv(bvec,vecql,vecqc)
! length and spherical coordinates of q-vector
call sphcrd(vecqc,qc,tp)
! generate the conjugate spherical harmonics of the q-vector
call genylm(lmaxvr,tp,ylm)
ylm(:)=conjg(ylm(:))
! get the eigenvector for k-point k
call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv1)
! find the matching coefficients for k-point k
call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm1)
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ik)+vecql(:)
! map vector components to [0,1) interval
call r3frac(epslat,vkql,iv)
! k+q-vector in Cartesian coordinates
call r3mv(bvec,vkql,vkqc)
! generate the G+k+q-vectors
call gengpvec(vkql,vkqc,ngkq,igkqig,vgkql,vgkqc,gkqc,tpgkqc)
! generate the structure factors
call gensfacgp(ngkq,vgkqc,ngkmax,sfacgkq)
! find the matching coefficients for k-point k+q
call match(ngkq,gkqc,tpgkqc,sfacgkq,apwalm2)
! get the eigenvector for k-point k+q
call getevecfv(vkql,vgkql,evecfv2)
! set the first-variational matrix element array to zero
em(:,:)=0.d0
!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!
do is=1,nspecies
! compute the spherical Bessel functions
  do irc=1,nrcmt(is)
    x=qc*rcmt(irc,is)
    call sbessel(lmaxvr,x,jlqr(:,irc))
  end do
  do ia=1,natoms(is)
    t1=dot_product(vecqc(:),atposc(:,ia,is))
    zt1=fourpi*cmplx(cos(t1),sin(t1),8)
    do l1=0,lmaxvr
      zl(l1)=zt1*zil(l1)
    end do
    do ist=1,nstfv
! calculate the wavefunction for k-point k+q
      call wavefmt(lradstp,lmaxvr,is,ia,ngkq,apwalm2,evecfv2(:,ist), &
       lmmaxvr,wfmt2(:,:,ist))
    end do
    do jst=1,nstfv
! calculate the wavefunction for k-point k
      call wavefmt(lradstp,lmaxvr,is,ia,ngk(1,ik),apwalm1,evecfv1(:,jst), &
       lmmaxvr,wfmt1)
! multiply wavefunction with exp(iq.r)
      do irc=1,nrcmt(is)
        zflm(:)=0.d0
        do l2=0,lmaxvr
          zt1=zl(l2)*jlqr(l2,irc)
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
            zt2=zt1*ylm(lm2)
            do lm3=1,lmmaxvr
              zt3=zt2*wfmt1(lm3,irc)
              do lm1=1,lmmaxvr
                zflm(lm1)=zflm(lm1)+gnt(lm1,lm2,lm3)*zt3
              end do
            end do
          end do
        end do
        wfmt3(:,irc)=zflm(:)
      end do
      do ist=1,nstfv
        em(ist,jst)=em(ist,jst)+zfmtinp(.true.,lmaxvr,nrcmt(is),rcmt(:,is), &
         lmmaxvr,wfmt2(:,:,ist),wfmt3)
      end do
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
! store q+k-k', where k' is the k+q-vector mapped to [0,1)
v1(:)=vecqc(:)+vkc(:,ik)-vkqc(:)
! compute exp(i(q+k-k').r) times by the characteristic function
ir=0
do i3=0,ngrid(3)-1
  v2(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    v2(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      v2(1)=dble(i1)/dble(ngrid(1))
      ir=ir+1
      call r3mv(avec,v2,v3)
      t1=dot_product(v1(:),v3(:))
      zfir1(ir)=cfunir(ir)*cmplx(cos(t1),sin(t1),8)
    end do
  end do
end do
! compute interstitial wavefunctions for k-point k
do jst=1,nstfv
  zfir2(:)=0.d0
  do igk=1,ngk(1,ik)
    ifg=igfft(igkig(igk,1,ik))
    zfir2(ifg)=evecfv1(igk,jst)
  end do
! Fourier transform wavefunction to real-space
  call zfftifc(3,ngrid,1,zfir2)
! multiply with the phase and characteristic function
  zfir2(:)=zfir2(:)*zfir1(:)
! Fourier transform back to G-space
  call zfftifc(3,ngrid,-1,zfir2)
! store in wfir
  do igk=1,ngkq
    ifg=igfft(igkqig(igk))
    wfir(igk)=zfir2(ifg)
  end do
! add to the first-variational matrix elements
  do ist=1,nstfv
    em(ist,jst)=em(ist,jst)+zdotc(ngkq,evecfv2(:,ist),1,wfir,1)
  end do
end do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
if (tevecsv) then
! get the second-variational eigenvectors
  call getevecsv(vkl(:,ik),evecsv1)
  call getevecsv(vkql,evecsv2)
  do i=1,nstsv
    do j=1,nstsv
      zsum=0.d0
      k=0
      do ispn=1,nspinor
        do ist=1,nstfv
          k=k+1
          l=(ispn-1)*nstfv
          do jst=1,nstfv
            l=l+1
            zsum=zsum+em(ist,jst)*conjg(evecsv2(k,i))*evecsv1(l,j)
          end do
        end do
      end do
      emat(i,j)=zsum
    end do
  end do
else
  emat(:,:)=em(:,:)
end if
deallocate(igkqig,gnt,jlqr,vgkql,vgkqc,gkqc,tpgkqc)
deallocate(sfacgkq,apwalm1,apwalm2,evecfv1,evecfv2)
if (tevecsv) deallocate(evecsv1,evecsv2)
deallocate(wfmt1,wfmt2,wfmt3,wfir,zfir1,zfir2,em)
return
end subroutine

