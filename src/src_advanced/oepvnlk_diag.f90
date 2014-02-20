
subroutine oepvnlk_diag(ikp,vnlvv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vnlvv(nstsv)
! local variables
integer ngknr,ik,ist2,ist3
integer lmax
integer iq,ig,iv(3),igq0
real(8) v(3),cfq
complex(8) zrho01,zrho02,zt1,zt2
! automatic arrays
real(8) zn(nspecies)
complex(8) sfacgq0(natmtot)
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:)
real(8), allocatable :: tpgqc(:,:)
real(8), allocatable :: gqc(:)
real(8), allocatable :: jlgqr(:,:,:)
real(8), allocatable :: jlgq0r(:,:,:)
real(8), allocatable :: evalsvp(:)
real(8), allocatable :: evalsvnr(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: ylmgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: wfcr1(:,:,:)
complex(8), allocatable :: wfcr2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
complex(8), allocatable :: zvcltp(:,:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax))
allocate(vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax))
allocate(tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec))
allocate(tpgqc(2,ngvec))
allocate(gqc(ngvec))
allocate(jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1,ngvec,nspecies))
allocate(jlgq0r(0:input%groundstate%lmaxvr,nrcmtmax,nspecies))
allocate(evalsvp(nstsv))
allocate(evalsvnr(nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec))
allocate(sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(wfcr1(lmmaxvr,nrcmtmax,2))
allocate(wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
allocate(zvcltp(lmmaxvr,nrcmtmax))
allocate(zfmt(lmmaxvr,nrcmtmax))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! set the nuclear charges to zero
zn(:)=0.d0
vnlvv(:)=0.d0
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(:,ikp),evalsvp)
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for all states for the input k-point
call genwfsv(.false.,ngk(1,ikp),igkig(:,1,ikp),evalsvp,apwalm,evecfv,evecsv, &
 wfmt1,wfir1)
! start loop over non-reduced k-point set
do ik=1,nkptnr
! generate G+k-vectors
  call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
   tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-points
  call getevalsv(vklnr(:,ik),evalsvnr)
  call getevecfv(vklnr(:,ik),vgklnr,evecfv)
  call getevecsv(vklnr(:,ik),evecsv)
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! determine q-vector
  iv(:)=ivk(:,ikp)-ivknr(:,ik)
  iv(:)=modulo(iv(:),input%groundstate%ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkcnr(:,ik)
  do ig=1,ngvec
! determine G+q vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vector
    call genylm(input%groundstate%lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  lmax=input%groundstate%lmaxvr+input%groundstate%npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
  call genjlgq0r(gqc(igq0),jlgq0r)
! calculate the wavefunctions for occupied states
  call genwfsv(.true.,ngknr,igkignr,evalsvnr,apwalm,evecfv,evecsv,wfmt2,wfir2)
! internal sum over the occupied states (k' in formula (2))
  do ist3=1,nstsv
    if (evalsvnr(ist3).lt.efermi) then
      do ist2=1,nstsv
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist2), &
           wfir2(:,:,ist3),wfir1(:,:,ist2),zrhomt,zrhoir)
! calculate the Coulomb potential
          call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq, &
           sfacgq,zn,zrhomt,zrhoir,zvclmt,zvclir,zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist2), &
           wfir2(:,:,ist3),wfir1(:,:,ist2),zrhomt,zrhoir)
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
! compute the density coefficient of the smallest G+q-vector
          call zrhogp(gqc(igq0),jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt, &
           zrhoir,zrho01)
          zt2=cfq*wiq2(iq)*(conjg(zrho01)*zrho02)
          vnlvv(ist2)=vnlvv(ist2)-(wkptnr(ik)*zt1+zt2)
! end loop over ist2
      end do
! end loop over ist3
    end if
  end do
! end loop over non-reduced k-point set
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(evalsvp,evalsvnr,evecfv,evecsv)
deallocate(apwalm,sfacgknr,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr1,wfcr2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zvcltp,zfmt)
return
end subroutine
!EOC

