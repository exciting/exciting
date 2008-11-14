
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getevecfv
! !INTERFACE:
subroutine getevecfv(vpl,vgpl,evecfv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   vpl    : p-point vector in lattice coordinates (in,real(3))
!   vgpl   : G+p-vectors in lattice coordinates (out,real(3,ngkmax))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv,nspnfv))
! !DESCRIPTION:
!   Reads in a first-variational eigenvector from file. If the input $k$-point,
!   ${\bf p}$, is not in the reduced set, then the eigenvector of the equivalent
!   point is read in and the required rotation/translation operations applied.
!
! !REVISION HISTORY:
!   Created Feburary 2007 (JKD)
!   Fixed transformation error, October 2007 (JKD, Anton Kozhevnikov)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer isym,lspl,ilspl
integer ilo,l,m,lm,i
integer ik,igp,igk,ig
integer is,ia,ja,ias,jas
integer recl,nmatmax_,nstfv_,nspnfv_
real(8) vkl_(3),v(3),t1
real(8) si(3,3),sc(3,3)
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: evecfvt(:,:)
complex(8), allocatable :: zflm1(:,:),zflm2(:,:)
! find the equivalent k-point number and crystal symmetry element
call findkpt(vpl,isym,ik)
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! find the record length
inquire(iolength=recl) vkl_,nmatmax_,nstfv_,nspnfv_,evecfv
!$OMP CRITICAL
open(70,file=trim(scrpath)//'EVECFV'//trim(filext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nmatmax_,nstfv_,nspnfv_,evecfv
close(70)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECFV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nmatmax.ne.nmatmax_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nmatmax for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nmatmax
  write(*,'(" EVECFV.OUT : ",I8)') nmatmax_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVECFV.OUT : ",I8)') nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVECFV.OUT : ",I8)') nspnfv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) return
if (spinsprl) then
  write(*,*)
  write(*,'("Error(getevec): code limitation - cannot rotate spin-spiral &
  &states")')
  write(*,'(" (first run one self-consistent loop with no k-point reduction)")')
  write(*,*)
  stop
end if
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:,:)=symlat(:,:,ilspl)
!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
allocate(evecfvt(nmatmax,nstfv))
do igk=1,ngk(1,ik)
  ig=igkig(igk,1,ik)
  v(:)=dble(ivg(:,ig))
  t1=-twopi*dot_product(v(:),vtlsymc(:,isym))
  zt1=cmplx(cos(t1),sin(t1),8)
  evecfvt(igk,:)=zt1*evecfv(igk,:,1)
end do
do igk=1,ngk(1,ik)
  call r3mtv(si,vgkl(:,igk,1,ik),v)
  do igp=1,ngk(1,ik)
    t1=abs(v(1)-vgpl(1,igp))+abs(v(2)-vgpl(2,igp))+abs(v(3)-vgpl(3,igp))
    if (t1.lt.epslat) then
      evecfv(igp,:,1)=evecfvt(igk,:)
      goto 10
    end if
  end do
10 continue
end do
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
if (nlotot.le.0) goto 20
allocate(zflm1(lolmmax,nstfv),zflm2(lolmmax,nstfv))
! make a copy of the local-orbital coefficients
do i=ngk(1,ik)+1,nmat(1,ik)
  evecfvt(i,:)=evecfv(i,:,1)
end do
! spatial rotation symmetry matrix in Cartesian coordinates
sc(:,:)=symlatc(:,:,lspl)
! rotate k-point by inverse symmetry matrix
call r3mtv(si,vkl(:,ik),v)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! equivalent atom for this symmetry
    ja=ieqatom(ia,is,isym)
    jas=idxas(ja,is)
! phase factor from translation
    t1=-twopi*dot_product(vkl(:,ik),atposl(:,ja,is))
    zt1=cmplx(cos(t1),sin(t1),8)
    t1=twopi*dot_product(v(:),atposl(:,ia,is))
    zt1=zt1*cmplx(cos(t1),sin(t1),8)
! rotate local orbitals
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      zflm1(:,:)=0.d0
      do m=-l,l
        lm=idxlm(l,m)
        i=ngk(1,ik)+idxlo(lm,ilo,jas)
        zflm1(lm,:)=evecfvt(i,:)
      end do
      call rotzflm(sc,l,nstfv,lolmmax,zflm1,zflm2)
      do m=-l,l
        lm=idxlm(l,m)
        i=ngk(1,ik)+idxlo(lm,ilo,ias)
        evecfv(i,:,1)=zt1*zflm2(lm,:)
      end do
    end do
  end do
end do
deallocate(zflm1,zflm2)
20 continue
deallocate(evecfvt)
return
end subroutine
!EOC

