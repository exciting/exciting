
! Copyright (C) 2005-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynsymapp(isym,vpl,dyn,dyns)
use modmain
implicit none
! arguments
integer, intent(in) :: isym
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: dyn(3*natmtot,3*natmtot)
complex(8), intent(inout) :: dyns(3*natmtot,3*natmtot)
! local variables
integer is,ia,ja,ias,jas
integer lspl,i,j,k,l,m,n,iv(3)
real(8) s(3,3),a(3,3),b(3,3),c(3,3)
real(8) v1(3),v2(3),v3(3),t1
complex(8) zt1
! automatic arrays
integer map(natmtot)
complex(8) zph(natmtot)
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! check if symmetry is the identity
if (lspl.eq.1) then
  dyns(:,:)=dyns(:,:)+dyn(:,:)
  return
end if
! symmetry in lattice coordinates
s(:,:)=dble(symlat(:,:,lspl))
! map vpl to the first Brillouin zone
v1(:)=vpl(:)
call vecfbz(epslat,bvec,v1,iv)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! equivalent atom with this symmetry
    ja=ieqatom(ia,is,isym)
    jas=idxas(ja,is)
    map(ias)=jas
! phase factor
    v2(:)=atposl(:,ja,is)+vtlsymc(:,isym)
    call r3mv(s,v2,v3)
    v3(:)=atposl(:,ia,is)-v3(:)
    t1=twopi*(v1(1)*v3(1)+v1(2)*v3(2)+v1(3)*v3(3))
    zph(ias)=cmplx(cos(t1),sin(t1),8)
  end do
end do
! symmetry in Cartesian coordinates
s(:,:)=symlatc(:,:,lspl)
! rotate and phase-shift dynamical matrix with symmetry
do ias=1,natmtot
  i=3*(ias-1)
  k=3*(map(ias)-1)
  do jas=1,natmtot
    j=3*(jas-1)
    l=3*(map(jas)-1)
    do m=1,3
      do n=1,3
        a(m,n)=dble(dyn(i+m,j+n))
        b(m,n)=aimag(dyn(i+m,j+n))
      end do
    end do
    call r3mtm(s,a,c)
    call r3mm(c,s,a)
    call r3mtm(s,b,c)
    call r3mm(c,s,b)
    zt1=zph(ias)*conjg(zph(jas))
    do m=1,3
      do n=1,3
        dyns(k+m,l+n)=dyns(k+m,l+n)+zt1*cmplx(a(m,n),b(m,n),8)
      end do
    end do
  end do
end do
return
end subroutine

