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
integer lspl,i,j,k,l,m,n
real(8) s(3,3),v(3),t1
real(8) d1(3,3),d2(3,3)
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
s(:,:)=dble(symlat(:,:,lspl))
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! equivalent atom with this symmetry
    ja=ieqatom(ia,is,isym)
    jas=idxas(ja,is)
    map(jas)=ias
    call r3mv(s,atposl(1,ja,is),v)
    v(:)=atposl(:,ia,is)-v(:)
    t1=twopi*(vpl(1)*v(1)+vpl(2)*v(2)+vpl(3)*v(3))
! phase factor
    zph(ias)=cmplx(cos(t1),sin(t1),8)
  end do
end do
! convert symmetry to Cartesian coordinates
call r3mm(s,ainv,s)
call r3mm(avec,s,s)
! rotate and phase-shift dynamical matrix with symmetry
do ias=1,natmtot
  i=3*(ias-1)
  k=3*(map(ias)-1)
  do jas=1,natmtot
    j=3*(jas-1)
    l=3*(map(jas)-1)
    do m=1,3
      do n=1,3
        d1(m,n)=dble(dyn(i+m,j+n))
        d2(m,n)=aimag(dyn(i+m,j+n))
      end do
    end do
    call r3mm(d1,s,d1)
    call r3mtm(s,d1,d1)
    call r3mm(d2,s,d2)
    call r3mtm(s,d2,d2)
    zt1=zph(ias)*conjg(zph(jas))
    do m=1,3
      do n=1,3
        dyns(k+m,l+n)=dyns(k+m,l+n)+zt1*cmplx(d1(m,n),d2(m,n),8)
      end do
    end do
  end do
end do
return
end subroutine

