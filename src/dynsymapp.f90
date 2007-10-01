subroutine dynsymapp(sym,vpl,dyn,tapp,dyns)
use modmain
implicit none
! arguments
integer, intent(in) :: sym(3,3)
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: dyn(3*natmtot,3*natmtot)
logical, intent(out) :: tapp
complex(8), intent(inout) :: dyns(3*natmtot,3*natmtot)
! local variables
integer is,ia1,ia2,ias1,ias2
integer isym,i,j,k,l,m,n
real(8) s(3,3),d1(3,3),d2(3,3),v(3),t1
complex(8) zt1
! automatic arrays
integer map(natmtot)
complex(8) zp(natmtot)
tapp=.true.
! check if symmetry is the identity
if ((sym(1,1).eq.1).and.(sym(1,2).eq.0).and.(sym(1,3).eq.0).and. &
    (sym(2,1).eq.0).and.(sym(2,2).eq.1).and.(sym(2,3).eq.0).and. &
    (sym(3,1).eq.0).and.(sym(3,2).eq.0).and.(sym(3,3).eq.1)) then
  dyns(:,:)=dyns(:,:)+dyn(:,:)
  return
end if
s(:,:)=dble(sym(:,:))
! find map from each atom to its equivalent under the symmetry
do is=1,nspecies
  do ia1=1,natoms(is)
    ias1=idxas(ia1,is)
    do ia2=1,natoms(is)
      ias2=idxas(ia2,is)
      do i=1,nsymeqat(ia2,ia1,is)
        isym=symeqat(i,ia2,ia1,is)
        do k=1,3
          do l=1,3
            if (sym(k,l).ne.symlat(k,l,isym)) goto 10
          end do
        end do
        goto 20
10 continue
      end do
    end do
! no equivalent atom found so symmetry cannot be applied to dynamical matrix
    tapp=.false.
    return
20 continue
! S maps atom 2 into atom 1
    map(ias2)=ias1
    call r3mv(s,atposl(1,ia2,is),v)
    v(:)=atposl(:,ia1,is)-v(:)
    t1=twopi*(vpl(1)*v(1)+vpl(2)*v(2)+vpl(3)*v(3))
    zp(ias1)=cmplx(cos(t1),sin(t1),8)
  end do
end do
! convert symmetry to Cartesian coordinates
call r3mm(s,ainv,s)
call r3mm(avec,s,s)
! rotate and phase-shift dynamical matrix with symmetry
do ias1=1,natmtot
  i=3*(ias1-1)
  k=3*(map(ias1)-1)
  do ias2=1,natmtot
    j=3*(ias2-1)
    l=3*(map(ias2)-1)
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
    zt1=zp(ias1)*conjg(zp(ias2))
    do m=1,3
      do n=1,3
        dyns(k+m,l+n)=dyns(k+m,l+n)+zt1*cmplx(d1(m,n),d2(m,n),8)
      end do
    end do
  end do
end do
return
end subroutine

