
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gensymdf
  implicit none
contains

  subroutine gensymdf(oc1,oc2)
    use modmain
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: oc1,oc2
    ! local variables
    integer :: isym,i,j,lspl
    real(8) :: sc(3,3),s(3,3),d(nsymcrys)

    ! pre-calculation for symmetrisation
    s(:,:)=0.d0
    do isym=1,nsymcrys
       lspl=lsplsymc(isym)
       sc(:,:)=dble(symlat(:,:,lspl))
       call r3mtm(sc,binv,sc)
       call r3mm(bvec,sc,sc)
       d(isym)=sc(1,1)*sc(2,2)-sc(1,2)*sc(2,1)
       do i=1,3
          do j=1,3
             s(i,j)=s(i,j) + sc(i,oc1)*sc(j,oc2)
          end do
       end do
    end do
    symdfq0(:,:)=s(:,:)/dble(nsymcrys)
    tsymdfq0dn=.false.
    ! occurrance of negative 3,3-minor
    if (any(d.lt.0.d0)) tsymdfq0dn=.true.

! /// old stuff  <= version 0.9.74
!!$    ! dielectric tensor components
!!$    s(:,:)=0.d0
!!$    do isym=1,nsymcrys
!!$       a(:,:)=dble(symcrys(:,:,isym))
!!$       call r3mtm(a,binv,a)
!!$       call r3mm(bvec,a,a)
!!$       do i=1,3
!!$          do j=1,3
!!$             s(i,j)=s(i,j)+a(i,oc1)*a(j,oc2)
!!$          end do
!!$       end do
!!$    end do
!!$    s(:,:)=s(:,:)/dble(nsymcrys)
!!$    symdfq0(:,:)=s(:,:)


!!$    ! trivial symmetry for the moment
!!$    symdfq0(:,:)=0.d0
!!$    if (nosym) then
!!$       symdfq0(optcomp(1,1),optcomp(1,1))=1.d0
!!$    else
!!$       forall (i=1:3)
!!$          symdfq0(i,i)=1.d0/3.d0
!!$       end forall
!!$    end if

  end subroutine gensymdf

end module m_gensymdf
