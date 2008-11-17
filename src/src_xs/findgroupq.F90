
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findgroupq
! !INTERFACE:
subroutine findgroupq(tfbz,vql,epslat,bvec,symlat,nsymcrys,lsplsymc,nsymcrysq, &
     scqmap,ivscwrapq)
! !DESCRIPTION:
!   Find the (little) group of {\bf q} (which includes finding the small group
!   of {\bf q}).
!   All symmetries, where the rotational part transforms {\bf q} into an
!   equivalent vector are collected for the small group of {\bf q}. Inclusion
!   of non-primitive translations yields the little group of {\bf q}.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  logical, intent(in) :: tfbz
  real(8), intent(in) :: vql(3)
  real(8), intent(in) :: epslat
  real(8), intent(in) :: bvec(3,3)
  integer, intent(in) :: symlat(3,3,48)
  integer, intent(in) :: nsymcrys
  integer, intent(in) :: lsplsymc(nsymcrys)
  integer, intent(out) :: nsymcrysq
  integer, intent(out) :: scqmap(nsymcrys)
  integer, intent(out) :: ivscwrapq(3,nsymcrys)
  ! local variables
  character(*), parameter :: thisnam = 'findgroupq'
  integer :: isym,lspl,iv(3),ivt(3)
  real(8) :: s(3,3),v1(3),v1t(3),v2(3),t1
  real(8), external :: r3taxi
  nsymcrysq=0
  ivscwrapq(:,:)=0
  ! loop over space group elements
  do isym=1,nsymcrys
     ! get lattice point group element
     lspl=lsplsymc(isym)
     ! rotation as real matrix in lattice coordinates
     s(:,:)=dble(symlat(:,:,lspl))
     ! Note: here we apply the rotation from the left side
     call r3mtv(s,vql,v1)
     ! save transformed vector
     v1t(:)=v1(:)
     ! convert v1 to equivalent point and wrapping vector
     call r3frac(epslat,v1,iv)
     ! mapt to first Brillouin zone
     if (tfbz) then
        call vecfbz(epslat,bvec,v1,ivt)
        iv(:)=ivt(:)-iv(:)
     end if
     iv=-iv
     ! check if new vector is equal to orinial vql vector
     t1=r3taxi(vql,v1)
     if (t1.lt.epslat) then
        ! check again if qL = sL^T qL + GL ( q = s^-1 q + G )
        v2(:)=vql(:)-dble(iv(:))
        t1=r3taxi(v1t,v2)
        ! +++ should be obsolescent if r3taxi is working properly +++
        if (t1.gt.epslat) then
           write(*,'(a)') 'Error('//thisnam//'): inconsistency in &
                &wrapping vector G from q1=q+G: q/q1/G:'
           write(*,'(a,3g18.10)') ' v2:',v2
           write(*,'(a,3g18.10)') ' q1:',v1t
           write(*,'(a,3g18.10)') ' q :',vql
           write(*,'(a,3i9)')     ' G :',iv
           call terminate
        end if
        ! rotation is in small group of q (G0(q))
        nsymcrysq=nsymcrysq+1
        ! map from little group of q (G(q)) to space group (G)
        scqmap(nsymcrysq)=isym
        ! wrapping vector (reciprocal lattice vector)
        ivscwrapq(:,isym)=iv(:)
     end if
  end do
  if (nsymcrysq.lt.1) then
     write(*,'(a,3g18.10)') 'Error('//thisnam//'): empty little group of q &
          &for q:', vql
     call terminate
  end if
end subroutine findgroupq
!EOC
