

! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsymi
! !INTERFACE:


subroutine findsymi(epslat, maxsymcrys, nsymcrys, symlat, lsplsymc, vtlsymc, isymlat, &
     scimap)
! !DESCRIPTION:
!   Throughout the code the symmetries are understood to be applied in a way 
!   $$ (\alpha_S|\alpha_R|{\bf t}) {\bf x} = \alpha_S\alpha_R
!   ({\bf x}+{\bf t})$$
!   which is different from the commonly used definition
!   $\{\alpha|\tau\}x=\alpha x+\tau$ -- see routine {\tt findsymcrys}.
!   This difference affects the inverse of the fractional translation
!   but has no effect on the inverse of the rotational part, so the inverse
!   spacegroup symmetry operations are the same for both definitions.
!
! !REVISION HISTORY:
!   Created April 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: epslat
  integer, intent(in) :: symlat(3, 3, 48)
  integer, intent(in) :: maxsymcrys, nsymcrys
  integer, intent(in) :: lsplsymc(nsymcrys)
  real(8), intent(in) :: vtlsymc(3, maxsymcrys)
  integer, intent(in) :: isymlat(48)
  integer, intent(out) :: scimap(maxsymcrys)
  ! local variables
  real(8) :: c(3, 3), si(3, 3), sj(3, 3), vtl(3)
  integer :: i, isym, jsym, lspli, lsplj, iv(3)
  scimap(:)=0
  do isym=1, nsymcrys
     lspli=lsplsymc(isym)
     si(:, :)=dble(symlat(:, :, lspli))
     do jsym=1, nsymcrys
	lsplj=lsplsymc(jsym)
	sj(:, :)=dble(symlat(:, :, lsplj))
        ! translation
	vtl(:)=vtlsymc(:, jsym)
	vtl=matmul(sj, vtl)
	vtl(:)=vtl(:)+vtlsymc(:, isym)
	call r3frac(epslat, vtl, iv)
        ! rotation
	call r3mm(si, sj, c)
        ! subract unit matrix
	forall (i=1:3)
	   c(i, i)=c(i, i)-1.d0
	end forall
	if ((sum(vtl).lt.epslat).and.(sum(abs(c)).lt.epslat)) then
           ! isym is inverse of jsym
	   scimap(isym)=jsym
	   goto 10
	end if	      
     end do
10   continue
     ! check if inverse symmetry is consistent with inverse lattice symmetry
     if (isymlat(lsplsymc(isym)).ne.lsplsymc(jsym)) then
	write(*, *)
	write(*, '("Error(findsymi): inconsistency with inverse lattice &
	     &symmetry")')
	write(*, '(" space group symmetry:", t40, i6)') isym
	write(*, '(" lattice symmetry:", t40, i6)') lsplsymc(isym)
	write(*, '(" inverse lattice symmetry:", t40, i6)') isymlat(lsplsymc(isym))
	write(*, '(" proposed inverse lattice symmetry:", t40, i6)') lsplsymc(jsym)
	write(*, *)
	call terminate
     end if
  end do
end subroutine findsymi
!EOC
