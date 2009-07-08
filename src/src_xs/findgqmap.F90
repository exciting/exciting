

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine findgqmap(iq, iqr, nsc, sc, ivgsc, n, isc, isci, ivgu, igqmap)
  use modmain
use modinput
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iq, iqr, nsc, sc(maxsymcrys), ivgsc(3, maxsymcrys), n
  integer, intent(out) :: isc, isci, ivgu(3), igqmap(n)
  ! local variables
  real(8) :: vqr(3), v2(3), t1
  integer :: iqrnr, j, isym, isymi, lspl, lspli, iv(3), ivg1(3), igq1
  ! find map from G-vectors to rotated G-vectors
  iqrnr=iqmap(ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))
  vqr(:)=vqlr(:, iqr)
  do j=1, nsc
     isym=sc(j)
     lspl=lsplsymc(isym)
     isymi=scimap(isym)
     lspli=lsplsymc(isymi)
     do igq1=1, n
	ivg1(:)=ivg(:, igqig(igq1, iq))
        ! G1 = si^-1 * ( G + G_s ) , where si is the inverse of s
	iv=matmul(transpose(symlat(:, :, lspli)), ivg1+ivgsc(:, j))
        ! |G1 + q|
	v2=matmul(bvec, iv+vqr)
	t1=sqrt(sum(v2**2))
	if ((n.gt.1).and.(t1.gt.input%xs%gqmax)) then
	   write(*, *)
	   write(*, '("Info(findgqmap): need one more symmetry operation")')
	   write(*, *)
	   goto 10
	end if
        ! locate G1 + q in G+q-vector set
	igqmap(igq1)=ivgigq(iv(1), iv(2), iv(3), iqrnr)
	if (igqmap(igq1).le.0) then
	   write(*, *)
	   write(*, '("Error(findgqmap): failed to map rotated G-vector")')
	   write(*, '(" non-reduced q-point		       :", i8)') iq
	   write(*, '(" reduced q-point 		       :", i8)') iqr
	   write(*, '(" reduced q-point in non-reduced set     :", i8)') iqrnr
	   write(*, '(" G+q-vector index (non-reduced q-point) :", i8)') igq1
	   write(*, '(" rotated G-vector		       :", 3i8)') iv
	   write(*, *)
	   call terminate
	end if
        ! end loop over G
     end do
     ! store G1 vector
     ivgu(:)=ivgsc(:, j)
     isc=isym
     isci=isymi
     goto 20
10   continue
     ! end loop over symmetry operations
  end do
  write(*, *)
  write(*, '("Error(findgqmap): failed to reduce q-point: ", i8)') iq
  write(*, *)
  call terminate
20 continue
end subroutine findgqmap
