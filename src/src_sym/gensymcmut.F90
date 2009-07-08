

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP


subroutine gensymcmut(eps, maxsymcrys, nsymcrys, symlat, lsplsymc, vtlsymc, scmut, &
     tabel)
! !DESCRIPTION:
!   Sets up the group multiplication table. The table is checked for consistency
!   in a way that it is required that every elements occurrs once and only once
!   in each row and column of the table. The first row and colmuns must consist
!   of the indentity since the first symmetry element is the identity by
!   convention.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: eps
  integer, intent(in) :: maxsymcrys, nsymcrys
  integer, intent(in) :: symlat(3, 3, 48)
  integer, intent(in) :: lsplsymc(nsymcrys)
  real(8), intent(in) :: vtlsymc(3, maxsymcrys)
  integer, intent(out) :: scmut(nsymcrys, nsymcrys)
  logical, intent(out) :: tabel
  ! local variables
  integer, parameter :: maxsymlat=48
  integer :: isym, jsym, asym, lspli, lsplj, lspla, iv(3)
  integer :: doner(maxsymlat), donec(maxsymlat)
  real(8) :: c(3, 3), ct(3, 3), si(3, 3), sj(3, 3), sa(3, 3), vtt(3), vtl(3), vtla(3)
  scmut(:, :)=0
  do isym=1, nsymcrys
     lspli=lsplsymc(isym)
     si(:, :)=dble(symlat(:, :, lspli))
     do jsym=1, nsymcrys
	lsplj=lsplsymc(jsym)
	sj(:, :)=dble(symlat(:, :, lsplj))
        ! calculate rotation
	c=matmul(si, sj)
        ! first translation
	vtl(:)=matmul(sj, vtlsymc(:, jsym))
        ! second translation
	vtl(:)=vtl(:)+vtlsymc(:, isym)
	call r3frac(eps, vtl, iv)
	vtl(:)=matmul(si, vtl)
	call r3frac(eps, vtl, iv)
	do asym=1, nsymcrys
	   lspla=lsplsymc(asym)
	   sa(:, :)=dble(symlat(:, :, lspla))
           ! third translation
	   vtla(:)=matmul(sa, vtlsymc(:, asym))
	   call r3frac(eps, vtla, iv)
           ! differece in rotation
	   ct(:, :)=c(:, :)-sa(:, :)
           ! difference in translation
	   vtt(:)=vtl(:)-vtla(:)
	   if ((sum(abs(ct)).lt.eps).and.(sum(abs(vtt)).lt.eps)) then
              ! add element to multiplication table
	      scmut(isym, jsym)=asym
	      cycle
	   end if
	end do
     end do
  end do
  ! check multiplication table for consistency
  do isym=1, nsymcrys
     donec(:)=0
     doner(:)=0
     do jsym=1, nsymcrys
	doner(scmut(isym, jsym))=doner(scmut(isym, jsym))+1
	donec(scmut(jsym, isym))=donec(scmut(jsym, isym))+1
     end do
     do jsym=1, nsymcrys
	if (doner(jsym).ne.1) then  
	   write(*, *)
	   write(*, '("Error(gensymcmut): error in multiplication table in &
		&row")')
	   write(*, '(" row number    : ", i6)') isym
	   write(*, '(" column number : ", i6)') jsym
	   write(*, '(" multiple occurrence : ", i6)') doner(jsym)
	   write(*, *)
	   stop
	end if
	if (donec(jsym).ne.1) then  
	   write(*, *)
	   write(*, '("Error(gensymcmut): error in multiplication table in &
		&column")')
	   write(*, '(" row number    : ", i6)') jsym
	   write(*, '(" column number : ", i6)') isym
	   write(*, '(" multiple occurrence : ", i6)') donec(jsym)
	   write(*, *)
	   stop
	end if
     end do
     ! check first row and column
     if ((scmut(1, isym).ne.isym).or.(scmut(isym, 1).ne.isym)) then
	write(*, *)
	write(*, '("Error(gensymcmut): error in multiplication table")')
	write(*, '(" first row or column have wrong property")')
	write(*, '(" position : ", i6)') isym
	write(*, *)
	stop
     end if
  end do
  ! check if group is Abelian
  tabel=.false.
  if (all((scmut-transpose(scmut)).eq.0)) tabel=.true.
end subroutine gensymcmut
!EOC
