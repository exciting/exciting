! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: findsymeqiv
! !INTERFACE:
subroutine findsymeqiv(tfbz, vpl, vplr, nsc, sc, ivgsc)
! !USES:
  use modmpi, only: terminate
  use modinput, only: input
  use mod_lattice, only: bvec
  use mod_symmetry, only: nsymcrys, lsplsymc, symlat, maxsymcrys
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   logical, tfbz    :  Use 1st Bz of [0,1) as unit cell
!   real(8), vpl(3)  :  Lattice coordinates of non-reduced k point
!   real(8), vplr(3) :  Lattice coordinates of reduced k point
! OUT:
!   integer(4), nsc  :  Number of symmetry operations that transform vplr into vpl
!   integer(4), sc(maxsymcrys) : The first nsc entries contain the 
!                                corresponding crystal symmetry indices
!   integer(4), ivgsc(3, maxsymcrys) : G vectors that shift the rotated vplr
!                                      back to the unit cell
!
! !DESCRIPTION:
!   Given one non-reduced k-point vector and a reduces one, the routine
!   determines the symmetry operations that rotate the reduced k point to
!   the non-reduced one.
!   Note: The routine terminates the program, if no symmetry operation is found.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich 2016)
!   Changed formatting and added comments. (Aurich 2016)
!
!EOP
!BOC

  implicit none

  ! I/O
  logical, intent(in) :: tfbz 
  real(8), intent(in) :: vpl(3)
  real(8), intent(in) :: vplr(3) 

  integer(4), intent(out) :: nsc
  integer(4), intent(out) :: sc(maxsymcrys) 
  integer(4), intent(out) :: ivgsc(3, maxsymcrys) 

  ! Local variables
  integer(4) :: isym, lspl, iv(3)
  real(8) :: s(3, 3), v1(3), t1

  ! External functions
  real(8), external :: r3taxi

  ! Symmetries that transform non-reduced q-point to reduced one, namely
  ! vpl = s^-1 * vplr + G_s
  nsc = 0
  ! Loop over all crystal symmetries
  do isym = 1, nsymcrys

    ! Index of spatial rotation element corresponding to 
    ! crystal symmetry isym
    lspl = lsplsymc(isym)

    ! Get the rotation matrix
    s(:, :) = dble(symlat(:, :, lspl))

    ! r-space rotation matrices act on k-space vectros with thier transposed
    ! Rotate reduced k-vector: s^T*k = k'
    call r3mtv(s, vplr, v1)

    if(tfbz) then
      ! Map back to 1st BZ
      call vecfbz(input%structure%epslat, bvec, v1, iv)
    else
      ! Map back to [0,1)
      call r3frac(input%structure%epslat, v1, iv)
    end if

    ! Check if transformed reduced k-vector is equal to non-reduced k-vector
    t1 = r3taxi(vpl, v1)
    if(t1 .lt. input%structure%epslat) then
       nsc = nsc + 1
       ! Save crystal symmetry index
       sc(nsc) = isym
       ! Save shift vector that maps the rotated k back to unit cell
       ivgsc(:, nsc) = -iv(:)
    end if

  end do

  if(nsc .eq. 0) then
    write(*,*)
    write(*, '("Error(findsymeqiv): p-points are not equivalent by symmetry")')
    write(*, '(" vpl  :", 3g18.10)') vpl
    write(*, '(" vplr :", 3g18.10)') vplr
    write(*,*)
    call terminate
  end if

end subroutine findsymeqiv
!EOC
