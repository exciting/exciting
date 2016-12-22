! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: findgqmap
! !INTERFACE:
subroutine findgqmap(iq, iqr, nsc, sc, ivgsc, n, isc, isci, ivgu, igqmap)
! !USES:
  use modmpi, only: terminate
  use modinput, only: input
  use mod_qpoint, only: iqmap
  use mod_lattice, only: bvec
  use modxs, only: vqlr, scimap, igqig, ivgigq, ivqr
  use mod_symmetry, only: lsplsymc, symlat, maxsymcrys
  use mod_Gvector, only: ivg
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   integer(4), iq : Index of non-reduced q-point
!   integer(4), iqr : Index of reduced q-point
!   integer(4), nsc : Number of symmetry operations that map iqr into iq
!   integer(4), sc(maxsymcrys) : The first nsc elements are the indices of the 
!                                corresponding crystal symmetries 
!   integer(4), ivgsc(3, maxsymcrys) : The G vectors that shift the rotated vqlr to 
!                                      the unit cell
!   integer(4), n   : Number of G+q vectors
! OUT:
!   integer(4), isc : Crystal symmety operation that maps qr to q
!   integer(4), isci : Inverse crystal symmetry operation that maps q to qr
!   integer(4), ivgu(3) : Contains -G_s of S^T_{isc}*qr = q + G_s
!   integer(4), igqmap(n) : Index mapping between the set of G+iq vectors and 
!                           and the G+iqrnr vectors
!
! !DESCRIPTION:
!   Given an non-reduced and a reduce q point index (iq, iqr) the routine 
!   tries to find a symmetry operation that mapps q onto qr and also
!   keeps all rotated $\vec{G}+\vec{q}$ vectors 
!   $\mathbf{S}*(\vec{G}+\vec{q}) = \vec{G}_1+\vec{q}_\text{r}$ shorter
!   than {\tt gqmax}.
!   If successful the routine returns an index map between the set of
!   $\left\{ \vec{G}+\vec{q} \right\}$ and the
!   $\left\{ \vec{G}+\vec{q}_r \right\}$ vectors.
!
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich 2016)
!   Changed formatting and added comments. (Aurich 2016)
!EOP
!BOC

  implicit none

  ! I/O
  integer(4), intent(in) :: iq, iqr 
  integer(4), intent(in) :: nsc 
  integer(4), intent(in) :: sc(maxsymcrys)
  integer(4), intent(in) :: ivgsc(3, maxsymcrys)
  integer(4), intent(in) :: n

  integer(4), intent(out) :: isc
  integer(4), intent(out) :: isci
  integer(4), intent(out) :: ivgu(3)
  integer(4), intent(out) :: igqmap(n)

  ! Local variables
  real(8) :: vqr(3), v2(3), t1
  integer(4) :: iqrnr, j, isym, isymi, lspl, lspli, iv(3), ivg1(3), igq1

  ! Find map from g-vectors to rotated g-vectors

  ! Get non-reduced 1d index form non-reduced 3d index of a reduced q point
  iqrnr = iqmap(ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))

  ! Get lattice coordinates of reduced q point
  vqr(:) = vqlr(:, iqr)

  ! Loop over symmetries that transform qr into q
  do j = 1, nsc

    ! Get crystal symmetry index
    isym = sc(j)
    ! Index of spatial rotation element corresponding to 
    ! crystal symmetry isym
    lspl = lsplsymc(isym)

    ! Get the corresponding symmetry operation that transforms q into qr
    ! Get index of inverse crystal symmetry
    isymi = scimap(isym)
    ! Index of spatial rotation element corresponding to 
    ! inverse crystal symmetry isym
    lspli = lsplsymc(isymi)

    ! Loop over G+q index
    do igq1 = 1, n

      ! Get integer lattice coordinates of the G corresponding to G+q
      ivg1(:) = ivg(:, igqig(igq1, iq))

      ! Note: The S are real space rotation matrices, which
      !       act with their transpose in k-space.
      ! Note: Rotation matrices are orthogonal S*S^T=1
      ! 
      ! S^T*qr = q + G_s  <--> qr = S*q + S*G_s
      ! S*(q+G) = qr + S*(G-G_s) = qr + G1
      !
      ! Note: findsymegqiv did not save G_s but -G_s in ivgsc

      ! G1 = Si^-1 *( G + G_s ) , where Si is the inverse of S
      iv = matmul(transpose(symlat(:, :, lspli)), ivg1+ivgsc(:, j))

      ! |g1 + qr|
      v2 = matmul(bvec, iv+vqr)
      t1 = sqrt(sum(v2**2))

      if( n .gt. 1 .and. t1 .gt. input%xs%gqmax ) then
        write(*,*)
        write(*, '("Info(findgqmap): Need one more symmetry operation")')
        write(*,*)
        go to 10
      end if

      ! Locate G1+qr in G+q-vector set
      igqmap(igq1) = ivgigq(iv(1), iv(2), iv(3), iqrnr)

      if(igqmap(igq1) .le. 0) then
        write(*,*)
        write(*, '("Error(findgqmap): Failed to map rotated g-vector")')
        write(*, '(" non-reduced q-point		       :", i8)') iq
        write(*, '(" reduced q-point 		       :", i8)') iqr
        write(*, '(" reduced q-point in non-reduced set     :", i8)') iqrnr
        write(*, '(" G+q-vector index (non-reduced q-point) :", i8)') igq1
        write(*, '(" rotated G-vector		       :", 3i8)') iv
        write(*,*)
        call terminate
      end if

    ! End loop over G+q
    end do

    ! Store G shift vector of the chosen 
    ! q reducing symmetry operation
    ivgu(:) = ivgsc(:, j)

    ! Save the symmetry operation index
    isc = isym

    ! Save the inverse symmetry operation index
    isci = isymi

    go to 20

10  continue

  ! end loop over symmetry operations
  end do

  write(*,*)
  write(*, '("Error(findgqmap): Failed to reduce q-point: ", i8)') iq
  write(*,*)
  call terminate

20    continue
end subroutine findgqmap
!EOC
