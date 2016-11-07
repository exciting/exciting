!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: findgroupq
! !INTERFACE:
!
!
Subroutine findgroupq (tfbz, vql, epslat, bvec, symlat, nsymcrys, &
& lsplsymc, nsymcrysq, scqmap, ivscwrapq)
  use modmpi
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
      Implicit None
  ! arguments
      Logical, Intent (In) :: tfbz
      Real (8), Intent (In) :: vql (3)
      Real (8), Intent (In) :: epslat
      Real (8), Intent (In) :: bvec (3, 3)
      Integer, Intent (In) :: symlat (3, 3, 48)
      Integer, Intent (In) :: nsymcrys
      Integer, Intent (In) :: lsplsymc (nsymcrys)
      Integer, Intent (Out) :: nsymcrysq
      Integer, Intent (Out) :: scqmap (nsymcrys)
      Integer, Intent (Out) :: ivscwrapq (3, nsymcrys)
  ! local variables
      Character (*), Parameter :: thisnam = 'findgroupq'
      Integer :: isym, lspl, iv (3)
      Real (8) :: s (3, 3), v1 (3), v1t (3), v2 (3), t1
      Real (8), External :: r3taxi
      nsymcrysq = 0
      ivscwrapq (:, :) = 0
  ! loop over space group elements
      Do isym = 1, nsymcrys
     ! get lattice point group element
         lspl = lsplsymc (isym)
     ! rotation as real matrix in lattice coordinates
         s (:, :) = dble (symlat(:, :, lspl))
     ! Note: here we apply the rotation from the left side
         Call r3mtv (s, vql, v1)
     ! save transformed vector
         v1t (:) = v1 (:)
     ! mapt to first Brillouin zone
         If (tfbz) Then
            Call vecfbz (epslat, bvec, v1, iv)
         Else
        ! convert v1 to equivalent point and wrapping vector
            Call r3frac (epslat, v1, iv)
         End If
         iv = - iv
     ! check if new vector is equal to orinial vql vector
         t1 = r3taxi (vql, v1)
         If (t1 .Lt. epslat) Then
        ! check again if qL = sL^T qL + GL ( q = s^-1 q + G )
            v2 (:) = vql (:) - dble (iv(:))
            t1 = r3taxi (v1t, v2)
        ! +++ should be obsolescent if r3taxi is working properly +++
            If (t1 .Gt. epslat) Then
               Write (*, '(a)') 'Error(' // thisnam // '): inconsistenc&
              &y in wrapping vector G from q1 = q + G: q/q1/G:'
               Write (*, '(a, 3g18.10)') ' v2:', v2
               Write (*, '(a, 3g18.10)') ' q1:', v1t
               Write (*, '(a, 3g18.10)') ' q :', vql
               Write (*, '(a, 3i9)') ' G :', iv
               Call terminate
            End If
        ! rotation is in small group of q (G0(q))
            nsymcrysq = nsymcrysq + 1
        ! map from little group of q (G(q)) to space group (G)
            scqmap (nsymcrysq) = isym
        ! wrapping vector (reciprocal lattice vector)
            ivscwrapq (:, isym) = iv (:)
         End If
      End Do
      If (nsymcrysq .Lt. 1) Then
         Write (*, '(a, 3g18.10)') 'Error(' // thisnam // '): empty lit&
        &tle group of q for q:', vql
         Call terminate
      End If
End Subroutine findgroupq
!EOC
