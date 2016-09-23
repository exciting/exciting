!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: findsym
! !INTERFACE:
!
!
Subroutine findsym (apl1, apl2, nsym, lspl, lspn, iea)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   apl1 : first set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   apl2 : second set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   nsym : number of symmetries (out,integer)
!   lspl : spatial rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   lspn : spin rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   iea  : equivalent atom index for each symmetry
!          (out,integer(iea(natmmax,nspecies,48))
! !DESCRIPTION:
!   Finds the symmetries which rotate one set of atomic positions into another.
!   Both sets of positions differ only by a translation vector and have the same
!   muffin-tin magnetic fields (stored in the global array {\tt bfcmt}). Any
!   symmetry element consists of a spatial rotation of the atomic position
!   vectors followed by a global magnetic rotation: $\{\alpha_S|\alpha_R\}$. In
!   the case of spin-orbit coupling $\alpha_S=\alpha_R$. The symmetries are
!   returned as indices of elements in the Bravais lattice point group. An
!   index to equivalent atoms is stored in the array {\tt iea}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!   Fixed use of proper rotations for spin, February 2008 (L. Nordstrom)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: apl1 (3, maxatoms, maxspecies)
      Real (8), Intent (In) :: apl2 (3, maxatoms, maxspecies)
      Integer, Intent (Out) :: nsym
      Integer, Intent (Out) :: lspl (48)
      Integer, Intent (Out) :: lspn (48)
      Integer, Intent (Out) :: iea (natmmax, nspecies, 48)
! local variables
      Integer :: isym, jsym, jsym0, jsym1
      Integer :: is, ia, ja, iv (3), md
      Real (8) :: sl (3, 3), v (3), t1
! automatic arrays
      Integer :: jea (natmmax, nspecies)
      Real (8) :: apl3 (3, natmmax)
! external functions
      Real (8) :: r3taxi
      External r3taxi
      nsym = 0
! loop over lattice symmetries (spatial rotations)
      Do isym = 1, nsymlat
! make real copy of lattice rotation symmetry
         sl (:, :) = dble (symlat(:, :, isym))
! loop over species
         Do is = 1, nspecies
! map apl1 coordinates to [0,1) and store in apl3
            Do ia = 1, natoms (is)
               apl3 (:, ia) = apl1 (:, ia, is)
               Call r3frac (input%structure%epslat, apl3(:, ia), iv)
            End Do
            Do ja = 1, natoms (is)
! apply lattice symmetry to atomic positions
               Call r3mv (sl, apl2(:, ja, is), v)
! map coordinates to [0,1)
               Call r3frac (input%structure%epslat, v, iv)
! check if atomic positions are invariant
               Do ia = 1, natoms (is)
                  t1 = r3taxi (apl3(:, ia), v)
                  If (t1 .Lt. input%structure%epslat) Then
! equivalent atom index
                     jea (ia, is) = ja
                     Go To 10
                  End If
               End Do
! not invariant so try new spatial rotation
               Go To 40
10             Continue
            End Do
         End Do
! all atomic positions invariant at this point
         jsym = 1
! spin polarised case
         If (associated(input%groundstate%spin)) Then
! check invariance of magnetic fields under global spin rotation
            If (isspinorb()) Then
! with spin-orbit coupling spin rotation equals spatial rotation
               jsym0 = isym
               jsym1 = isym
            Else
! without spin-orbit coupling spin rotation independent of spatial rotation
               jsym0 = 1
               jsym1 = nsymlat
            End If
            Do jsym = jsym0, jsym1
! determinant of the symmetry matrix
               md = symlatd (jsym)
! rotate global field and check invariance using proper part of symmetry matrix
               Call r3mv (symlatc(:, :, jsym), &
              & input%groundstate%spin%bfieldc, v)
               v (:) = v (:) * dble (md)
               t1 = r3taxi (input%groundstate%spin%bfieldc, v)
! if not invariant try a different global spin rotation
               If (t1 .Gt. input%structure%epslat) Go To 20
! rotate muffin-tin magnetic fields and check invariance
               Do is = 1, nspecies
                  Do ia = 1, natoms (is)
! equivalent atom
                     ja = jea (ia, is)
                     Call r3mv (symlatc(:, :, jsym), input%structure%speciesarray(is)%species%atomarray(ja)%atom%bfcmt(:), v)
                     v (:) = v (:) * dble (md)
                     t1 = r3taxi (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:), v)
! if not invariant try a different global spin rotation
                     If (t1 .Gt. input%structure%epslat) Go To 20
                  End Do
               End Do
! all fields invariant
               Go To 30
20             Continue
! end loop over global spin rotations
            End Do
! magnetic fields not invariant so try different spatial rotation
            Go To 40
         End If
30       Continue
! everything invariant so add symmetry to set
         nsym = nsym + 1
         lspl (nsym) = isym
         lspn (nsym) = jsym
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               iea (ia, is, nsym) = jea (ia, is)
            End Do
         End Do
40       Continue
! end loop over spatial rotations
      End Do
      Return
End Subroutine
!EOC
