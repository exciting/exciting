!
!
!
!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: findsymcrys
! !INTERFACE:
!
!
Subroutine findsymcrys
! !USES:
      Use modinput
      Use modmain
#ifdef XS
      Use modxs
#endif
! !DESCRIPTION:
!   Finds the complete set of symmetries which leave the crystal structure
!   (including the magnetic fields) invariant. A crystal symmetry is of the
!   form $\{\alpha_S|\alpha_R|{\bf t}\}$, where ${\bf t}$ is a translation
!   vector, $\alpha_R$ is a spatial rotation operation and $\alpha_S$ is a
!   global spin rotation. Note that the order of operations is important and
!   defined to be from right to left, i.e. translation followed by spatial
!   rotation followed by spin rotation. In the case of spin-orbit coupling
!   $\alpha_S=\alpha_R$. In order to determine the translation vectors, the
!   entire atomic basis is shifted so that the first atom in the smallest set of
!   atoms of the same species is at the origin. Then all displacement vectors
!   between atoms in this set are checked as possible symmetry translations. If
!   the global variable {\tt tshift} is set to {\tt .false.} then the shift is
!   not performed. See L. M. Sandratskii and P. G. Guletskii, {\it J. Phys. F:
!   Met. Phys.} {\bf 16}, L43 (1986) and the routine {\tt findsym}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ia, ja, is, js, i, n
      Integer :: isym, nsym, iv (3)
      Integer :: lspl (48), lspn (48)
      Real (8) :: v (3), t1
      Real (8) :: apl (3, maxatoms, maxspecies), aplt (3, maxatoms, &
     & maxspecies)
!
! allocatable arrays
      Integer, Allocatable :: iea (:, :, :)
      Real (8), Allocatable :: vtl (:, :)
! allocate local array
      Allocate (iea(natmmax, nspecies, 48))
! allocate equivalent atom arrays
      If (allocated(ieqatom)) deallocate (ieqatom)
      Allocate (ieqatom(natmmax, nspecies, maxsymcrys))
      If (allocated(eqatoms)) deallocate (eqatoms)
      Allocate (eqatoms(natmmax, natmmax, nspecies))
! find the smallest set of atoms
      is = 1
      Do js = 1, nspecies
         If (natoms(js) .Lt. natoms(is)) is = js
      End Do
      If ((input%structure%tshift) .And. (natmtot .Gt. 0)) Then
! shift basis so that the first atom in the smallest atom set is at the origin
         v (:) = input%structure%speciesarray(is)%species%atomarray(1)%atom%coord(:)
         Do js = 1, nspecies
            Do ia = 1, natoms (js)
! shift atom
               input%structure%speciesarray(js)%species%atomarray(ia)%atom%coord(:) = &
              & input%structure%speciesarray(js)%species%atomarray(ia)%atom%coord(:) - v (:)
! map lattice coordinates back to [0,1)
               Call r3frac (input%structure%epslat, input%structure%speciesarray(js)%species%atomarray(ia)%atom%coord(:), iv)
! determine the new Cartesian coordinates
               Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(js)%species%atomarray(ia)%atom%coord(:), &
              & atposc(:, ia, js))
            End Do
         End Do
      End If
! determine possible translation vectors from smallest set of atoms
      n = Max (natoms(is)*natoms(is), 1)
      Allocate (vtl(3, n))
      n = 1
      vtl (:, 1) = 0.d0
      Do ia = 1, natoms (is)
         Do ja = 2, natoms (is)
            v (:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:) - &
           & input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord(:)
            Call r3frac (input%structure%epslat, v, iv)
            Do i = 1, n
               t1 = Abs (vtl(1, i)-v(1)) + Abs (vtl(2, i)-v(2)) + Abs &
              & (vtl(3, i)-v(3))
               If (t1 .Lt. input%structure%epslat) Go To 10
            End Do
            n = n + 1
            vtl (:, n) = v (:)
10          Continue
         End Do
      End Do
      eqatoms (:, :, :) = .False.
      nsymcrys = 0
! loop over all possible translations
      Do i = 1, n
! construct new array with translated positions
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               apl (:, ia, is) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:) + vtl (:, i)
               aplt (:, ia, is) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            End Do
         End Do
! find the symmetries for current translation
         Call findsym (aplt, apl, nsym, lspl, lspn, iea)
         Do isym = 1, nsym
#ifdef XS
     ! exclude non-zero translations and check if spin-rotations is equal to spatial one
            If (associated(input%groundstate)) Then
               if (input%groundstate%symmorph) then
                  if (sum(Abs(vtl(:, i))) .Gt. input%structure%epslat) Go To 20
                  if (associated(input%groundstate%spin) .and. &
                     (lspnsymc(isym) .ne. lsplsymc(isym))) goto 20
               end if
            End If
#endif
            nsymcrys = nsymcrys + 1
            If (nsymcrys .Gt. maxsymcrys) Then
               Write (*,*)
               Write (*, '("Error(findsymcrys): too many symmetries")')
               Write (*, '(" Adjust maxsymcrys in modmain and recompile&
              & code")')
               Write (*,*)
               Stop
            End If
            vtlsymc (:, nsymcrys) = vtl (:, i)
            lsplsymc (nsymcrys) = lspl (isym)
            lspnsymc (nsymcrys) = lspn (isym)
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ja = iea (ia, is, isym)
                  ieqatom (ia, is, nsymcrys) = ja
                  eqatoms (ia, ja, is) = .True.
                  eqatoms (ja, ia, is) = .True.
               End Do
            End Do
#ifdef XS
20          Continue
#endif
         End Do
      End Do
      Deallocate (iea, vtl)
      Return
End Subroutine
!EOC
