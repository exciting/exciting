!
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Subroutine findsymgenr (nc, mt, ngen, nsgen, gen, orbgen)
      Use modsym
      Implicit None
  ! arguments
      Integer, Intent (In) :: nc
      Integer, Intent (In) :: mt (nc, nc)
      Integer, Intent (Out) :: ngen
      Integer, Intent (Out) :: nsgen (nc)
      Integer, Intent (Out) :: gen (nc)
      Integer, Intent (Out) :: orbgen (nc, nc)
  ! local variables
      Integer :: i, is, j, js, norb (nc), orb (nc, nc), idx (nc), done &
     & (nc)
  ! set up orbits of group elements
      norb (:) = 0
      norb (1) = 1
      Do i = 1, nc
         orb (i, 1) = i
         Do j = 2, nc
            orb (i, j) = mt (orb(i, j-1), i)
            If (orb(i, j) .Eq. 1) Then
               norb (i) = j - 1
               Exit
            End If
         End Do
      End Do
  ! sort orbits according to their number of elements
      Call sortidx (nc, dble(norb), idx)
  ! add largest generator to set
      ngen = 1
      gen (:) = 0
      gen (ngen) = orb (idx(nc), 1)
      nsgen (:) = 0
      nsgen (1) = norb (idx(nc))
      orbgen (:, :) = 0
      orbgen (1, :nsgen(1)) = orb (idx(nc), :nsgen(1))
      Do i = nc - 1, 1, - 1
         is = idx (i)
         Do j = nc, i + 1, - 1
            js = idx (j)
        ! discard orbit if one of its elements is equal to previous generators
            If (any(orb(js, :norb(js)) .Eq. orb(is, 1))) Go To 10
        ! discard trivial generator
            If (orb(is, 1) .Eq. 1) Go To 10
         End Do
     ! add new generator to set
         ngen = ngen + 1
         nsgen (ngen) = norb (is)
         gen (ngen) = orb (is, 1)
         orbgen (ngen, :nsgen(ngen)) = orb (is, :nsgen(ngen))
10       Continue
      End Do
  ! check if orbits cover the symmetry group
      done (:) = 0
      done (1) = ngen
      Do i = 1, ngen
         done (orbgen(i, :nsgen(i))) = done (orbgen(i, :nsgen(i))) + 1
      End Do
      If (any(done .Eq. 0)) Then
         Write (*,*)
         Write (*, '("Error(findsymgenr): Generators do not cover the s&
        &ymmetry group")')
         Write (*,*)
         Stop
      End If
  ! add identity to orbits
      Do i = 1, ngen
         If (gen(i) .Ne. 1) Then
            nsgen (i) = nsgen (i) + 1
            orbgen (i, nsgen(i)) = 1
         End If
      End Do
End Subroutine findsymgenr
