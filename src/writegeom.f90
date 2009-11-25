!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writegeom
! !INTERFACE:
!
!
Subroutine writegeom (topt)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   topt : if .true. then the filename will be {\tt GEOMETRY_OPT.OUT}, otherwise
!          {\tt GEOMETRY.OUT} (in,logical)
! !DESCRIPTION:
!   Outputs the lattice vectors and atomic positions to file, in a format
!   which may be then used directly in {\tt exciting.in}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: topt
! local variables
      Integer :: is, ia
      Real (8) :: v (3)
      If (topt) Then
         Open (50, File='GEOMETRY_OPT'//trim(filext), Action='WRITE', &
        & Form='FORMATTED')
      Else
         Open (50, File='GEOMETRY'//trim(filext), Action='WRITE', &
        & Form='FORMATTED')
      End If
      Write (50,*)
      Write (50, '("scale")')
      Write (50, '(" 1.0")')
      Write (50,*)
      Write (50, '("scale1")')
      Write (50, '(" 1.0")')
      Write (50,*)
      Write (50, '("scale2")')
      Write (50, '(" 1.0")')
      Write (50,*)
      Write (50, '("scale3")')
      Write (50, '(" 1.0")')
      Write (50,*)
      Write (50, '("avec")')
      Write (50, '(3G18.10)') input%structure%crystal%basevect(:, 1)
      Write (50, '(3G18.10)') input%structure%crystal%basevect(:, 2)
      Write (50, '(3G18.10)') input%structure%crystal%basevect(:, 3)
      If (input%structure%molecule) Then
         Write (50,*)
         Write (50, '("molecule")')
         Write (50, '(" ", L1)') input%structure%molecule
      End If
      Write (50,*)
      Write (50, '("atoms")')
      Write (50, '(I4, T40, " : nspecies")') nspecies
      Do is = 1, nspecies
         Write (50, '(" ''", A, "''", T40, " : spfname")') trim &
        & (input%structure%speciesarray(is)%species%speciesfile)
         Write (50, '(I4, T40, " : natoms; atpos, bfcmt below")') &
        & natoms (is)
         Do ia = 1, natoms (is)
            If (input%structure%molecule) Then
! write Cartesian coordinates for the molecular case
               Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
              & v)
            Else
! otherwise write lattice coordinates
               v (:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            End If
            Write (50, '(3F14.8, "  ", 3F12.8)') v (:), input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
         End Do
      End Do
      Close (50)
      Return
End Subroutine
!EOC
