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
      Integer :: is, ia, ias
      Real (8) :: v (3)
      If (topt) Then
! GB 04.10.2012 new file with geometrical optimization history (xyz format)
     If (input%structureoptimization%history) Then
         Open (51, File='history.xyz', Action='WRITE', &
        & POSITION='APPEND')
! GB 04.10.2012   END
      End If
      End If

      Write (51,*) natmtot
     If (input%structureoptimization%history) Then
      Write (51,*) "Total Energy =",  engytot*27.211396641344194, " eV "
     End If
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
              & v)
           Write (51,'(A6, 3F14.8, 3F14.8)')(input%structure%speciesarray(is)%species%chemicalSymbol), ((v (:)))*0.529177249, &
                & (forcetot (:, ias))*51.42208245
! GB 4.10.2012 END
         End Do
      End Do
      Close (51)
      Return
End Subroutine
!EOC
