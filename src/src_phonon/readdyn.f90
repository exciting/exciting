!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine readdyn (tsym,dynq)
      Use modmain
      Use modinput
      Implicit None
! arguments
      logical, intent(in) :: tsym
      Complex (8), Intent (Out) :: dynq (3*natmtot, 3*natmtot, nqpt)
! local variables
      Logical :: exist
      Integer :: iq, is, js, ia, ja
      Integer :: ip, jp, i, j
      Real (8) :: a, b
      Character (256) :: fext
! external functions
      Integer :: gcd
      External gcd
      Do iq = 1, nqpt
         i = 0
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               Do ip = 1, 3
                  i = i + 1
                  Call phfext (iq, is, ia, ip, fext)
                  Inquire (File='DYN'//trim(fext), Exist=Exist)
                  If ( .Not. exist) Then
                     Write (*,*)
                     Write (*, '("Error(readdyn): file not found :")')
                     Write (*, '(A)') ' DYN' // trim (fext)
                     Write (*,*)
                     Stop
                  End If
                  Open (50, File='DYN'//trim(fext), Action='READ', &
                 & Status='OLD', Form='FORMATTED')
                  j = 0
                  Do js = 1, nspecies
                     Do ja = 1, natoms (js)
                        Do jp = 1, 3
                           j = j + 1
                           Read (50,*) a, b
                           dynq (i, j, iq) = cmplx (a, b, 8)
                        End Do
                     End Do
                  End Do
                  Close (50)
               End Do
! end loops over atoms and species
            End Do
         End Do
! symmetrise the dynamical matrix
         if (tsym) Call dynsym (vql(:, iq), dynq(:, :, iq))
! end loop over q-vectors
      End Do
      Return
End Subroutine
