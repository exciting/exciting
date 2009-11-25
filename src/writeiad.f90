!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writeiad
! !INTERFACE:
!
!
Subroutine writeiad (topt)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   topt : if .true. then the filename will be {\tt IADIST_OPT.OUT}, otherwise
!          {\tt IADIST.OUT} (in,logical)
! !DESCRIPTION:
!   Outputs the interatomic distances to file.
!
! !REVISION HISTORY:
!   Created May 2005 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: topt
! local variables
      Integer :: is, js, ia, ja
      Integer :: i1, i2, i3
      Real (8) :: d, dmin, v (3)
! external functions
      Real (8) :: r3dist
      External r3dist
      If (topt) Then
         Open (50, File='IADIST_OPT'//trim(filext), Action='WRITE', &
        & Form='FORMATTED')
      Else
         Open (50, File='IADIST'//trim(filext), Action='WRITE', Form='F&
        &ORMATTED')
      End If
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Write (50,*)
            Write (50, '("Distance between is = ", I4, " (", A, "), ia &
           &= ", I4, " and")') is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol), &
           & ia
            Do js = 1, nspecies
               Do ja = 1, natoms (js)
                  dmin = 1.d8
                  Do i1 = - 1, 1
                     Do i2 = - 1, 1
                        Do i3 = - 1, 1
                           v (:) = dble (i1) * &
                          & input%structure%crystal%basevect(:, 1) + &
                          & dble (i2) * &
                          & input%structure%crystal%basevect(:, 2) + &
                          & dble (i3) * &
                          & input%structure%crystal%basevect(:, 3) + &
                          & atposc (:, ja, js)
                           d = r3dist (atposc(:, ia, is), v)
                           dmin = Min (d, dmin)
                        End Do
                     End Do
                  End Do
                  Write (50, '(" is = ", I4, " (", A, "), ia = ", I4, " : ", G18.10)') js, trim &
                 & (input%structure%speciesarray(js)%species%chemicalSymbol), ja, dmin
               End Do
            End Do
         End Do
      End Do
      Close (50)
      Return
End Subroutine
!EOC
