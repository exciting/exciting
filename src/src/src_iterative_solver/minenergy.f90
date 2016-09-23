
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine minenergy (sigma)
      Use modmain
      Implicit None
      Complex (8), Intent (Out) :: sigma
!
      Real (8) :: lorbearray (maxlorbord*maxlorb*maxspecies)
      Integer :: en, i, j, k, m
      Integer, External :: idamax
      Real (8), Parameter :: one = 1
      m = 1
      lorbearray = 0
      en = maxlorbord * maxlorb * maxspecies
      Do i = 1, maxlorbord
         Do j = 1, maxlorb
            Do k = 1, maxspecies
               If (Abs(lorbe0(i, j, k)) .Lt. 20) lorbearray (m) = &
              & lorbe0 (i, j, k)
               m = m + 1
            End Do
         End Do
      End Do
      Call dscal (en,-one, lorbearray, 1)
      sigma = dcmplx (Min(-Abs(lorbearray(idamax(en, lorbearray, &
     & 1))),-1.d0))
!write(*,*)sigma
End Subroutine
