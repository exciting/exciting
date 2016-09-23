!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genlocmt (ngp, isti, istf, evecfv, wfcmt)
      Use modmain
      Implicit None
  ! arguments
      Integer, Intent (In) :: ngp, isti, istf
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (Out) :: wfcmt (istf-isti+1, &
     & nlomax,-lolmax:lolmax, natmtot)
  ! local variables
      Integer :: ist, istc, is, ia, ias, i, ilo, l, m, lm
      Do istc = isti, istf
         ist = istc - isti + 1
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ilo = 1, nlorb (is)
                  l = lorbl (ilo, is)
                  Do m = - l, l
                     lm = idxlm (l, m)
                     i = idxlo (lm, ilo, ias)
                     wfcmt (ist, ilo, m, ias) = evecfv (ngp+i, istc, 1)
                  End Do
               End Do
            End Do
         End Do
      End Do
End Subroutine genlocmt
