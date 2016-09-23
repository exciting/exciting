!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmdexcdc (ikp, evecsv, dedc)
! calculate the derivative of exchange-correlation energy w.r.t. evecsv
      Use modinput
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: ikp
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Complex (8), Intent (Inout) :: dedc (nstsv, nstsv)
! local variables
      Integer :: ik, jk, iv (3)
      Integer :: ist1, ist2, ist3, ist4
      Real (8) :: t1, t2
! allocatable arrays
      Complex (8), Allocatable :: vnl (:, :, :, :)
! external functions
      Real (8) :: r3taxi
      External r3taxi
      If (input%groundstate%RDMFT%rdmxctype .Eq. 0) Return
! calculate the prefactor
      If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
         t1 = 1.d0 / occmax
      Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) Then
         If (associated(input%groundstate%spin)) Then
            t1 = 1.d0
         Else
            t1 = 2.d0 * (0.25d0) ** input%groundstate%RDMFT%rdmalpha
         End If
      Else
         Write (*,*)
         Write (*, '("Error(rdmdexcdc): rdmxctype not defined : ", I8)') input%groundstate%RDMFT%rdmxctype
         Write (*,*)
         Stop
      End If
      Allocate (vnl(nstsv, nstsv, nstsv, nkptnr))
! calculate non-local matrix elements of the type (l-jj-k)
      Call rdmvnlc (ikp, vnl)
! start loop over non-reduced k-points
      Do ik = 1, nkptnr
! copy the matrix elements of the type i-jj-i to vnlrdm
         Do ist1 = 1, nstsv
            Do ist2 = 1, nstsv
               vnlrdm (ist1, ikp, ist2, ik) = dble (vnl(ist1, ist1, &
              & ist2, ik))
            End Do
         End Do
! find the equivalent reduced k-point
         iv (:) = ivknr (:, ik)
         jk = ikmap (iv(1), iv(2), iv(3))
         Do ist1 = 1, nstsv
            Do ist2 = 1, nstsv
               Do ist3 = 1, nstsv
                  Do ist4 = 1, nstsv
                     If (input%groundstate%RDMFT%rdmxctype .Eq. 1) Then
! Hartree-Fock functional
                        t2 = t1 * occsv (ist3, ikp) * occsv (ist4, jk)
                     Else If (input%groundstate%RDMFT%rdmxctype .Eq. 2) &
                    & Then
! SDLG functional
                        If ((ist3 .Eq. ist4) .And. (r3taxi(vkl(1, ikp), &
                       & vklnr(1, jk)) .Lt. input%structure%epslat)) &
                       & Then
                           t2 = (1.d0/occmax) * occsv (ist4, jk) ** 2
                        Else
                           t2 = t1 * (occsv(ist3, ikp)*occsv(ist4, jk)) &
                          & ** input%groundstate%RDMFT%rdmalpha
                        End If
                     End If
                     dedc (ist2, ist3) = dedc (ist2, ist3) - t2 * &
                    & evecsv (ist2, ist1) * vnl (ist1, ist3, ist4, ik)
                  End Do
               End Do
            End Do
         End Do
! end loop over non-reduced k-points
      End Do
      Deallocate (vnl)
      Return
End Subroutine
