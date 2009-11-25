!
!
!
!
Subroutine phfext (iq, is, ia, ip, fext)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
      Character (256), Intent (Out) :: fext
! local variables
      Integer :: i, j, m (3), n (3)
! external functions
      Integer :: gcd
      External gcd
      Do i = 1, 3
         If (ivq(i, iq) .Ne. 0) Then
            j = gcd (ivq(i, iq), ngridq(i))
            m (i) = ivq (i, iq) / j
            n (i) = ngridq (i) / j
         Else
            m (i) = 0
            n (i) = 0
         End If
      End Do
      Write (fext, '("_Q", 2I2.2, "_", 2I2.2, "_", 2I2.2, "_S", I2.2, "&
     &_A", I3.3, "_P", I1, ".OUT")') m (1), n (1), m (2), n (2), m (3), &
     & n (3), is, ia, ip
      Return
End Subroutine
