!
!
!
!
Subroutine phfext (iq, is, ia, ip, iph, istepph, fext, fextdyn, dirname)
      Use modmain
!     Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
      Integer, Intent (In) :: iph
      Integer, Intent (In) :: istepph
      Character (256), Intent (Out) :: fext
      Character (256), Intent (Out) :: fextdyn
      Character (256), Intent (Out) :: dirname
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
      ! for gndstate output use full suffix incl cos, sin, step
      if (iph .eq. 0) then
         ! cos displacement
         Write (fext, '("_Q", 2I2.2, "_", 2I2.2, "_", 2I2.2, "_S", I2.2, "&
        &_A", I3.3, "_P", I1, "_C", I1, ".OUT")') m (1), n (1), m (2), n (2), m (3), &
        & n (3), is, ia, ip, istepph
      else
         ! sin displacement
         Write (fext, '("_Q", 2I2.2, "_", 2I2.2, "_", 2I2.2, "_S", I2.2, "&
        &_A", I3.3, "_P", I1, "_S", I1, ".OUT")') m (1), n (1), m (2), n (2), m (3), &
        & n (3), is, ia, ip, istepph
      endif
      ! for DYN files use only q, species, atom, polarization
      Write (fextdyn, '("_Q", 2I2.2, "_", 2I2.2, "_", 2I2.2, "_S", I2.2, "&
     &_A", I3.3, "_P", I1, ".OUT")') m (1), n (1), m (2), n (2), m (3), &
     & n (3), is, ia, ip
      ! subdirectory name = filextdyn without .OUT
      Write (dirname, '("Q", 2I2.2, "_", 2I2.2, "_", 2I2.2, "_S", I2.2, "&
     &_A", I3.3, "_P", I1, "/")') m (1), n (1), m (2), n (2), m (3), &
     & n (3), is, ia, ip
      Return
End Subroutine
