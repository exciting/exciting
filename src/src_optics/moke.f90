!
!
!
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine moke
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: iw, iostat
      Complex (8) zt1, zt2, zt3
! allocatable arrays
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: sig1 (:, :)
      Real (8), Allocatable :: sig2 (:, :)
      Complex (8), Allocatable :: kerr (:)
! calculate dielectric function for the 11 and 12 components
      noptcomp = 2
      if (associated(input%properties%dielectric%optcomp)) &
        deallocate(input%properties%dielectric%optcomp)
      allocate(input%properties%dielectric%optcomp(3,noptcomp))
      input%properties%dielectric%optcomp(1, 1) = 1
      input%properties%dielectric%optcomp(2, 1) = 1
      input%properties%dielectric%optcomp(1, 2) = 1
      input%properties%dielectric%optcomp(2, 2) = 2
      Call dielectric
! allocate local arrays
      Allocate (w(input%properties%dos%nwdos))
      Allocate (sig1(input%properties%dos%nwdos, 2), &
     & sig2(input%properties%dos%nwdos, 2))
      Allocate (kerr(input%properties%dos%nwdos))
! read diagonal contribution to optical conductivity
      Open (50, File='SIGMA_11.OUT', Action='READ', Status='OLD', &
     & Form='FORMATTED', IoStat=IoStat)
      If (iostat .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(moke): error opening SIGMA_11.OUT")')
         Write (*,*)
         Stop
      End If
      Do iw = 1, input%properties%dos%nwdos
         Read (50, '(2G18.10)') w (iw), sig1 (iw, 1)
      End Do
      Read (50,*)
      Do iw = 1, input%properties%dos%nwdos
         Read (50, '(2G18.10)') w (iw), sig2 (iw, 1)
      End Do
      Close (50)
! read off-diagonal contribution to optical conductivity
      Open (50, File='SIGMA_12.OUT', Action='READ', Status='OLD', &
     & Form='FORMATTED', IoStat=IoStat)
      If (iostat .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(moke): error opening SIGMA_12.OUT")')
         Write (*,*)
         Stop
      End If
      Do iw = 1, input%properties%dos%nwdos
         Read (50, '(2G18.10)') w (iw), sig1 (iw, 2)
      End Do
      Read (50,*)
      Do iw = 1, input%properties%dos%nwdos
         Read (50, '(2G18.10)') w (iw), sig2 (iw, 2)
      End Do
      Close (50)
! calculate the complex Kerr angle
      Do iw = 1, input%properties%dos%nwdos
         If (w(iw) .Gt. 0.d0) Then
            zt1 = cmplx (sig1(iw, 1), sig2(iw, 1), 8)
            zt2 = cmplx (sig1(iw, 2), sig2(iw, 2), 8)
            zt3 = zt1 * Sqrt (1.d0+fourpi*zi*zt1/w(iw))
            If (Abs(zt3) .Gt. 1.d-8) Then
               kerr (iw) = - zt2 / zt3
            Else
               kerr (iw) = 0.d0
            End If
         Else
            kerr (iw) = 0.d0
         End If
      End Do
      Open (50, File='KERR.OUT', Action='WRITE', Form='FORMATTED')
      Do iw = 1, input%properties%dos%nwdos
         Write (50, '(2G18.10)') w (iw), dble (kerr(iw)) * 180.d0 / pi
      End Do
      Write (50, '("	  ")')
      Do iw = 1, input%properties%dos%nwdos
         Write (50, '(2G18.10)') w (iw), aimag (kerr(iw)) * 180.d0 / pi
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(moke):")')
      Write (*, '(" complex Kerr angle in degrees written to KERR.OUT")&
     &')
      Write (*,*)
      Deallocate (w, sig1, sig2, kerr)
      Return
End Subroutine
