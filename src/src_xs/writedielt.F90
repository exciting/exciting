!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Subroutine writedielt (filtag, nw, w, dt, switch)
      Use modmain
      Use m_getunit
      Implicit None
  ! arguments
      Character (*), Intent (In) :: filtag
      Integer, Intent (In) :: nw
      Real (8), Intent (In) :: w (nw)
      Complex (8), Intent (In) :: dt (3, 3, nw)
      Integer, Intent (In) :: switch
  ! local variables
      Integer :: un, oct, iw
      Call getunit (un)
      Open (un, File=trim(filtag)//trim(filext), Form='formatted', &
     & Action='write', Status='replace')
      Write (un,*)
      If (switch .Eq. 0) Then
         Write (un, '(" (dielectric tensor, independent particle approx&
        &imation)")')
      Else
         Write (un, '(" (dielectric tensor, random phase approximation)&
        &")')
      End If
      Write (un,*)
      Do iw = 1, nw
         Write (un, '(" frequency index and value: ",i6,f14.8)') iw, w &
        & (iw)
         Write (un, '(" real part, imaginary part below")')
         Write (un, '(3f14.8,5x,3f14.8)') (dble(dt(oct, :, iw)), &
        & aimag(dt(oct, :, iw)), oct=1, 3)
         Write (un,*)
      End Do
      Close (un)
End Subroutine writedielt
