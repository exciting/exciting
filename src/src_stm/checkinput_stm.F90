!
!
! Copyright (C) 2014 S. Rigamonti and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine checkinput_stm()
  Use modmain
  Use modinput
  use modplotlabels

  ! local variables
  Integer :: i, count, iz
  Real(8) :: xc,yc,zc

  If(associated(input%properties%stm%region)) Then
     count = 0
     Do i = 1, 3
        zc = input%structure%crystal%basevect(3, i)
        If (Abs(zc).gt.input%structure%epslat) Then
           count=count+1
           iz = i
        End If
     End do
     If (count .Gt. 1) Then
        write(*,*)
        write(*,'("Error(stm): More than one lattice vector with nonzero z component detected.")')
        write(*,'("When using stm/plane element, it is assumed that")')
        write(*,'("only one lattice vector has nonzero z component. Stopping")')
        write(*,*)
        stop
     End If
     xc =  input%structure%crystal%basevect(1, iz)
     yc =  input%structure%crystal%basevect(2, iz)
     If ((Abs(xc).gt.input%structure%epslat).or.(Abs(yc).gt.input%structure%epslat)) Then
        write(*,*)
        write(*,'("Error(stm): Lattice vector defining the normal to the surface has non-null x or y component.")')
        write(*,'("When using stm/plane element, it is assumed that")')
        write(*,'("the lattice vector defining the normal to the surface is orthogonal to the xy plane. Stopping")')
        write(*,*)
        stop
     End If
  End If
End Subroutine checkinput_stm



