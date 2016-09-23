!
!
! Copyright (C) 2014 S. Rigamonti and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genplot3d(fname)
  Use modmain
  Use modinput
  use modplotlabels
  ! arguments
  Character (*), Intent(In)  :: fname
  ! local variables
  Type (plot3d_type), pointer :: plot3d_
  type(plotlabels),pointer ::labels
  Real (8) :: zminl, zmaxl, zc
  Integer :: i

     
  If(associated(input%properties%stm%region)) Then
     
     zc = input%structure%crystal%basevect(3, 3)
     
     !zminl = input%properties%stm%region%zmin/zc
     !zmaxl = input%properties%stm%region%zmax/zc
     zminl = input%properties%stm%region%zrange(1)/zc
     zmaxl = input%properties%stm%region%zrange(2)/zc
     
     allocate(plot3d_)
     allocate(plot3d_%box)
     allocate(plot3d_%box%origin)
     allocate(plot3d_%box%pointarray(3))
     allocate(plot3d_%box%pointarray(1)%point)
     allocate(plot3d_%box%pointarray(2)%point)
     allocate(plot3d_%box%pointarray(3)%point)
     
     plot3d_%box%grid(1)=input%properties%stm%region%grid3d(1)
     plot3d_%box%grid(2)=input%properties%stm%region%grid3d(2)
     plot3d_%box%grid(3)=input%properties%stm%region%grid3d(3)
     plot3d_%box%origin%coord=(/0.0d0, 0.0d0, zminl/)
     plot3d_%box%pointarray(1)%point%coord=(/1.0d0, 0.0d0, zminl/)
     plot3d_%box%pointarray(2)%point%coord=(/0.0d0, 1.0d0, zminl/)
     plot3d_%box%pointarray(3)%point%coord=(/0.0d0, 0.0d0, zmaxl/)
     !plot3d_%box%origin%coord=(/0.0d0, 0.0d0, 0.0d0/)
     !plot3d_%box%pointarray(1)%point%coord=(/1.0d0, 0.0d0, 0.0d0/)
     !plot3d_%box%pointarray(2)%point%coord=(/0.0d0, 1.0d0, 0.0d0/)
     !plot3d_%box%pointarray(3)%point%coord=(/0.0d0, 0.0d0, 1.0d0/)

  Else
     write(*,*)
     write(*,'("Error(stm): stm element in topographic mode should")')
     write(*,'("contain one region element.")')
     write(*,*)
     stop
  End If

  !Call gencube (input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir, plot3d_)
  labels=>create_plotlablels("STM",fname,3)
  call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,4,"Wave Function Norm Squared","","graceunit")

  Call plot3d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir, plot3d_)

  call destroy_plotlablels(labels)

  Deallocate(plot3d_%box%pointarray(1)%point,plot3d_%box%pointarray(2)%point)
  Deallocate(plot3d_%box%pointarray)
  Deallocate(plot3d_%box%origin)
  Deallocate(plot3d_%box)
  Deallocate(plot3d_)

End Subroutine genplot3d

