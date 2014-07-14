!
!
! Copyright (C) 2014 S. Rigamonti and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genplot2d(fname)
  Use modmain
  Use modinput
  use modplotlabels
  ! arguments
  Character (*), Intent(In)  :: fname
  ! local variables
  Type (plot2d_type), pointer :: plot2d_
  type(plotlabels),pointer ::labels
  Real (8) :: heightl, zc
  Integer :: i

  If(associated(input%properties%stm%region)) Then
     
     ! only one of the z components must be different from zero
     Do i = 1, 3
        zc = input%structure%crystal%basevect(3, i)
        If (zc.gt.input%structure%epslat) exit
     End do

     heightl = input%properties%stm%region%height/zc

     allocate(plot2d_)
     allocate(plot2d_%parallelogram)
     allocate(plot2d_%parallelogram%origin)
     allocate(plot2d_%parallelogram%pointarray(2))
     allocate(plot2d_%parallelogram%pointarray(1)%point)
     allocate(plot2d_%parallelogram%pointarray(2)%point)

     plot2d_%parallelogram%grid(1)=input%properties%stm%region%grid2d(1)
     plot2d_%parallelogram%grid(2)=input%properties%stm%region%grid2d(2)
     plot2d_%parallelogram%pointarray(1)%point%coord=(/1.0d0, 0.0d0, heightl/)
     plot2d_%parallelogram%pointarray(2)%point%coord=(/0.0d0, 1.0d0, heightl/)
     plot2d_%parallelogram%origin%coord=(/0.0d0, 0.0d0, heightl/)
  Else If (associated(input%properties%stm%plot2d)) Then
     allocate(plot2d_)
     plot2d_=input%properties%stm%plot2d
  Else
     write(*,*)
     write(*,'("Error(stm): stm element in constantHeight mode should")')
     write(*,'("contain one region element or one plot2d element.")')
     write(*,*)
     stop
  End If

  labels=>create_plotlablels("2D STM image",fname,2)
  call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
  call set_plotlabel_axis(labels,3,"STM","","graceunit")
  !Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir, input%properties%stm%plot2d)
  write(*,*) "enter plot2d"
  Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir, plot2d_)
  write(*,*) "exit plot2d"
  call destroy_plotlablels(labels)

  If(associated(input%properties%stm%region)) Then
     Deallocate(plot2d_%parallelogram%pointarray(1)%point,plot2d_%parallelogram%pointarray(2)%point)
     Deallocate(plot2d_%parallelogram%pointarray)
     Deallocate(plot2d_%parallelogram%origin)
     Deallocate(plot2d_%parallelogram)
     Deallocate(plot2d_)
  End If

End Subroutine genplot2d

