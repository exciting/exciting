!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine backup0
      Use modmain
      Use modinput
      Use modxs
      Implicit None
      filext_b = trim (filext)
      nosym_b = input%groundstate%nosym
      swidth_b = input%groundstate%swidth
      lmaxapw_b = input%groundstate%lmaxapw
      lmaxmat_b = input%groundstate%lmaxmat
      maxscl_b = input%groundstate%maxscl
      If (associated(input%groundstate%spin)) Then
         bfieldc_b = input%groundstate%spin%bfieldc
      End If
End Subroutine backup0

