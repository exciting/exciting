!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine restore0
      Use modmain
      Use modinput
      Use modxs
      Implicit None
      filext = trim (filext_b)
      input%groundstate%nosym = nosym_b
      input%groundstate%swidth = swidth_b
      input%groundstate%lmaxapw = lmaxapw_b
      input%groundstate%lmaxmat = lmaxmat_b
      input%groundstate%maxscl = maxscl_b
      If (associated(input%groundstate%spin)) Then
         input%groundstate%spin%bfieldc = bfieldc_b
      End If
End Subroutine restore0
