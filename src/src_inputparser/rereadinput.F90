
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine rereadinput
         use inputdom
         use modinput
         implicit none
         Integer :: is, ia         
         ! read in input again to reset atomic positions and lattice vectors in particular
         call loadinputDOM("input.xml")
         call setdefault
         input=getstructinput(inputnp)
         call ifparseerrorstop()
         call destroyDOM()
         call checkinput()
         call initatomcounters()
         call initlattice()
         call initlibxc()
         call initldapu
         call initsolver()
         
         ! Assuming that this routine (rereadinput) is not called before a scf cycle,
         ! the b-fields should be zero, if the user intended to reduce them every scf iteration. 
         If (associated(input%groundstate%spin)) Then
            If (input%groundstate%spin%reducebf .Lt. 1.d0) Then
               input%groundstate%spin%bfieldc(:) = 0.d0
               Do is = 1, size(input%structure%speciesarray(:))
                  Do ia = 1, size(input%structure%speciesarray(is)%species%atomarray(:))
                     input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) = 0.d0
                  End Do
               End Do
            End If
         End If
         
end subroutine

