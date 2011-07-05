
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine rereadinput
         use inputdom
         use modinput
         implicit none
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
end subroutine

