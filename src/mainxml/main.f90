! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program main
use inputdom
use modinput
use scl_xml_out_Module
use modmpi
use mod_misc
implicit none
! set the version string
call versionfromdate()
call initMPI()
call loadinputDOM("input.xml")
!initialize default values that are not definet in the input schema
call setdefault()
!construct the input datastructure
input=getstructinput(inputnp)
call ifparseerrorstop()
call destroyDOM()
!some consistency checks of the input params
call checkinput()
!some initializations
call initatomcounters()
call initlattice()
call initlibxc()
call initldapu
call initsolver()
call readspeciesxml()
call scl_xml_out_create()
call tasklauncher()
call scl_xml_out_close()
call finitMPI()
end program
