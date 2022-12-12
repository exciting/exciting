! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> Run exciting
!>
!> Unit test driver is run instead, if specified by passing the command
!> line argument to exciting.
!>
!> When the build system gets updated to CMake, unit tests will be compiled
!> as executables that link to the appropriate exciting modules, removing
!> the need to call the unit test driver in this manner.
program main
   use inputdom
   use modinput
   use scl_xml_out_Module
   use modmpi
   use mod_misc
   use cmd_line_args, only: cmd_line_args_type
   use unit_test_drivers, only: unit_test_driver
   implicit none

   !> Command line arguments
   type(cmd_line_args_type) :: args

   call initmpi()
   call versionfromdate()
   call args%parse(mpiglobal)
   if (args%run_unit_tests) then
      call unit_test_driver(mpiglobal, args%kill_on_failure)
   else
      call loadinputDOM("input.xml")
      ! Initialise default values that are not defined in the input schema
      call setdefault()
      ! Construct the input datastructure
      input = getstructinput(inputnp)
      call ifparseerrorstop()
      call destroyDOM()
      ! Some consistency checks of the input params
      call checkinput()
      ! Some initializations
      call initatomcounters()
      call initlattice()
      call initlibxc()
      call initldapu
      call initsolver()
      call readspeciesxml()
      call scl_xml_out_create()
      !See if there'anything useful in testingfun
      !call testingfun
      call tasklauncher()
      call scl_xml_out_close()
   end if
   call finitmpi()

end program
