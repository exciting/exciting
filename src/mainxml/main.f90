! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> Run exciting
!>
!> Unit test driver is run instead, if specified by passing the command
!> line argument to exciting.
!>
program main
   use inputdom
   use modinput
   use scl_xml_out_Module
   use modmpi
   use mod_device_offload, only: init_device_world, finish_device_world
   use mod_misc
   use cmd_line_args, only: cmd_line_args_type
   use unit_test_drivers, only: unit_test_driver
#ifdef USEOMP
   use omp_lib
#endif
   implicit none

   !> Command line arguments
   type(cmd_line_args_type) :: args

   call initmpi()
#ifdef USEOMP
   ! Old OpenMP standards allow for 
   ! nested parallelism. Since OpenMP 5.0 nested parallelism is
   ! no longer the default by the standard but the default is left to 
   ! the compiler. This call to the API ensures that 
   ! exciting does not use that, sometimes that leads 
   ! to frozen runs or garbage numbers in the results.
   call omp_set_max_active_levels(1)
#endif
   call init_device_world(mpiglobal%comm)
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
   call finish_device_world()
   call finitmpi()

end program
