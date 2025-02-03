! Copyright (C) 2024 Deviceacc developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

!> @file
!> This file contains a test of a pure openmp offload to the GPU
!> notice that the library in this case is only used to
!> properly associate each host process to a device.
program test_do_loop_gpu

   use iso_fortran_env,         only: i32=>int32, r32=> real32, r64=>real64
   use iso_c_binding
   use m_device_world_t,        only: device_world_t
#if defined(DEVICEOFFLOAD)
   use mpi
#endif
   implicit none

   integer(i32), parameter :: n = 6000
   integer(i32)            :: i
   real(r32)               :: a(n), b(n), c(n)
   type(device_world_t)    :: device_world
   integer(i32)            :: mpi_world, err
   real(r32)               :: rdiff
   real(r32), parameter    :: tolerance = 1.0e-06_r64

#if defined(DEVICEOFFLOAD)
   ! Init MPI world
   call mpi_init(err)
   mpi_world = MPI_COMM_WORLD
#endif

   ! Init GPU stuff
    call device_world%init(mpi_world)

   ! Set data in the device
   ! Contrary to the register stuff map fails for declared type elements
   ! as well as sometimes between objects in different files
   ! Nevertheless, here we see a pure openmp way

   do i = 1, n
      a(i) = i
      b(i) = i * 2
   end do

   !$omp target teams distribute parallel do map(to: a, b) map(from: c)
   do i = 1, N
      c(i) = sum_gpu(a(i), b(i))
   end do
   !$omp end target teams distribute parallel do

   ! Check result
   write(*,*) '[TEST : test_do_loop_gpu]'
   write(*,'(A, E16.7)')  '  * Regression result (relative difference): ', rdiff
   rdiff = norm2( c - a - b) / norm2(a+b)
   if ( rdiff .lt. tolerance) then
      write(*,*)  '[TEST : test_do_loop_gpu : PASSED]'
   else
      write(*,*)  '[TEST : test_do_loop_gpu : FAILED]'
      stop 1
   end if

#if defined(DEVICEOFFLOAD)
   call mpi_finalize(err)
#endif

contains

   ! Example of how external function in GPU should look
   ! Notice that for gfortran most of the intrinsics are implemented
   ! For ifx that is not the case, few are implemented
   ! Also keep in mind that you should not allocate memory within a device
   ! block
   function sum_gpu(a1, a2) result(a3)
      real(r32), intent(in) :: a1
      real(r32), intent(in) :: a2

      !! When inside a procedure block compiler will try to also compile it for the GPU
      !$omp declare target

      real(r32) :: a3

      a3 = a1 + a2

   end function sum_gpu

end program test_do_loop_gpu

