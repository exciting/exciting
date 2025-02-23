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
!> This file contains a test of zgemm using GPU accelerated routines
program test_zgemm
    
    use iso_fortran_env,         only: i32=>int32, r64=>real64 
    use iso_c_binding
    use m_device_world_t,        only: device_world_t
    use device_linalg_common_interface, only: zgemm_gpu
#if defined(DEVICEOFFLOAD)
    use mpi
#endif
    implicit none

    integer(i32), parameter   :: n = 4000_i32
    complex(r64), parameter   :: zzero = cmplx(0.0, 0.0, kind=r64)
    complex(r64), parameter   :: zone  = cmplx(1.0, 0.0, kind=r64)
    complex(r64), target, allocatable :: A(:,:), B(:,:), C(:,:), C_ref(:,:)
    type(device_world_t)      :: device_world
    integer(i32)              :: i, j
    real(r64)                 :: r(4)
    real(r64)                 :: rdiff
    real(r64), parameter      :: tolerance = 1.0e-12_r64
    integer(i32)              :: mpi_world, err
    
    ! LAPACK
    type(c_ptr) :: dA_ptr, dB_ptr, dC_ptr
#if defined(DEVICEOFFLOAD)
    ! Init MPI world
    call mpi_init(err)
    mpi_world = MPI_COMM_WORLD
#endif
    ! Allocating the matrix A
    allocate(A(n,n), B(n,n), C(n,n),C_ref(n,n))

    do i = 1, n
        do j = 1, n
            call random_number(r)
            A(j,i) = cmplx(r(1),r(2), r64)
            B(j,i) = cmplx(r(3),r(4), r64)
        end do
    end do 

    ! Compute reference
    call zgemm('t', 'c', n, n, n, zone, A, n, B, n, zzero, C_ref, n)

    ! Init GPU stuff
    call device_world%init(mpi_world)

    !Allocate in device A, associate to host A, and copy from host to device
    call device_world%register%alloc("A", size(A) * c_sizeof(zzero), device_world%get_device())
    call device_world%register%assoc("A", C_loc(A))
    call device_world%register%to_device("A")
    ! Do the same with B
    call device_world%register%alloc("B", size(B) * c_sizeof(zzero), device_world%get_device())
    call device_world%register%assoc("B", C_loc(B))
    call device_world%register%to_device("B")
    ! Idem C but not need to copy anything
    call device_world%register%alloc("C", size(C) * c_sizeof(zzero), device_world%get_device())
    call device_world%register%assoc("C", C_loc(C))

    ! Obtain device pointers
    dA_ptr = device_world%register%device_ptr("A")
    dB_ptr = device_world%register%device_ptr("B")
    dC_ptr = device_world%register%device_ptr("C")
    
    ! Performing the actual zgemm in the GPU
    call zgemm_gpu('t', 'c', n, n, n, zone, dA_ptr, n, dB_ptr, n, zzero, dC_ptr, n, device_world)

    ! Not needed but example of sync
    call device_world%synchronize()

    ! Retrieve C from the device, and free there the A, B, C space
    call device_world%register%from_device("C")
    call device_world%register%remove("A")
    call device_world%register%remove("B")
    call device_world%register%remove("C")

    call device_world%finish()
#if defined(DEVICEOFFLOAD)
    call mpi_finalize(err)
#endif
    ! Check
    rdiff = sum(abs(C - C_ref)) / sum(abs(C))  
    write(*,*) '[TEST : test_zgemm]'   
    write(*,'(A, E16.7)')  '  * Regression result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : test_zgemm : PASSED]'
    else
        write(*,*)  '[TEST : test_zgemm : FAILED]'
        stop 1
    end if

    stop 0

end program test_zgemm
