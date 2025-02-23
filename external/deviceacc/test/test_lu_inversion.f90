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
!> This file contains a test of the LU inversion using GPU accelerated routines
program test_lu_inversion
    
    use iso_fortran_env,         only: i32=>int32, r64=>real64 
    use iso_c_binding
    use m_device_world_t,        only: device_world_t
    use device_linalg_common_interface, only: get_zgetri_nb_gpu,zgetrf_gpu, zgetri_gpu
#if defined(DEVICEOFFLOAD)
    use mpi
#endif

    implicit none

    integer(i32), parameter   :: n = 8000_i32
    complex(r64), parameter   :: zzero = cmplx(0.0, 0.0, kind=r64)
    complex(r64), target, allocatable :: A(:,:), inverse_A(:,:)
    type(device_world_t)      :: device_world
    integer(i32)              :: i
    real(r64)                 :: rdiff
    real(r64), parameter      :: tolerance = 1.0e-12_r64
    integer(i32)              :: mpi_world, err
    
    ! LAPACK
    integer(i32) :: info, lwork, nb
    integer(i32), allocatable :: ipiv(:)
    complex(r64), target, allocatable :: work(:)
    type(c_ptr) :: dwork_ptr, dA_ptr

#if defined(DEVICEOFFLOAD)
    ! Init MPI world
    call mpi_init(err)
    mpi_world = MPI_COMM_WORLD
#endif
    
    ! Allocating the matrix A
    allocate(A(n,n), source=zzero)

    do i = 1, n
        A(i,i) = cmplx(i, 0, r64)
    end do

    ! Init GPU stuff
    call device_world%init(mpi_world)
    
    ! Doing the inverse
    nb = get_zgetri_nb_gpu(n)
    lwork = nb * size(A,1)
    allocate(ipiv(size(A,1)))

    allocate(inverse_A, source=A)

    !Allocate in device inverse_A, associate to host inverse_A, and copy from host to device
    call device_world%register%alloc("inverse_A", size(inverse_A) * c_sizeof(zzero), device_world%get_device())
    call device_world%register%assoc("inverse_A", C_loc(inverse_A))
    call device_world%register%to_device("inverse_A")
    ! Allocate in device work array
    call device_world%register%alloc("work", lwork*c_sizeof(zzero), device_world%get_device())

    ! Obtain device pointers
    dA_ptr    = device_world%register%device_ptr("inverse_A")
    dwork_ptr = device_world%register%device_ptr("work")
    
    call zgetrf_gpu(n, n, dA_ptr, n, ipiv, info, device_world)
    if (info /= 0) error stop "inverse_complex_LU: error calling zgetrf gpu"
    call zgetri_gpu(n, dA_ptr, n, ipiv, dwork_ptr, lwork, info, device_world)
    if (info /= 0) error stop "inverse_complex_LU: error calling zgetri gpu"

    call device_world%register%from_device("inverse_A")
    call device_world%register%remove("inverse_A")
    call device_world%register%remove("work")

    deallocate(ipiv)

    call device_world%synchronize()

    ! Checking the results
    write(*,*) '[TEST : test_lu_inversion]' 
    do i = 1, n
        write(*,*) 1.0_r64 / A(i,i), inverse_A(i,i)
        rdiff = abs(1.0_r64 / A(i,i) - inverse_A(i,i))/ abs(1.0_r64 / A(i,i))
        write(*,'(A,I0,A, E16.7)')  '  * Regression (inverse,', i ,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
            write(*,*)  '[TEST : test_lu_inversion (inverse,', i,'): PASSED]'
        else
            write(*,*)  '[TEST : test_lu_inversion (inverse,', i,'): FAILED]'
            stop 1
        end if
    end do

    call device_world%finish()
#if defined(DEVICEOFFLOAD)
    call mpi_finalize(err)
#endif

    stop 0

end program test_lu_inversion
