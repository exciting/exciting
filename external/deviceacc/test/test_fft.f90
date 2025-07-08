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
!> This file contains a test of the FFT using this library
program test_fft

    use iso_fortran_env,         only: i32=>int32, r64=>real64 
    use iso_c_binding
    use m_device_world_t,        only: device_world_t
    use m_fft_device,            only: fft_device_t
#if defined(DEVICEOFFLOAD)
    use mpi
#endif

    implicit none

    type(fft_device_t)        :: fft
    type(device_world_t)      :: device_world
    integer(i32)              :: mpi_world, err
    real(r64)                 :: res(2)

    integer(i32), parameter :: m=768, n=512
    complex(r64), target :: a(m,n) 
    type(c_ptr) :: da 

#if defined(DEVICEOFFLOAD)
    ! Init MPI world
    call mpi_init(err)
    mpi_world = MPI_COMM_WORLD
#endif
    ! Init GPU stuff
    call device_world%init(mpi_world)

    ! Solve FFT
    a(:,:) = cmplx(1.0, 0.0, kind=r64) 
    
    ! Update a to the gpu
    call device_world%register%alloc("a", c_sizeof(a), device_world%get_device())
    call device_world%register%assoc("a", c_loc(a))
    call device_world%register%to_device("a")

    da = device_world%register%device_ptr("a")

    ! Init FFT 
    call fft%initialize(shape(a), -1, da, .true., device_world)

    ! Perform the FFT
    call fft%execute(da, .false., .true.)

    ! Clean FFT
    call fft%delete()

    ! Get from the device
    call device_world%register%from_device("a")

    ! End device and MPI, the register destructor frees all device memory and deassociates it from the host 
    ! so individual calls to deassoc and remove is not necessary
    call device_world%finish()
#if defined(DEVICEOFFLOAD)
    call mpi_finalize(err)
#endif

    ! Check results
    write(*,*) '[TEST : fft]' 
    res(1) = maxval( real(a) ) - sum( real(a) )
    res(2) = maxval( aimag(a) )

    if (maxval( res ) .lt. 1.0e-12_r64) then
        write(*,*)  '[TEST : fft : PASSED]'
    else
        write(*,*)  '[TEST : fft : FAILED]'
        stop 1
    end if

    stop 0

end program test_fft
