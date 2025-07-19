! Copyright (C) 2020-2024 GreenX library
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
!> Contains procedures to test averages in 2D cases (MoS2)

!> This module contains the procedures for the test of averaging of the dielectric matrix in 2D cases
program test_2d

    use idiel_constants, only: aip, i32, pi, twopi, zzero
    use idiel, only: idiel_t


#if defined(DEVICEOFFLOAD)
    use mpi
#endif

    implicit none

    ! GW data
    integer, parameter :: mbsize = 438 ! Mixed basis size
    integer, parameter :: nomega = 12  ! Number of frequencies
    integer(i32) :: iom

    ! Mesh data
    integer(i32), parameter :: ngrid(3) = [2,2,1]

    ! Silicon crystal data
    integer(i32), parameter :: natoms = 3_i32
    real(aip) :: a(3), b(3), c(3), lattice(3,3)
    real(aip), allocatable :: redpos(:,:)
    integer(i32), allocatable :: types(:)

    ! Dielectric data
    complex(aip), allocatable :: head(:,:,:), wingL(:,:,:)
    complex(aip), allocatable :: wingU(:,:,:), body(:,:,:)

    ! Reference data to compare
    real(aip) :: head_ref(4) = [ -0.27165343E+03_aip, -0.97552958E+02_aip, -0.49406608E+02_aip, -0.46006706E+01_aip]
    complex(aip) :: body_ref(4) = [ (-1.98341339762901914E-002_aip,-3.38140317306334658E-019_aip), &
                                    (-9.59818027259629059E-003_aip,-4.06834974195416980E-020_aip), &
                                    (-6.16528641837310598E-003_aip, 1.23527194124926345E-019_aip), &
                                    (-9.44568563628567226E-004_aip,-2.00846029587093419E-021_aip)] 
    real(aip) :: rdiff
#ifdef USE_SINGLE_PRECISION
    real(aip), parameter    :: tolerance = 0.07_aip
#else
    real(aip), parameter    :: tolerance = 0.05_aip
#endif
    ! Computation object
    type(idiel_t) :: inv_diel

    ! In case that no SPG we set only Identity (i.e. no symmetry)
    integer(i32) :: nsym = 12
    real(aip)    :: crot(3,3,12)

    ! MPI
    integer(i32) :: mpi_world, err


#if defined(DEVICEOFFLOAD)
    ! Init MPI world
    call mpi_init(err)
    mpi_world = MPI_COMM_WORLD
#endif

    ! Lattice vectors
    a(:) = [0.00000000_aip, 6.020289060_aip, 0.00000000_aip]
    b(:) = [5.213723260_aip, 3.010144530_aip, 0.00000000_aip]
    c(:) = [0.0_aip,  0.0_aip, 20.0_aip]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [1_i32, 2_i32, 2_i32]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.00_aip, 0.00_aip, 0.00_aip]
    redpos(:,2) = [2.0_aip/3.0_aip , 2.0_aip/3.0_aip,  0.14764992_aip]
    redpos(:,3) = [2.0_aip/3.0_aip , 2.0_aip/3.0_aip, -0.14764992_aip]

    ! Init common objects
#ifdef USE_SPGLIB
#ifdef DEVICEOFFLOAD
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i32, host_world=mpi_world)
#else
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i32)
#endif
#else
    crot(:,:, 1 ) = reshape(real([   1.00,        0.0000000000000000,        0.0,        0.0000000000000000,        1.00, 0.0, 0.0 , 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:, 2 ) = reshape(real([ -0.500,      -0.86602540446306686,        0.0,       0.86602540310581044,      -0.500, 0.0, 0.0 , 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:, 3 ) = reshape(real([ -0.500,       0.86602540446306686,        0.0,      -0.86602540310581044,      -0.500, 0.0, 0.0 , 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:, 4 ) = reshape(real([   1.00,        0.0000000000000000,        0.0,        0.0000000000000000,        1.00, 0.0, 0.0 , 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:, 5 ) = reshape(real([ -0.500,      -0.86602540446306686,        0.0,       0.86602540310581044,      -0.500, 0.0, 0.0 , 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:, 6 ) = reshape(real([ -0.500,       0.86602540446306686,        0.0,      -0.86602540310581044,      -0.500, 0.0, 0.0 , 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:, 7 ) = reshape(real([ -0.500,      -0.86602540446306686,        0.0,      -0.86602540310581044,       0.500, 0.0, 0.0 , 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:, 8 ) = reshape(real([   1.00,        0.0000000000000000,        0.0,        0.0000000000000000,       -1.00, 0.0, 0.0 , 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:, 9 ) = reshape(real([ -0.500,       0.86602540446306686,        0.0,       0.86602540310581044,       0.500, 0.0, 0.0 , 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:,10 ) = reshape(real([ -0.500,      -0.86602540446306686,        0.0,      -0.86602540310581044,       0.500, 0.0, 0.0 , 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:,11 ) = reshape(real([   1.00,        0.0000000000000000,        0.0,        0.0000000000000000,       -1.00, 0.0, 0.0 , 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:,12 ) = reshape(real([ -0.500,       0.86602540446306686,        0.0,       0.86602540310581044,       0.500, 0.0, 0.0 , 0.0,  1.0 ], kind = aip), [3,3])
#ifdef DEVICEOFFLOAD
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i32, nsym, crot, host_world=mpi_world)
#else
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i32, nsym, crot)
#endif
#endif

    ! Read the dielectric data from previous G0W0 run
    call load_from_file('head.2d.dat',  head)
    call load_from_file('wings.2d.dat', wingL, wingU)
    call load_from_file('body.2d.dat',  body)

    ! Do the average
    write(*,*) '[TEST: idiel_t (2d)]'
    do iom = 1, size(head,3)
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.true.)

        ! Check the head
        rdiff = abs(inv_diel%idiel_head - head_ref(iom))/abs(head_ref(iom))
        write(*,'(A,I2,A,e20.13)')  '  * Regression (HEAD,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
             write(*,*)  '[TEST: idiel_t (2d) (HEAD,',iom,'): PASSED]', abs(inv_diel%idiel_head), abs(head_ref(iom))
        else
             write(*,*)  '[TEST: idiel_t (2d) (HEAD,',iom,'): FAILED]', abs(inv_diel%idiel_head), abs(head_ref(iom))
             stop 1
        end if

        ! Check the body (only first element)
        rdiff = abs(inv_diel%idiel_body(1,1) - body_ref(iom))/abs(body_ref(iom))
        write(*,'(A,I2,A,e20.13)')  '  * Regression (BODY,',iom,') result (relative difference): ', rdiff
        if ( rdiff .lt. tolerance) then
             write(*,*)  '[TEST: idiel_t (2d) (BODY,',iom,'): PASSED]', abs(inv_diel%idiel_body(1,1)), abs(body_ref(iom))
        else
             write(*,*)  '[TEST: idiel_t (2d) (BODY,',iom,'): FAILED]', abs(inv_diel%idiel_body(1,1)), abs(body_ref(iom))
             stop 1
        end if

    end do

#if defined(DEVICEOFFLOAD)
    call mpi_finalize(err)
#endif

    stop 0

contains 

    !> Process a file and load to anray
    !> @param[in] fname - the file name   
    !> @param[in] data_shape - the shape of the data to reads
    !> @param[out] data - the data
    subroutine load_from_file(fname, data, data2)

        character(len=*), intent(in) :: fname
        complex(aip), allocatable, intent(out) :: data(:,:,:)
        complex(aip), allocatable, optional, intent(out) :: data2(:,:,:)

        integer :: fin
        real(aip), allocatable  :: dreal(:,:,:), dimag(:,:,:)
        integer(i32) :: data_shape(3)
        
        open(file=fname, newunit=fin, status='old', action='read')

        read(fin, *) data_shape

        if (present(data2)) then 
            allocate(data(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dreal(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dimag, mold=dreal)
            allocate(data2, mold=data)
            read(fin, *) dreal
            read(fin, *) dimag
            data = cmplx(dreal, dimag, aip) 
            read(fin, *) dreal
            read(fin, *) dimag
            data2 = cmplx(dreal, dimag, aip) 
        else
            allocate(data(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dreal(data_shape(1),data_shape(2),data_shape(3)))
            allocate(dimag, mold=dreal)
            read(fin, *) dreal
            read(fin, *) dimag
            data = cmplx(dreal, dimag, aip) 
        end if

        close(fin)

    end subroutine load_from_file

end program test_2d
