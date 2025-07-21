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
!> Contains program to see procedure 2d for a real case (phosphorene)

!> Program showcasing the 2D routines
program test_2d

    use idiel_constants, only: aip, i32, pi, twopi, zzero
    use idiel, only: idiel_t

    implicit none

    ! GW data
    integer(i32), parameter :: nomega = 32_i32  ! Number of frequencies
    integer(i32) :: iom

    ! Mesh data
    integer(i32), parameter :: ngrid(3) = [18_i32, 14_i32, 1_i32]

    ! Silicon crystal data
    integer(i32), parameter :: natoms = 4_i32
    real(aip) :: a(3), b(3), c(3), lattice(3,3)
    real(aip), allocatable :: redpos(:,:)
    integer(i32), allocatable :: types(:)

    ! Dielectric data
    complex(aip), allocatable :: head(:,:,:), wingL(:,:,:)
    complex(aip), allocatable :: wingU(:,:,:), body(:,:,:)

    ! Computation object
    type(idiel_t) :: inv_diel

    ! In case that no SPG we set only Identity (i.e. no symmetry)
    integer(i32) :: nsym = 8
    real(aip)    :: crot(3,3,8)

    ! Lattice vectors
    a(:) = [6.26122957273549_aip, 0.00000000000000_aip, 0.000000000000000_aip]
    b(:) = [0.00000000000000_aip, 8.26566207441072_aip, 0.000000000000000_aip]
    c(:) = [0.00000000000000_aip, 0.00000000000000_aip, 19.79110171589015_aip]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [15_i32, 15_i32, 15_i32, 15_i32]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.50000000000000_aip, 0.08060000000000_aip, 0.60340000000000_aip]
    redpos(:,2) = [0.00000000000000_aip, 0.58060000000000_aip, 0.39660000000000_aip]
    redpos(:,3) = [0.00000000000000_aip, 0.41940000000000_aip, 0.60340000000000_aip]
    redpos(:,4) = [0.50000000000000_aip, 0.91940000000000_aip, 0.39660000000000_aip]

    ! Reading files
    call load_from_file('diel2d/head',  head)
    call load_from_file('diel2d/wingL', wingL)
    call load_from_file('diel2d/wingU', wingU)
    call load_from_file('diel2d/body',  body)

    !!!!!!!!!!!! Here it goes the main programm !!!!!!!!!!!!!!!!!!!
#ifdef USE_SPGLIB
    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i32)
#else 
    crot(:,:,1 ) = reshape(real([  1.0,  0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:,2 ) = reshape(real([ -1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:,3 ) = reshape(real([ -1.0,  0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:,4 ) = reshape(real([  1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:,5 ) = reshape(real([  1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0 ], kind = aip), [3,3])
    crot(:,:,6 ) = reshape(real([ -1.0,  0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:,7 ) = reshape(real([ -1.0,  0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0,  1.0 ], kind = aip), [3,3])
    crot(:,:,8 ) = reshape(real([  1.0,  0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, -1.0 ], kind = aip), [3,3])

    call inv_diel%init_common(lattice, redpos, types, ngrid, 2_i32, nsym, crot)
#endif

    ! Do the average (Hermitian)
    do iom = 1, nomega
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.true.)
    end do

    ! Do the average (General case)
    do iom = 1, nomega
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_scrcoulomb_2d(.false.)
    end do

contains 

    !> Process a file and load to an array
    !> @param[in] fname - the file name   
    !> @param[in] data_shape - the shape of the data to reads
    !> @param[out] data - the data
    subroutine load_from_file(fname, data)

        character(len=*), intent(in) :: fname
        complex(aip), allocatable, intent(out) :: data(:,:,:)

        integer :: fin
        real(aip), allocatable  :: dreal(:,:,:), dimag(:,:,:)
        integer(i32) :: data_shape(3)
        
        open(file=fname, newunit=fin, status='old', action='read')

        read(fin, *) data_shape

        allocate(data(data_shape(1),data_shape(2),data_shape(3)))
        allocate(dreal(data_shape(1),data_shape(2),data_shape(3)))
        allocate(dimag, mold=dreal)
        read(fin, *) dreal
        read(fin, *) dimag
        data = cmplx(dreal, dimag, aip) 
        
        close(fin)

    end subroutine load_from_file

end program test_2d
