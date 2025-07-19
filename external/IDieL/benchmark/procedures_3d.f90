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
!> Contains program to see procedure 3d for a real case (black phosphorous)

!> Program showcasing the 3D routines
program procedures_3d

    use idiel_constants, only: aip, i32, pi, twopi, zzero
    use idiel, only: idiel_t

    implicit none

    ! GW data
    integer(i32), parameter :: nomega = 32_i32  ! Number of frequencies
    integer(i32) :: iom

    ! Mesh data
    integer(i32), parameter :: ngrid(3) = [6_i32, 6_i32, 8_i32]

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
    a(:) = [3.13061478636775_aip, -9.89555085794508_aip, 0.00000000000000_aip]
    b(:) = [3.13061478636775_aip,  9.89555085794508_aip, 0.00000000000000_aip]
    c(:) = [0.00000000000000_aip,  0.00000000000000_aip, 8.26566207441072_aip]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [15_i32, 15_i32, 15_i32, 15_i32]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.8966_aip, 0.1034_aip, 0.0806_aip]
    redpos(:,2) = [0.6034_aip, 0.3966_aip, 0.5806_aip]
    redpos(:,3) = [0.3966_aip, 0.6034_aip, 0.4194_aip]
    redpos(:,4) = [0.1034_aip, 0.8966_aip, 0.9194_aip]

    ! Reading files
    call load_from_file('diel3d/head',  head)
    call load_from_file('diel3d/wingL', wingL)
    call load_from_file('diel3d/wingU', wingU)
    call load_from_file('diel3d/body',  body)

    !!!!!!!!!!!! Here it goes the main programm !!!!!!!!!!!!!!!!!!!

    ! Init common objects
#ifdef USE_SPGLIB
    call inv_diel%init_common(lattice, redpos, types, ngrid, 3_i32)
#else 
    crot(:,:,1 ) = reshape(real([  1.0,  0.0, 0.0,  0.0,  1.0,  0.0, 0.0, 0.0,  1.0 ], kind = aip ), [3,3])
    crot(:,:,2 ) = reshape(real([ -1.0,  0.0, 0.0,  0.0, -1.0,  0.0, 0.0, 0.0, -1.0 ], kind = aip ), [3,3])
    crot(:,:,3 ) = reshape(real([ -1.0,  0.0, 0.0,  0.0, -1.0,  0.0, 0.0, 0.0,  1.0 ], kind = aip ), [3,3])
    crot(:,:,4 ) = reshape(real([  1.0,  0.0, 0.0,  0.0,  1.0,  0.0, 0.0, 0.0, -1.0 ], kind = aip ), [3,3])
    crot(:,:,5 ) = reshape(real([  1.0,  0.0, 0.0,  0.0, -1.0,  0.0, 0.0, 0.0, -1.0 ], kind = aip ), [3,3])
    crot(:,:,6 ) = reshape(real([ -1.0,  0.0, 0.0,  0.0,  1.0,  0.0, 0.0, 0.0,  1.0 ], kind = aip ), [3,3])
    crot(:,:,7 ) = reshape(real([ -1.0,  0.0, 0.0,  0.0,  1.0,  0.0, 0.0, 0.0, -1.0 ], kind = aip ), [3,3])
    crot(:,:,8 ) = reshape(real([  1.0,  0.0, 0.0,  0.0, -1.0,  0.0, 0.0, 0.0,  1.0 ], kind = aip ), [3,3])

    call inv_diel%init_common(lattice, redpos, types, ngrid, 3_i32, nsym, crot)
#endif

    ! First do cycle for each frequency by considering the Hermitian matrix (actual case)
    do iom = 1, nomega
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_inversedielectric_3d(.true.)
    end do

    ! First do cycle for each frequency by considering the general case (not the case but we want to check performance)
    do iom = 1, nomega
        ! Invert body
        call inv_diel%invert_body(body(:,:,iom))
        ! Load the data to the worker
        call inv_diel%set_dielectric_blocks(head(:,:,iom), wingL(:,:,iom), wingU(:,:,iom))
        ! Compute the average
        call inv_diel%compute_anisotropic_avg_inversedielectric_3d(.false.)
    end do

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

end program procedures_3d
