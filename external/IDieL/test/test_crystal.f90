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
!> Tests the generated crystal structure as well as the symmetry and compares it with reference results
program test_crystal
    
    use idiel_constants,         only: i32, aip, pi, twopi
    use idiel_crystal_cell,      only: cell_t
    use idiel_crystal_symmetry,  only: symmetry_t

    implicit none

#ifdef USE_SINGLE_PRECISION
    real(aip), parameter    :: tolerance = 1.0e-6_aip
#else
    real(aip), parameter    :: tolerance = 1.0e-12_aip
#endif
    real(aip) :: rdiff

    ! Elements
    type(cell_t)     :: my_cell
    type(symmetry_t) :: my_symmetry

    ! cell_t elements
    integer(i32), parameter :: natoms = 4_i32
    real(aip) :: a(3) , b(3), c(3)
    real(aip) :: lattice(3,3) 
    real(aip), allocatable :: redpos(:,:)
    integer(i32), allocatable :: types(:)

    ! ref_values
    real(aip), parameter :: vuc_ref = 85.990737984502_aip
    real(aip) :: invlat(3,3)

    ! symmetry_t reference values
    integer(i32), parameter :: space_group = 186_i32
    integer(i32), parameter :: nsym = 12_i32
    integer(i32) :: rot_ref(3,3,12)
    real(aip)    :: crot_ref(3,3,12)

    ! Locals
    integer(i32) :: ii

    invlat(1,:) = [1.601845996473_aip, 0.924826217264_aip, 0.000000000000_aip]
    invlat(2,:) = [0.000000000000_aip, 1.849652434528_aip, 0.000000000000_aip]
    invlat(3,:) = [0.000000000000_aip, 0.000000000000_aip, 0.973592098042_aip]

    rot_ref(1,:,1) = [1_i32, 0_i32, 0_i32]
    rot_ref(2,:,1) = [0_i32, 1_i32, 0_i32]
    rot_ref(3,:,1) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,2) = [1_i32, -1_i32, 0_i32]
    rot_ref(2,:,2) = [1_i32, 0_i32, 0_i32]
    rot_ref(3,:,2) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,3) = [0_i32, -1_i32, 0_i32]
    rot_ref(2,:,3) = [1_i32, -1_i32, 0_i32]
    rot_ref(3,:,3) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,4) = [-1_i32, 0_i32, 0_i32]
    rot_ref(2,:,4) = [0_i32, -1_i32, 0_i32]
    rot_ref(3,:,4) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,5) = [-1_i32, 1_i32, 0_i32]
    rot_ref(2,:,5) = [-1_i32, 0_i32, 0_i32]
    rot_ref(3,:,5) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,6) = [0_i32, 1_i32, 0_i32]
    rot_ref(2,:,6) = [-1_i32, 1_i32, 0_i32]
    rot_ref(3,:,6) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,7) = [0_i32, 1_i32, 0_i32]
    rot_ref(2,:,7) = [1_i32, 0_i32, 0_i32]
    rot_ref(3,:,7) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,8) = [1_i32, 0_i32, 0_i32]
    rot_ref(2,:,8) = [1_i32, -1_i32, 0_i32]
    rot_ref(3,:,8) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,9) = [1_i32, -1_i32, 0_i32]
    rot_ref(2,:,9) = [0_i32, -1_i32, 0_i32]
    rot_ref(3,:,9) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,10) = [0_i32, -1_i32, 0_i32]
    rot_ref(2,:,10) = [-1_i32, 0_i32, 0_i32]
    rot_ref(3,:,10) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,11) = [-1_i32, 0_i32, 0_i32]
    rot_ref(2,:,11) = [-1_i32, 1_i32, 0_i32]
    rot_ref(3,:,11) = [0_i32, 0_i32, 1_i32]
    rot_ref(1,:,12) = [-1_i32, 1_i32, 0_i32]
    rot_ref(2,:,12) = [0_i32, 1_i32, 0_i32]
    rot_ref(3,:,12) = [0_i32, 0_i32, 1_i32]

    crot_ref(1,:,1) = [0.9999999999999999_aip, -1.2798967289084123e-16_aip, 0.0_aip]
    crot_ref(2,:,1) = [0.0_aip, 0.9999999999999999_aip, 0.0_aip]
    crot_ref(3,:,1) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,2) = [0.49999999999999994_aip, -0.8660254037844387_aip, 0.0_aip]
    crot_ref(2,:,2) = [0.8660254037844385_aip, 0.4999999999999999_aip, 0.0_aip]
    crot_ref(3,:,2) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,3) = [-0.49999999999999994_aip, -0.8660254037844386_aip, 0.0_aip]
    crot_ref(2,:,3) = [0.8660254037844385_aip, -0.5_aip, 0.0_aip]
    crot_ref(3,:,3) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,4) = [-0.9999999999999999_aip, 1.2798967289084123e-16_aip, 0.0_aip]
    crot_ref(2,:,4) = [0.0_aip, -0.9999999999999999_aip, 0.0_aip]
    crot_ref(3,:,4) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,5) = [-0.49999999999999994_aip, 0.8660254037844387_aip, 0.0_aip]
    crot_ref(2,:,5) = [-0.8660254037844385_aip, -0.4999999999999999_aip, 0.0_aip]
    crot_ref(3,:,5) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,6) = [0.49999999999999994_aip, 0.8660254037844386_aip, 0.0_aip]
    crot_ref(2,:,6) = [-0.8660254037844385_aip, 0.5_aip, 0.0_aip]
    crot_ref(3,:,6) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,7) = [-0.49999999999999994_aip, 0.8660254037844387_aip, 0.0_aip]
    crot_ref(2,:,7) = [0.8660254037844385_aip, 0.4999999999999999_aip, 0.0_aip]
    crot_ref(3,:,7) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,8) = [0.49999999999999994_aip, 0.8660254037844386_aip, 0.0_aip]
    crot_ref(2,:,8) = [0.8660254037844385_aip, -0.5_aip, 0.0_aip]
    crot_ref(3,:,8) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,9) = [0.9999999999999999_aip, -1.2798967289084123e-16_aip, 0.0_aip]
    crot_ref(2,:,9) = [0.0_aip, -0.9999999999999999_aip, 0.0_aip]
    crot_ref(3,:,9) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,10) = [0.49999999999999994_aip, -0.8660254037844387_aip, 0.0_aip]
    crot_ref(2,:,10) = [-0.8660254037844385_aip, -0.4999999999999999_aip, 0.0_aip]
    crot_ref(3,:,10) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,11) = [-0.49999999999999994_aip, -0.8660254037844386_aip, 0.0_aip]
    crot_ref(2,:,11) = [-0.8660254037844385_aip, 0.5_aip, 0.0_aip]
    crot_ref(3,:,11) = [0.0_aip, 0.0_aip, 1.0_aip]
    crot_ref(1,:,12) = [-0.9999999999999999_aip, 1.2798967289084123e-16_aip, 0.0_aip]
    crot_ref(2,:,12) = [0.0_aip, 0.9999999999999999_aip, 0.0_aip]
    crot_ref(3,:,12) = [0.0_aip, 0.0_aip, 1.0_aip]


    write(*,*) '[TEST : cell_t and symmetry_t]' 
    write(*,*) ' * kind : regression test against precomputed data'
    write(*,'(A, e20.13)') '  * tolerance : ', tolerance

    ! Data for WZ ZnSe

    ! Lattice vectors
    a(:) = 3.9224652813161520_aip * [1.0_aip,  0.0_aip,               0.0_aip]
    b(:) = 3.9224652813161520_aip * [-0.5_aip, sqrt(3.0_aip)/2.0_aip, 0.0_aip]
    c(:) = 3.9224652813161520_aip * [0.0_aip,  0.0_aip,               1.6452947797075332_aip]

    lattice(1,:) = a(:)
    lattice(2,:) = b(:)
    lattice(3,:) = c(:)

    ! Atomic types (they might differ that the atomic number)
    allocate(types(natoms))
    types(:) = [30_i32, 30_i32, 34_i32, 34_i32]

    ! Reduced positions
    allocate(redpos(3,natoms))
    redpos(:,1) = [0.3333333333333357_aip, 0.6666666666666643_aip, 0.0000000000000000_aip]
    redpos(:,2) = [0.6666666666666643_aip, 0.3333333333333357_aip, 0.5000000000000000_aip]
    redpos(:,3) = [0.3333333333333357_aip, 0.6666666666666643_aip, 0.3736232276372607_aip]
    redpos(:,4) = [0.6666666666666643_aip, 0.3333333333333357_aip, 0.8736232276372607_aip]

    ! Initialize the crystal structure
    call my_cell%initialize(lattice, redpos, types)

    ! Check computed values
    rdiff = abs(my_cell%vuc - vuc_ref)/vuc_ref
    write(*,'(A, e20.13)')  '  * Regression (vuc) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : cell_t%vuc : PASSED]'
    else
        write(*,*)  '[TEST : cell_t%vuc : FAILED]'
        stop 1
    end if

    rdiff = relative_diference_m(my_cell%rlattice,invlat)
    write(*,'(A, e20.13)')  '  * Regression (rlattice) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : cell_t%rlattice : PASSED]'
    else
        write(*,*)  '[TEST : cell_t%rlattice : FAILED]'
        stop 1
    end if
#ifdef USE_SPGLIB
    ! Computing symmetry stuff
    call my_symmetry%initialize(my_cell)

    ! Check the space group
    if ( my_symmetry%spgid .eq. space_group) then
        write(*,*)  '[TEST : symmetry_t%spgid : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%spgid : FAILED]'
        stop 1
    end if

    !Check the rotations
    if ( my_symmetry%nsym .eq. nsym) then
        write(*,*)  '[TEST : symmetry_t%nsym : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%nsym : FAILED]'
        stop 1
    end if

    rdiff = relative_diference_i3(my_symmetry%rot, rot_ref)
    write(*,'(A, e20.13)')  '  * Regression (rotations) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : symmetry_t%rot : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%rot : FAILED]'
        stop 1
    end if

    rdiff = relative_diference_r3(my_symmetry%crot, crot_ref)
    write(*,'(A, e20.13)')  '  * Regression (Cartesian rotations) result (relative difference): ', rdiff
    if ( rdiff .lt. tolerance) then
        write(*,*)  '[TEST : symmetry_t%crot : PASSED]'
    else
        write(*,*)  '[TEST : symmetry_t%crot : FAILED]'
        stop 1
    end if
#endif
    stop 0
contains
        !> Computes the Frobenious norm of a matrix
        pure function relative_diference_m(array,reference) result(rd)
            !> The array for which the Frobenious norm will be computed
            real(aip), intent(in) :: array(:,:)
            real(aip), intent(in) :: reference(:,:)
            real(aip) :: rd
            rd = sqrt(real(sum( (array-reference)**2 ))) / sqrt(real(sum((reference)**2 )))
        end function relative_diference_m
        
        !> Computes the Frobenious norm of a matrix
        pure function relative_diference_i3(array,reference) result(rd)
            !> The array for which the Frobenious norm will be computed
            integer(i32), intent(in) :: array(:,:,:)
            integer(i32), intent(in) :: reference(:,:,:)
            real(aip) :: rd
            rd = sqrt(real(sum( (array-reference)**2 ))) / sqrt(real(sum((reference)**2 )))
        end function relative_diference_i3

        !> Computes the Frobenious norm of a matrix
        pure function relative_diference_r3(array,reference) result(rd)
            !> The array for which the Frobenious norm will be computed
            real(aip), intent(in) :: array(:,:,:)
            real(aip), intent(in) :: reference(:,:,:)
            real(aip) :: rd
            rd = sqrt(real(sum( (array-reference)**2 ))) / sqrt(real(sum((reference)**2 )))
        end function relative_diference_r3

end program test_crystal
