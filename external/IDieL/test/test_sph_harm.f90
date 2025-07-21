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
!> Tests the spherical harmonics generation
program test_sph_harm

    use idiel_constants, only: aip, i32
    use idiel_sph_quadrature, only: compute_angular_mesh_lebedev_21
    use idiel_spherical_harmonics,  only: sph_harm

    implicit none 

    ! Spherical harmonics
    integer(i32), parameter :: lmax = 5_i32
    integer(i32), parameter :: nsph = (lmax+1)**2
    complex(aip), allocatable :: ylm(:,:)
    integer :: i
    real(aip) :: real_ylm(nsph), imag_ylm(nsph)
    complex(aip) :: ylm_ref(170, nsph)

    ! Tolerance
#ifdef USE_SINGLE_PRECISION
    real(aip), parameter :: tolerance = 1.0e-6_aip
#else
    real(aip), parameter :: tolerance = 1.0e-12_aip
#endif
    real(aip) :: rdiff

    ! Weights and points
    integer(i32), parameter :: np = 170_i32
    real(aip), allocatable  :: ang_mesh(:,:)
    real(aip), allocatable  :: xyz_mesh(:,:)
    real(aip), allocatable  :: weights(:)

    ! Print info
    write(*,*) '[TEST : sph_harm]' 
    write(*,*) ' * kind : regression test against precomputed data'

    ! Compute the integral points and weights
    call compute_angular_mesh_lebedev_21(ang_mesh, weights, xyz_mesh)

    ! Compute the spherical harmonics
    call sph_harm(lmax, ang_mesh, ylm)

    ! Read data to compare to
    open(unit=101, file='real_part.txt', status='old', action='read')
    open(unit=102, file='imag_part.txt', status='old', action='read')
    do i = 1, np
        read(101, *) real_ylm(:)
        read(102, *) imag_ylm(:)
        ylm_ref(i,:) = cmplx( real_ylm, imag_ylm, aip)
    end do
    close(101)
    close(102)

    rdiff = sum(abs(ylm - ylm_ref))/sum(abs(ylm_ref))
    write(*,'(A, e20.13)')  '  * Regression spherical harmonics result (relative difference): ', rdiff
    if ( rdiff < tolerance ) then
         write(*,*)  '[TEST : sph_harm: PASSED]'
    else
         write(*,*)  '[TEST : sph_harm: FAILED]'
         stop 1
    end if

    stop 0

end program test_sph_harm
