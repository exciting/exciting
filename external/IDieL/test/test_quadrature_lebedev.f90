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
!> Tests the spherical quadrature using Lebedev-mesh works by comparing to objective functions

program test_quadrature_lebedev
    
    use idiel_constants, only: aip, i32, pi, twopi
    use idiel_sph_quadrature, only: compute_angular_mesh_lebedev_131, compute_angular_mesh_lebedev_41, &
                                    compute_angular_mesh_lebedev_21

    implicit none

    ! Tolerance
    real(aip), parameter :: tolerance = 1.0e-5_aip
    ! The alpha parameter
    real(aip), parameter :: alpha = 12.0_aip
    ! The exact integrals of test functions 
    real(aip), parameter :: I1 = 216_aip * pi / 35.0_aip
    real(aip), parameter :: I2 = 6.6961822200736179523_aip
    real(aip), parameter :: I3 = 4.0_aip * pi / alpha
    real(aip), parameter :: I4 = 4.0_aip * pi / alpha
    real(aip), parameter :: I5 = 4.0_aip * pi / alpha

    ! Weights and points
    real(aip), allocatable :: ang_mesh(:,:)
    real(aip), allocatable :: xyz_mesh(:,:)
    real(aip), allocatable :: weights(:)

    ! Storage
    real(aip), allocatable :: f(:)
    real(aip) :: numerical_integral, rdiff

    write(*,*) '[TEST : compute_ang_lebedev]' 
    write(*,*) ' * kind : regression test against precomputed data'

    ! Compute the integral points and weights
    call compute_angular_mesh_lebedev_131(ang_mesh, weights, xyz_mesh)

    ! Compute integral 1 and check
    f = f1(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I1 ) / I1
    write(*,'(A, e20.13)')  '  * Regression Lebedev 131 (f1) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f1: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f1: FAILED]'
        stop 1
    end if

    !Compute integral 2 and check (this one has lower convergence for all quadratures, thus we reduce convergence criteria)
    f = f2(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I2 ) / I2
    write(*,'(A, e20.13)')  '  * Regression Lebedev 131 (f2) result (relative difference): ', rdiff
    if (rdiff .lt. 1.0e+3_aip * tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f2: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f2: FAILED]'
        stop 1
    end if

    ! Compute integral 3 and check
    f = f3(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I3 ) / I3
    write(*,'(A, e20.13)')  '  * Regression Lebedev 131 (f3) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f3: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f3: FAILED]'
        stop 1
    end if

    ! Compute integral 4 and check
    f = f4(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I4 ) / I4
    write(*,'(A, e20.13)')  '  * Regression Lebedev 131 (f4) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f4: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f4: FAILED]'
        stop 1
    end if

    ! Compute integral 5 and check
    f = f5(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I5 ) / I5
    write(*,'(A, e20.13)')  '  * Regression Lebedev 131 (f5) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f5: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f5: FAILED]'
        stop 1
    end if

    deallocate(xyz_mesh, weights, ang_mesh)

    ! Compute the integral points and weights
    call compute_angular_mesh_lebedev_21(ang_mesh, weights, xyz_mesh)

    ! Compute integral 1 and check
    f = f1(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I1 ) / I1
    write(*,'(A, e20.13)')  '  * Regression Lebedev 21 (f1) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f1: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f1: FAILED]'
        stop 1
    end if

    !Compute integral 2 and check (this one has lower convergence for all quadratures, thus we reduce convergence criteria)
    f = f2(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I2 ) / I2
    write(*,'(A, e20.13)')  '  * Regression Lebedev 21 (f2) result (relative difference): ', rdiff
    if (rdiff .lt. 1.0e+7_aip * tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f2: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f2: FAILED]'
        stop 1
    end if

    ! Compute integral 3 and check
    f = f3(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I3 ) / I3
    write(*,'(A, e20.13)')  '  * Regression Lebedev 21 (f3) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f3: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f3: FAILED]'
        stop 1
    end if

    ! Compute integral 4 and check
    f = f4(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I4 ) / I4
    write(*,'(A, e20.13)')  '  * Regression Lebedev 21 (f4) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f4: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f4: FAILED]'
        stop 1
    end if

    ! Compute integral 5 and check
    f = f5(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I5 ) / I5
    write(*,'(A, e20.13)')  '  * Regression Lebedev 21 (f5) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f5: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f5: FAILED]'
        stop 1
    end if

    deallocate(xyz_mesh, weights, ang_mesh)

    ! Compute the integral points and weights
    call compute_angular_mesh_lebedev_41(ang_mesh, weights, xyz_mesh)

    ! Compute integral 1 and check
    f = f1(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I1 ) / I1
    write(*,'(A, e20.13)')  '  * Regression Lebedev 41 (f1) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f1: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f1: FAILED]'
        stop 1
    end if

    !Compute integral 2 and check (this one has lower convergence for all quadratures, thus we reduce convergence criteria)
    f = f2(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I2 ) / I2
    write(*,'(A, e20.13)')  '  * Regression Lebedev 41 (f2) result (relative difference): ', rdiff
    if (rdiff .lt. 1.0e+7_aip * tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f2: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f2: FAILED]'
        stop 1
    end if

    ! Compute integral 3 and check
    f = f3(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I3 ) / I3
    write(*,'(A, e20.13)')  '  * Regression Lebedev 41 (f3) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f3: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f3: FAILED]'
        stop 1
    end if

    ! Compute integral 4 and check
    f = f4(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I4 ) / I4
    write(*,'(A, e20.13)')  '  * Regression Lebedev 41 (f4) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f4: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f4: FAILED]'
        stop 1
    end if

    ! Compute integral 5 and check
    f = f5(xyz_mesh)
    numerical_integral = sum(weights * f)
    rdiff = abs( numerical_integral - I5 ) / I5
    write(*,'(A, e20.13)')  '  * Regression Lebedev 41 (f5) result (relative difference): ', rdiff
    if (rdiff .lt. tolerance) then
        write(*,*)  '[TEST : compute_ang_lebedev f5: PASSED]'
    else
        write(*,*)  '[TEST : compute_ang_lebedev f5: FAILED]'
        stop 1
    end if


    stop 0
contains

    !> Test function 1
    pure function f1(r)
        !> Mesh in which to compute
        real(aip), intent(in)  :: r(:,:)
        !> Result
        real(aip), allocatable :: f1(:)
        !> Elements
        real(aip), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f1 = 1.0_aip + x + y**2 + y * x**2 + x**4 + y**5 + (x*y*z)**2 

    end function f1

    !> Test function 2
    pure function f2(r)
        !> Mesh in which to compute
        real(aip), intent(in)  :: r(:,:)
        !> Result
        real(aip), allocatable :: f2(:)
        !> Elements
        real(aip), allocatable  :: x(:), y(:), z(:)
        real(aip), allocatable, dimension(:) :: exp1, exp2, exp3, exp4

        x = 9.0_aip * r(:,1)
        y = 9.0_aip * r(:,2)
        z = 9.0_aip * r(:,3)

        exp1 = 0.75 * exp( -0.25 * (x - 2.0)**2)     * exp(-0.25*(y - 2.0)**2)  * exp(-0.25*(z - 2.0)**2)
        exp2 = 0.75 * exp( -1.0/49.0 * (x + 1.0)**2) * exp(-0.10*(y + 1.0))     * exp(-0.10*(z + 1.0))
        exp3 = 0.50 * exp( -0.25 * (x - 7.0)**2)     * exp(-0.25*(y - 3.0)**2)  * exp(-0.25*(z - 5.0)**2)
        exp4 = 0.20 * exp( -(x - 4.0)**2)            * exp(-(y - 7.0)**2)       * exp(-(z - 5.0)**2)

        f2 = exp1 + exp2 + exp3 - exp4

    end function f2

    !> Test function 3
    pure function f3(r)
        !> Mesh in which to compute
        real(aip), intent(in) :: r(:,:)
        !> Result
        real(aip), allocatable :: f3(:)
        !> Elements
        real(aip), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f3 = (1.0_aip + tanh(-alpha*(x+y-z)))/alpha

    end function f3

    !> Sgn function
    pure function sgn(value)
        
        real(aip), intent(in) :: value(:)
        real(aip), allocatable :: sgn(:)
        integer(i32) :: i
        
        allocate(sgn,source=value)

        do i = 1, size(value,1)
            if (sgn(i) .lt. 0.0_aip) then
                sgn(i) = -1.0_aip
            else if (sgn(i) .eq. 0.0_aip) then
                sgn(i) =  0.0_aip
            else 
                sgn(i) =  1.0_aip
            end if
        end do

    end function sgn

    !> Test function 4
    pure function f4(r)
        !> Mesh in which to compute
        real(aip), intent(in) :: r(:,:)
        !> Result
        real(aip), allocatable :: f4(:)
        !> Elements
        real(aip), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f4 = (1.0_aip - sgn(x+y-z))/alpha
        
    end function f4

    !> Test function 4
    pure function f5(r)
        !> Mesh in which to compute
        real(aip), intent(in) :: r(:,:)
        !> Result
        real(aip), allocatable :: f5(:)
        !> Elements
        real(aip), allocatable :: x(:), y(:), z(:)

        x = r(:,1)
        y = r(:,2)
        z = r(:,3)

        f5 = (1.0_aip - sgn(real(pi,aip)*x + y))/alpha
        
    end function f5
    

end program test_quadrature_lebedev
