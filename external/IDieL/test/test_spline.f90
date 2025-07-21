! Copyright 2023 EXCITING developers
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
!> Tests the computation of 1D splines
program test_spline_1d

    use idiel_constants, only: r64, i32, twopi, pi
    use idiel_cubic_spline

    implicit none
    
    integer(i32), parameter :: n = 60_i32
    real(r64), parameter :: tolerance  = 5.0e-3_r64
    real(r64), parameter :: rtolerance = 5.0e-5_r64
    real(r64), parameter :: ana_integral = 2.35718_r64 ! from 0.5 to pi
    real(r64), parameter :: ana_derivative = -0.550858_r64 ! for 0.3*pi
    complex(r64) :: y(n), yref(10*n), yint
    real(r64)    :: x(n), xdenser(10*n)
    integer(i32) ::  i
    real(r64) :: dx, dx2, diff, rdiff, ics, dcs
    type(cubic_spline_t) :: cs

    ! Print info
    write(*,*) '[TEST : cubic_spline_t]' 
    write(*,*) ' * kind : regression test against precomputed data'
    write(*,*) ' * Test function : y(x) = cos(x)**2 * x '  

    dx  = twopi / (n-1)
    dx2 = twopi / (10*n-1) 

    do i = 1, n
        x(i) = (i - 1) * dx
        write(*,*) x(i)
        y(i) = cmplx(cos(x(i))**2 * x(i) , 0.0, r64) 
    end do

    do i = 1, 10*n
        xdenser(i) = (i - 1) * dx2
        yref(i) = cmplx(cos(xdenser(i))**2 * xdenser(i) , 0.0, r64) 
    end do

    call cs%initialize(x, y)
    

    write(*,'(A, e20.13)')  '  * Compare vs exact function:  cubic_spline_t%interpolate'
    do i = 1, 9*n
        yint = cs%interpolate(xdenser(i))
        if (real(yref(i)) == 0) cycle
        diff = abs(real(yint)-real(yref(i)))
        if ( diff > tolerance) then
           write(*,*)  '[TEST : cubic_spline_t%interpolate: FAILED] , diff : ', diff
           stop 1
        else 
           write(*,*)  '[TEST : cubic_spline_t%interpolate: PASSED] , diff : ', diff
        end if 
    end do

    write(*,'(A, e20.13)')  '  * Compare vs exact result:  cubic_spline_t%derivative'
    dcs = real(cs%derivative(0.3_r64*pi))
    rdiff = abs(dcs - ana_derivative)/abs(ana_derivative)
    if ( rdiff > 100*rtolerance) then
        write(*,*)  '[TEST : cubic_spline_t%derivative: FAILED] , rdiff : ', rdiff
        stop 1
     else
        write(*,*)  '[TEST : cubic_spline_t%derivative: PASSED] , rdiff : ', rdiff
     end if

    write(*,'(A, e20.13)')  '  * Compare vs exact result:  cubic_spline_t%integrate'
    ics = real(cs%integrate(0.5_r64,pi))
    rdiff = abs(ics - ana_integral)/abs(ana_integral)
    if ( rdiff > rtolerance) then
        write(*,*)  '[TEST : cubic_spline_t%integrate: FAILED] , rdiff : ', rdiff
        stop 1
     else 
        write(*,*)  '[TEST : cubic_spline_t%integrate: PASSED] , rdiff : ', rdiff
     end if 
    
end program test_spline_1d

