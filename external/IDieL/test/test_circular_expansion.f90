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
!> Tests the circular expansion and quadrature
program test_circular_expansion

    use idiel_constants, only: aip, i32, twopi, pi, fourpi, zzero
    use idiel_circular_harmonics
    use idiel_circle_quadrature 

    implicit none

    integer(i32), parameter :: msize = 75_i32
    integer(i32), parameter :: lmax  = 20_i32
#ifdef USE_SINGLE_PRECISION
    real(aip), parameter    :: tolerance = 1.0e-4_aip
#else
    real(aip), parameter    :: tolerance = 1.0e-12_aip
#endif
    real(aip), allocatable :: rphi(:,:), w(:), xyz(:,:)
    complex(aip), allocatable :: blm(:,:)
    complex(aip), allocatable :: f2expand(:), fexact(:), fint(:)
    complex(aip) :: clm(2*lmax+1)
    real(aip) :: integral, rdiff
    integer(i32) :: i


    write(*,*) '[TEST : idiel_circular_harmonics and idiel_circle_quadrature]' 
    write(*,*)  '  * Check exact integration of circular basis'
    ! We compute the angular mesh
    call compute_angular_mesh_gauss_legendre(msize, rphi, w, xyz)

    ! We compute the circular harmonics up to 22
    call circ_harm(msize, rphi(:,2), blm)

    ! Check that the quadrature integrates the circular harmonics
    integral = real(sum(w * blm(:,1)))
    if ( abs(integral-2.50662827463100_aip)/2.50662827463100_aip <= tolerance) then
        write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: PASSED]', rdiff
    else   
        write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: FAILED]', rdiff
        stop 1
    end if

    do i = 2, 2*lmax+1
        integral = real(sum(w * blm(:,i)))
        if ( integral <= tolerance) then
            write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: PASSED]', rdiff
        else   
            write(*,*)  '[TEST : idiel_circular_harmonics and idiel_circle_quadrature: FAILED]', rdiff
            stop 1
        end if
    end do

    write(*,*)  '[TEST : circ_harm_expansion]'
    f2expand = f1(xyz)
    call circ_harm_expansion(lmax, f2expand, w, blm, clm)
    deallocate(rphi, w, xyz, blm)
    
    call compute_angular_mesh_gauss_legendre(10_i32*msize, rphi, w, xyz)
    call circ_harm(msize, rphi(:,2), blm)
    fexact = f1(xyz)
    allocate(fint, mold=fexact)
    fint = zzero
    
    ! Compute the expansion
    do i = 1, 2*lmax+1
        fint(:) = fint(:) + clm(i) * blm(:,i) 
    end do

    ! Compute the function
    do i = 1, 10*msize
        rdiff = abs(real(fint(i))-real(fexact(i))) / abs(real(fexact(i)))
        if ( rdiff <= tolerance) then
            write(*,*)  '[TEST : circ_harm_expansion: PASSED]', rdiff
        else   
            write(*,*)  '[TEST : circ_harm_expansion: FAILED]', rdiff
            stop 1
        end if
    end do 

    stop 0

contains

    pure function f1(r)
        !> Mesh in which to compute
        real(aip), intent(in)  :: r(:,:)
        !> Result
        complex(aip), allocatable :: f1(:)
        !> Elements
        real(aip), allocatable :: x(:), y(:)

        x = r(:,1)
        y = r(:,2)
        
        allocate(f1(size(r,1)))

        f1 = cmplx(1.0_aip + x + y**2 + y * x**2 + x**4 + y**5 + (x*y)**3, 0.0_aip, aip)

    end function f1

end program test_circular_expansion
