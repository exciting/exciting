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
!> Contains elements to compute Gauss-Legendre grid for a unitary
!> circle
module idiel_circle_quadrature
        
    use idiel_constants,   only: i32, r64, aip, pi, twopi
    use iso_c_binding,     only: C_int, C_ptr, C_loc

    implicit none
    
    private 
    public :: compute_angular_mesh_gauss_legendre

contains 

    !> Polar to cartesian coordinates
    !> @param[in] rphi - polar mesh (r,phi)
    !> @result    xyz - the Cartesian coordinates 
    pure function polar_2_cartesian(rphi) result(xyz)
        !> The angular mesh (theta,phi)
        real(aip), intent(in) :: rphi(:,:)
        !> The Cartesian mesh
        real(aip), allocatable :: xyz(:,:)

        allocate(xyz(size(rphi,1),3))

        xyz(:,1) = rphi(:,1) * cos(rphi(:,2))
        xyz(:,2) = rphi(:,1) * sin(rphi(:,2))
        xyz(:,3) = 0.0_aip

    end function polar_2_cartesian

    !> Cartesian to polar coordinates
    !> @param[in] xyz - the Cartesian coordinates
    !> @result    rphi - polar mesh (r,phi) 
    pure function cartesian_2_polar(xyz) result(rphi)
        !> The Cartesian mesh
        real(aip), intent(in)  :: xyz(:,:)
        !> The angular mesh (theta,phi)
        real(aip), allocatable :: rphi(:,:)

        allocate(rphi(size(xyz,1),2), source=0.0_aip)

        rphi(:,1) = hypot(xyz(:,1),xyz(:,2))
        rphi(:,2) = atan2(xyz(:,2),xyz(:,1))

    end function cartesian_2_polar

    !> Gauss-Legendre quadrature using Golub–Welsch algorithm
    !> matches in precision accurate implementation of Boost C++. 
    !> @param[in]   n - the number of points
    !> @param[out]  x - the abscisa points
    !> @param[out]  w - weights
    subroutine gauss_legendre_golub_welsch(n, x, w)

        integer(i32), intent(in)  :: n
        real(aip),    intent(out) :: x(:)
        real(aip),    intent(out) :: w(:)

        real(r64), allocatable :: jacobi(:,:)
        real(r64), allocatable :: work(:)
        real(r64), allocatable :: eigenvalues(:)
        integer(i32)  :: i, info

        external :: dsyev

        ! Initialize the Jacobi matrix with beta coefficients
        allocate(jacobi(n,n), source = 0.0_r64)
        do i = 1, n - 1
            jacobi(i,i+1) = i / sqrt(4.0_r64 * i * i - 1.0_r64)
            jacobi(i+1,i) = jacobi(i,i+1)
        end do

        ! Compute eigenvalues (knots) and eigenvectors
        allocate(eigenvalues(n), work(3*n-1))

        call dsyev('v', 'u', n, jacobi, n, eigenvalues, work, 3*n-1, info)
        if (info /= 0) error stop "Error(gauss_legendre_golub_welsch) dsyev failed"

        x(:) = real(eigenvalues(:), kind=aip)
        w(:) = real(2_r64*jacobi(1,:)**2, kind=aip)

        deallocate(eigenvalues, jacobi, work)

    end subroutine gauss_legendre_golub_welsch

    !> Gauss-Legendre quadrature using G. Rybicki approach
    !> @param[in]   n - the number of points
    !> @param[out]  x - the abscisa points
    !> @param[out]  w - weights
    subroutine gauss_legendre_rybicki(n, x, w)

        integer(i32), intent(in)  :: n
        real(aip),    intent(out) :: x(:)
        real(aip),    intent(out) :: w(:)

        integer(i32) :: i, j
        real(r64), parameter :: tol = 1.0e-13_r64
        real(r64) :: root, root_old
        real(r64) :: p1, p2, p3, pp

        ! We exploit the symmetry of the quadrature
        do i=1, (n+1)/2

            ! Use Newton's method
            root     = cos(pi*(i-0.25_r64)/(n+0.5_r64))
            root_old = huge(root)

            do while(abs(root-root_old) >= tol)

                p1=1.0_r64
                p2=0.0_r64

                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0_r64*j - 1.0_r64) * root * p2 - (j-1.0_r64) * p3) / j
                end do
                ! compute the derivative

                pp = n * (p2 - root * p1)/(1.0_r64-root**2)
                root_old = root
                root  = root_old - p1 / pp

            end do

            x(i)     = real(-root, kind = aip)
            x(n+1-i) = real( root, kind = aip)
            w(i)     = real(2.0_r64 / ((1.0_r64 - root**2) * pp**2), kind = aip)
            w(n+1-i) = w(i)

        end do

    end subroutine gauss_legendre_rybicki

    !> Compute the gauss legendre mesh for arbitrary range
    !> @param[in]   n - the number of points 
    !> @param[in]   a - the lower bound of the integral
    !> @param[in]   b - the upper bound of the integral
    !> @param[out]  x - the abscisa points
    !> @param[out]  w - weights
    subroutine compute_gauss_legendre(n, a, b, x, w)

        integer(i32), intent(in) :: n
        real(aip), intent(in) :: a
        real(aip), intent(in) :: b
        real(aip), allocatable, target, intent(out) :: x(:)
        real(aip), allocatable, target, intent(out) :: w(:)

        real(aip) :: dx, shift
        integer   :: i
         
        allocate(x(n),w(n))

        call gauss_legendre_golub_welsch(n, x, w)

        ! Remap
        dx    = 0.5_aip * ( b - a )
        shift = 0.5_aip * ( a + b )
        x(:)  = dx * x(:) + shift
        w(:)  = dx * w(:)

    end subroutine compute_gauss_legendre

    !> Compute an angular mesh Gauss-Legendre for the angle
    !> @param[in]   mesh_size - the mesh size
    !> @param[out]  rphi      - the angles (r,phi)
    !> @param[out]  w         - weights
    !> @param[out]  xyz       - the mesh in cartesian coordinates
    subroutine compute_angular_mesh_gauss_legendre(mesh_size, rphi, w, xyz)
        
        integer(i32), intent(in)            :: mesh_size
        real(aip), allocatable, intent(out) :: rphi(:,:)
        real(aip), allocatable, intent(out) :: w(:)
        real(aip), allocatable, intent(out) :: xyz(:,:)

        ! Phi mesh
        real(aip), allocatable :: x_phi(:)

        ! Idx
        integer(i32) :: i, idx 

        ! Build theta mesh
        call compute_gauss_legendre(mesh_size, 0.0_aip, real(twopi, aip), x_phi, w)

        allocate(rphi(mesh_size,2), source=1.0_aip)

        ! Save in proper format
        rphi(:,2) = x_phi(:)

        ! Go to cartesian
        xyz =  polar_2_cartesian(rphi)

    end subroutine compute_angular_mesh_gauss_legendre

end module idiel_circle_quadrature
