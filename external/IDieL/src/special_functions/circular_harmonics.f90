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
!> Contains elements to compute circular harmonics and related elements

!> This module contains procedures to compute circular harmonics and functions expansion in terms of such a basis
module idiel_circular_harmonics
   
   use idiel_constants,   only: i32, aip, pi, twopi, zone, zzero

   implicit none 

   private
   public :: circ_harm, circ_harm_expansion


contains 
    
   !>   Generates a sequence of circular harmonics
   !>   @param[in] lmax - the maximum momentum for which the circular harmonic are computed
   !>   @param[in] ang  - the angular points for which the circular harmonics are computed 
   !>   @param[out] blm - the spherical harmonics at ang points
   subroutine circ_harm(lmax, ang, blm)
   
      integer(i32), intent(in)  :: lmax
      real(aip),    intent(in)  :: ang(:)
      complex(aip), allocatable, intent(out) :: blm(:,:)
      
      ! local variables
      integer(i32) :: l, m, npoints
      
      npoints = size(ang)
      
      allocate(blm(npoints,(2*lmax+1)), source = cmplx(1.0/sqrt(pi), 0.0, aip))
      
      blm(:,1) = cmplx(1.0/sqrt(twopi), 0.0, aip)

      !$omp parallel default(shared) private(l, m)
      !$omp do schedule(dynamic)
      do l = 1, lmax
         blm(:, 2*l    ) = blm(:, 2*l    ) * cmplx(sin(l * ang(:)), 0.0, aip)
         blm(:, 2*l + 1) = blm(:, 2*l + 1) * cmplx(cos(l * ang(:)), 0.0, aip)
      end do
      !$omp end do
      !$omp end parallel 
   
   end subroutine circ_harm
    
   !> It expands the function f in terms of circular harmonics (up to l = 20)
   !> Notice that this function is intended to be used with a Gauss-Legendre grid of size 20 which integrates all harmonics up to 10
   !> expansion coefficients should be accurate up to l = 20
   !> @param[in] lexp     - the Fourier expansion order
   !> @param[in] f        - the function to expand (the values for different radius are provided in the radii)
   !> @param[in] weights  - the mesh points weights
   !> @param[in] blm      - the spherical harmonics in the mesh
   !> @param[out] clm     - the expansion coefficients
   subroutine circ_harm_expansion(lexp, f, weights, blm, clm)
#if defined(USE_GPU) && defined(HAVEOMP5)
      !$omp declare target
#endif  
      integer(i32), intent(in)   :: lexp
      complex(aip), intent(in)   :: f(:)
      real(aip),    intent(in)   :: weights(:)
      complex(aip), intent(in)   :: blm(:,:)
      complex(aip), intent(out)  :: clm(2*lexp+1)


      complex(aip), allocatable :: fw(:)
      integer(i32) :: i
      ! Multiply the function with the weights so the matrix-vector products solves the integral for all the (l,m)
      allocate(fw, source=f)
      fw(:) = fw(:) * weights(:)
      ! This is a small calculation so it is done in the CPU (Notice that we use only T not H since blm imaginary part is zero)
      ! clm = blm**T \cdot f
      ! We need to write it by hand (GPU not supporting matmul intrinscic)
      do i = 1, size(clm)
        clm(i) = sum(blm(:,i) * fw(:))
      end do
   end subroutine circ_harm_expansion

end module idiel_circular_harmonics
