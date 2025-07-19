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
!> This module contains parameters and constans
!>
module idiel_constants
    
    use iso_fortran_env, only: real32, real64, int32
    
    implicit none
    
    public 
    
    !> Integer single precision
    integer, parameter :: i32  = int32
    !> Real single precision
    integer, parameter :: r32  = real32
    !> Real double precision
    integer, parameter :: r64  = real64

#ifdef USE_SINGLE_PRECISION
    !> Real precision used by the ab initio code
    integer, parameter :: aip  = real32
#else
    !> Real precision used by the ab initio code
    integer, parameter :: aip  = real64
#endif

    !> Pi 
    real(r64), parameter :: pi       = 4.0_r64 * atan(1.0_r64)
    !> Two pi
    real(r64), parameter :: twopi    = 2.0_r64 * pi
    !> Four pi
    real(r64), parameter :: fourpi   = 4.0_r64 * pi
    !> 1.0/3.0
    real(r64), parameter :: onethird = 1.0_r64 / 3.0_r64
    !> The Y_{0}^{0} spherical harmonic
    complex(r64), parameter :: y00   = (0.28209479177387814_r64,  0.0_r64)
    
    !> Imaginary unit
    complex(aip), parameter :: iunit = cmplx(0.0_aip, 1.0_aip, kind=aip)
    !> Complex one
    complex(aip), parameter :: zone  = cmplx(1.0_aip, 0.0_aip, kind=aip)
    !> Complex zero
    complex(aip), parameter :: zzero = cmplx(0.0_aip, 0.0_aip, kind=aip)
    
end module idiel_constants
