! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 

!> Data precision module 
module precision
    use, intrinsic :: iso_fortran_env, Only: int16, int32, int64, real32, real64, real128
    implicit none
    private

    !> Single precision
    integer, parameter, public :: sp = real32
    !> Double precision
    integer, parameter, public :: dp = real64
    !> Quadruple precision
    integer, parameter, public :: qp = real128

    !> Working precision
    !> Gives the flexibility to recompile complete code blocks
    !> with different precision, if memory is an issue
    integer, parameter, public :: wp = dp

    !> Normal 4 byte integer
    integer, parameter, public :: i32 = int32
    !> Long integer
    integer, parameter, Public :: long_int = int64
    !> Short integer
    integer, parameter, Public :: short_int = int16

    !> 16 character string
    integer, parameter, public :: str_16 = 16
    !> 32 character string
    integer, parameter, public :: str_32 = 32
    !> 64 character string
    integer, parameter, public :: str_64 = 64
    !> 128 character string
    integer, parameter, public :: str_128 = 128
    !> 256 character string
    integer, parameter, public :: str_256 = 256
    !> 512 character string
    integer, parameter, public :: str_512 = 512
    !> length of a long string
    integer, parameter, public :: str_1024 = 1024

end module
