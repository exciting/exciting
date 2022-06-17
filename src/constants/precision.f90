! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 

!> Data precision module 
module precision
    implicit none    
    !> Single precision 
    integer, public,parameter :: sp = selected_real_kind(6, 37)
    !> Double precision 
    integer, public,parameter :: dp = selected_real_kind(15, 307)
    !> Quadruple precision 
    integer, public,parameter :: qp = selected_real_kind(33, 4931)
    !> 16 character string
    integer, parameter :: str_16 = 16
    !> 16 character string
    integer, parameter :: str_32 = 32
    !> 64 character string
    integer, parameter :: str_64 = 64
    !> 128 character string
    integer, parameter :: str_128 = 128
    !> 256 character string
    integer, parameter :: str_256 = 256
    !> 512 character string
    integer, parameter :: str_512 = 512
    !> length of a long string
    integer, parameter :: str_1024 = 1024
    !> format character for converting a **real(sp)** variable to a string
    character(*), parameter :: sp_to_char = '(ES13.6)'
    !> format character for converting a **real(dp)** variable to a string
    character(*), parameter :: dp_to_char = '(ES21.14)'
end module