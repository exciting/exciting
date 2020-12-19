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
end module