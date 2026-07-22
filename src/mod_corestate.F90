! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


!> Core state variables, if in modinput frozencore is set to true
!> core state wavefunctions, densities and
!> energies calculated only in the first iteration
module mod_corestate

      use precision, only: i32, dp

      implicit none

      private

      public :: save_corestate, load_corestate

      !> eigenvalues for core states
      real(dp), public, allocatable :: evalcr (:, :)
      !> radial wavefunctions for core states
      real(dp), public, allocatable :: rwfcr (:, :, :, :)
      !> radial charge density for core states
      real(dp), public, allocatable :: rhocr (:, :)

contains

      !> Saves core state wavefunctions, densities and
      !> energies
      subroutine save_corestate()

            integer(i32) :: unit

            open(newunit=unit, file="core.basis", form='unformatted', access='stream', status='replace')

            ! ===== evalcr =====
            write(unit) lbound(evalcr,1,i32), ubound(evalcr,1,i32)
            write(unit) lbound(evalcr,2,i32), ubound(evalcr,2,i32)
            write(unit) evalcr

            ! ===== rwfcr =====
            write(unit) lbound(rwfcr,1,i32), ubound(rwfcr,1,i32)
            write(unit) lbound(rwfcr,2,i32), ubound(rwfcr,2,i32)
            write(unit) lbound(rwfcr,3,i32), ubound(rwfcr,3,i32)
            write(unit) lbound(rwfcr,4,i32), ubound(rwfcr,4,i32)
            write(unit) rwfcr

            ! ===== rhocr =====
            write(unit) lbound(rhocr,1,i32), ubound(rhocr,1,i32)
            write(unit) lbound(rhocr,2,i32), ubound(rhocr,2,i32)
            write(unit) rhocr

            close(unit)
      end subroutine save_corestate

      !> Loads core state wavefunctions, densities and
      !> energies
      subroutine load_corestate()

            integer(i32) :: unit
            integer(i32) :: l1, u1, l2, u2
            integer(i32) :: l3, u3, l4, u4

            open(newunit=unit, file="core.basis", form='unformatted', access='stream', status='old')

            ! ===== evalcr =====
            read(unit) l1, u1
            read(unit) l2, u2
            if (allocated(evalcr)) deallocate(evalcr)
            allocate(evalcr(l1:u1, l2:u2))
            read(unit) evalcr

            ! ===== rwfcr =====
            read(unit) l1, u1
            read(unit) l2, u2
            read(unit) l3, u3
            read(unit) l4, u4
            if (allocated(rwfcr)) deallocate(rwfcr)
            allocate(rwfcr(l1:u1, l2:u2, l3:u3, l4:u4))
            read(unit) rwfcr

            ! ===== rhocr =====
            read(unit) l1, u1
            read(unit) l2, u2
            if (allocated(rhocr)) deallocate(rhocr)
            allocate(rhocr(l1:u1, l2:u2))
            read(unit) rhocr

            close(unit)
      end subroutine load_corestate

end module

