!> This module contains wrappers around 'inquire` and `open` statements with large 
!> buffer sizes. It should be removed as soon as large integers in those statements will 
!> be included in Fortran standard. The motivation of using wrappers instead of the 
!> in-place calls is hiding the default integer behaviour into a single file. 
!> The standard dictates inquire to use the default integer size. 
!> For compilers following this strictly (Cray) we need to add a change into 
!> the flags for specific files including inquire, this simply makes it easier 
!> for us we only need to do it for a single file.
module mod_large_io
#include "asserts.fpp"
  use mod_kpoint, only: nkpt, vkl
  use mod_names, only: filetag_occsv
  use modinput, only: input
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32, long_int, str_256

  implicit none
  private

  public :: inquire_large, open_direct_unformatted_large

  interface inquire_large
    module procedure :: &
      inquire_large_cdp, &
      inquire_large_rdp_rdp, &
      inquire_large_rdp_i32, &
      inquire_large_rdp_i32_cdpr2, &
      inquire_large_rdp_i32_cdp, &
      inquire_large_i32_rdp_i32_rdp_cdp_cdp_cdp, &
      inquire_large_i32_rdp_i32_rdp_cdp, &
      inquire_large_i32_cdp, &
      inquire_large_rdp_i32_rdp, &
      inquire_large_rdp_rdp_i32, &
      inquire_large_rdp_rdp_i32_cdp, &
      inquire_large_i32_cdp_cdp_cdp_cdp, &
      inquire_large_i32_long, &
      inquire_large_i32_short, &
      inquire_large_i32_i32_cdp, &
      inquire_large_cdpr1, &
      inquire_large_cdpr2, &
      inquire_large_i32r1, &
      inquire_large_i32_rdp_cdp, &
      inquire_large_i32_rdp_cdp_cdp_cdp, &
      inquire_large_i32_rdp_rdp_rdp_rdp, &
      inquire_large_cdp_cdp_cdp
  end interface

contains
  
  !> Wrapper of `open` statement with unformatted direct access. 
  subroutine open_direct_unformatted_large( io_unit, fname, action, large_iolength, file_status )
    !> IO unit
    integer(i32), intent(out) :: io_unit
    !> File name
    character(len=*), intent(in) :: fname
    !> IO action
    character(len=*), intent(in) :: action
    !> IO large recl
    integer(long_int), intent(in) :: large_iolength
    !> File status
    character(len=*), intent(in) :: file_status

    integer(i32) :: stat

    open (newunit=io_unit, file=fname, action = action, access="direct", &
      recl=large_iolength, form="unformatted", iostat = stat, status = file_status )
    call terminate_if_false( logical( stat == 0, kind = i32), "Error opening file: " // trim( fname ) )
  end subroutine

  subroutine inquire_large_cdpr1( large_iolength, array_1 )
    integer(long_int), intent(out) :: large_iolength
    complex(dp), intent(in) :: array_1(:)

    inquire( iolength = large_iolength ) array_1
  end subroutine

  subroutine inquire_large_cdpr2( large_iolength, array_1 )
    integer(long_int), intent(out) :: large_iolength
    complex(dp), intent(in) :: array_1(:, :)

    inquire( iolength = large_iolength ) array_1
  end subroutine

  subroutine inquire_large_i32r1( large_iolength, array_1 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)

    inquire( iolength = large_iolength ) array_1
  end subroutine

  subroutine inquire_large_cdp( large_iolength, array_1 )
    integer(long_int), intent(out) :: large_iolength
    complex(dp), intent(in) :: array_1(:, :, :)

    inquire( iolength = large_iolength ) array_1
  end subroutine

  subroutine inquire_large_rdp_rdp( large_iolength, array_1, array_2 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)

    inquire( iolength = large_iolength ) array_1, array_2  
  end subroutine

  subroutine inquire_large_rdp_i32( large_iolength, array_1, array_2 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)

    inquire( iolength = large_iolength ) array_1, array_2  
  end subroutine

  subroutine inquire_large_rdp_i32_cdpr2( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)
    complex(dp), intent(in) :: array_3(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

  subroutine inquire_large_i32_cdp( large_iolength, array_1, array_2 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    complex(dp), intent(in) :: array_2(:, :)

    inquire( iolength = large_iolength ) array_1, array_2  
  end subroutine

  subroutine inquire_large_i32_i32_cdp( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)
    complex(dp), intent(in) :: array_3(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

  subroutine inquire_large_i32_rdp_cdp( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    complex(dp), intent(in) :: array_3(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

  subroutine inquire_large_i32_rdp_cdp_cdp_cdp( large_iolength, array_1, array_2, &
      array_3, array_4, array_5 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    complex(dp), intent(in) :: array_3(:, :)
    complex(dp), intent(in) :: array_4(:, :, :)
    complex(dp), intent(in) :: array_5(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5
  end subroutine

  subroutine inquire_large_rdp_i32_cdp( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)
    complex(dp), intent(in) :: array_3(:, :, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

  subroutine inquire_large_rdp_i32_rdp( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)
    real(dp), intent(in) :: array_3(:)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

  subroutine inquire_large_rdp_rdp_i32( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    integer(i32), intent(in) :: array_3(:)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

  subroutine inquire_large_rdp_rdp_i32_cdp( large_iolength, array_1, array_2, array_3, array_4 )
    integer(long_int), intent(out) :: large_iolength
    real(dp), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    integer(i32), intent(in) :: array_3(:)
    complex(dp), intent(in) :: array_4(:, :, :, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4
  end subroutine

  subroutine inquire_large_i32_cdp_cdp_cdp_cdp( large_iolength, array_1, array_2, array_3, array_4, array_5 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    complex(dp), intent(in) :: array_2(:, :)
    complex(dp), intent(in) :: array_3(:, :)
    complex(dp), intent(in) :: array_4(:, :)
    complex(dp), intent(in) :: array_5(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5
  end subroutine

  subroutine inquire_large_i32_short( large_iolength, array_1, array_2, array_3, array_4, array_5, &
    array_6, array_7, array_8, array_9, array_10 )
    integer(long_int), intent(out) :: large_iolength
    logical(i32), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)
    integer(i32), intent(in) :: array_3(:)
    real(dp), intent(in) :: array_4(:)
    logical(i32), intent(in) :: array_5(:)
    real(dp), intent(in) :: array_6(:)
    real(dp), intent(in) :: array_7(:)
    integer(i32), intent(in) :: array_8(:)
    real(dp), intent(in) :: array_9(:)
    integer(i32), intent(in) :: array_10(:)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5, &
      array_6, array_7, array_8, array_9, array_10
  end subroutine

  subroutine inquire_large_i32_long( large_iolength, array_1, array_2, array_3, array_4, array_5, &
    array_6, array_7, array_8, array_9, array_10, array_11, array_12, array_13, array_14, array_15 )
    integer(long_int), intent(out) :: large_iolength
    logical(i32), intent(in) :: array_1(:)
    integer(i32), intent(in) :: array_2(:)
    integer(i32), intent(in) :: array_3(:)
    real(dp), intent(in) :: array_4(:)
    logical(i32), intent(in) :: array_5(:)
    real(dp), intent(in) :: array_6(:)
    real(dp), intent(in) :: array_7(:)
    integer(i32), intent(in) :: array_8(:)
    real(dp), intent(in) :: array_9(:)
    integer(i32), intent(in) :: array_10(:)
    integer(i32), intent(in) :: array_11(:)
    integer(i32), intent(in) :: array_12(:)
    integer(i32), intent(in) :: array_13(:, :)
    integer(i32), intent(in) :: array_14(:)
    integer(i32), intent(in) :: array_15(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5, &
      array_6, array_7, array_8, array_9, array_10, array_11, array_12, array_13, &
      array_14, array_15
  end subroutine

  subroutine inquire_large_i32_rdp_i32_rdp_cdp_cdp_cdp( large_iolength, array_1, array_2, array_3, array_4, array_5, array_6, array_7 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    integer(i32), intent(in) :: array_3(:)
    real(dp), intent(in) :: array_4(:)
    complex(dp), intent(in) :: array_5(:, :)
    complex(dp), intent(in) :: array_6(:, :, :)
    complex(dp), intent(in) :: array_7(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5, array_6, array_7
  end subroutine

  subroutine inquire_large_i32_rdp_i32_rdp_cdp( large_iolength, array_1, array_2, array_3, array_4, array_5 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    integer(i32), intent(in) :: array_3(:)
    real(dp), intent(in) :: array_4(:)
    complex(dp), intent(in) :: array_5(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5
  end subroutine

  subroutine inquire_large_i32_rdp_rdp_rdp_rdp( large_iolength, array_1, array_2, &
      array_3, array_4, array_5 )
    integer(long_int), intent(out) :: large_iolength
    integer(i32), intent(in) :: array_1(:)
    real(dp), intent(in) :: array_2(:)
    real(dp), intent(in) :: array_3(:)
    real(dp), intent(in) :: array_4(:)
    real(dp), intent(in) :: array_5(:)

    inquire( iolength = large_iolength ) array_1, array_2, array_3, array_4, array_5
  end subroutine

  subroutine inquire_large_cdp_cdp_cdp( large_iolength, array_1, array_2, array_3 )
    integer(long_int), intent(out) :: large_iolength
    complex(dp), intent(in) :: array_1(:, :)
    complex(dp), intent(in) :: array_2(:, :)
    complex(dp), intent(in) :: array_3(:, :)

    inquire( iolength = large_iolength ) array_1, array_2, array_3
  end subroutine

end module
