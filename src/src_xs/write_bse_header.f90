! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> module for writing a common header to BSE out files
module m_write_bse_header
  use modinput, only: input_type
  use modxs, only: escale, ivgmt, vqlmt, vgcmt, vqcmt, sptclg, ivgigq
  use mod_lattice, only: omega
  use modbse, only: nk_bse
  use to_char_conversion, only: to_char
  use string_utils, only: append_line

  implicit none

  contains

    !> routine that creates a common header to insert in a BSE out file
    function generate_bse_header(input, iq, iop1, iop2) result(header_string)

      ! Arguments
      !> input file container
      type(input_type), intent(in) :: input
      !> indices for QMT point and oscillator direction
      !> oscillator indices are ignored if QMT > 1
      integer, intent(in) :: iq, iop1, iop2

      !> result
      character(:), allocatable :: header_string

      ! Local variables
      !> Gmt+qmt index
      integer :: igqmt
      !> temp string
      character(200) :: tmp

      ! Note: If you change the number of lines here, adjust also in `readoscillator` ncommentlines

      ! Get Gmt+qmt index
      igqmt = ivgigq(ivgmt(1,iq),ivgmt(2,iq),ivgmt(3,iq),iq)

      header_string = '#'
      call append_line(header_string, '# Momentum transfer Q=G+q in lattice coordinates')
      write(tmp, '("# G:",3i4)') ivgmt(1:3,iq)
      call append_line(header_string, tmp)
      write(tmp, '("# q:",3f12.7)') vqlmt(1:3,iq)
      call append_line(header_string, tmp)
      call append_line(header_string, "# Momentum transfer Q=G+q in Cartesian coordinates")
      write(tmp, '("# G:",3f12.7)') vgcmt(1:3,iq)
      call append_line(header_string, tmp)
      write(tmp, '("# q:",3f12.7)') vqcmt(1:3,iq)
      call append_line(header_string, tmp)
      write(tmp, '("# Norm2(G+q)",f12.7)') norm2(vgcmt(:,iq)+vqcmt(:,iq))
      call append_line(header_string, tmp)
      call append_line(header_string, "#")
      write(tmp, '("# Energy scale=", f12.7)') escale
      call append_line(header_string, tmp)
      write(tmp, '("# Used broadening in scaled energy units:", f12.6)')&
        & escale*input%xs%broad
      call append_line(header_string, tmp)
      call append_line(header_string, "#")
      write(tmp, '("# Number of k-points=", i8)') nk_bse
      call append_line(header_string, tmp)
      write(tmp, '("# Unit cell volume [au]=", f12.7)') omega
      call append_line(header_string, tmp)
      write(tmp, '("# Coulomb potential v(Q) [au]=", f18.7)') sptclg(igqmt,iq)**2
      call append_line(header_string, tmp)
      call append_line(header_string, "#")
      ! Note that both branches must have the same length in the header: for `read_oscillator` number of commentlines
      ! is hard-coded
      if (iq == 1) then
        write(tmp, '("# First optical index=", i4)') iop1
        call append_line(header_string, tmp)
        write(tmp, '("# Second optical index=", i4)') iop2
        call append_line(header_string, tmp)
      else
        write(tmp, '("# QMT index=", i4)') iq
        call append_line(header_string, tmp)
        write(tmp, '("# QMT vector= ", a)') to_char(input%xs%qpointset%qpoint(:, iq))
        call append_line(header_string, tmp)
      end if
      write(tmp, '("# BSE type=", a)') trim(input%xs%bse%bsetype)
      call append_line(header_string, tmp)
      write(tmp, '("# TDA used: ", a)') to_char(.not. input%xs%bse%coupling)
      call append_line(header_string, tmp)
      write(tmp, '("# Screening type=", a)') trim(input%xs%screening%screentype)
      call append_line(header_string, tmp)
      write(tmp, '("# Antiresonant part used: ", a)') to_char(input%xs%bse%aresbse)
      call append_line(header_string, tmp)
      write(tmp, '("# Chi bar finite q: ", a)') to_char(input%xs%bse%chibarq)
      call append_line(header_string, tmp)
      write(tmp, '("# Chi bar q=0: ", a)') to_char(input%xs%bse%chibar0)
      call append_line(header_string, tmp)
      call append_line(header_string, "#")

   end function generate_bse_header

end module m_write_bse_header
