! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeeps

  implicit none

  contains

    subroutine writeeps(iq, iop1, iop2, w, eps, fn)
      use modinput
      use modmpi
      use modxs, only: escale, unitout, unit1
      use m_write_bse_header, only: generate_bse_header

      implicit none

      ! Arguments
      integer, intent(in) :: iq, iop1, iop2
      real(8), intent(in) :: w(:)
      complex(8), intent(in) :: eps(:)
      character(*), intent(in) :: fn

      ! Local variables
      character(*), parameter :: thisnam = 'writeeps'
      integer :: n, iw
      real(8), allocatable :: imeps(:), kkeps(:)
      if(any(shape(w) .ne. shape(eps))) then
        write(unitout, '(a)') 'Error(' // thisnam // '): input&
          & arrays have different shapes'
        call terminate
      end if

      n = size(w)

      allocate(imeps(n), kkeps(n))

      ! Kramers-Kronig transform imaginary part
      imeps(:) = aimag(eps(:))

      Call kramkron(iop1, iop2, 1.d-8, n, w, imeps, kkeps)

      Open(newunit=unit1, File=trim(fn), Action='write')
      write(unit1, '("# Macroscopic dielectric function epsm")')
      write(unit1, '(a)') generate_bse_header(input, iq, iop1, iop2)
      write(unit1, '("#",a22,1x,a23,1x,a23,1x,a23)')&
        & "Frequency", "Re(epsm)", "Im(epsm)", "Re(epsm) form KKT"
      write(unit1, '(SP,E23.16,1x,E23.16,1x,E23.16,1x,E23.16)')&
        & (w(iw)*escale, eps(iw), kkeps(iw), iw=1, n)
      Close(unit1)

      deallocate(imeps, kkeps)
   end subroutine writeeps

end module m_writeeps
