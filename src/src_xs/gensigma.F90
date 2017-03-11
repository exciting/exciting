! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gensigma

  implicit none

  contains

    subroutine gensigma(w, eps, oc, sigma)
      use mod_constants, only: pi, zi
      use modmpi
      use modxs
      implicit none

      ! Arguments
      real(8), intent(in) :: w(:)
      complex(8), intent(in) :: eps(:)
      integer, intent(in) :: oc(2)
      complex(8), intent(out) :: sigma(:)

      ! Local variables
      character(*), parameter :: thisnam = 'gensigma'
      real(8) :: delt

      if(any(shape(eps) .ne. shape(sigma))) then
        write(unitout, '(a)') 'Error(' // thisnam // '): input and&
          & output arrays have diffenrent shape'
        call terminate
      end if

      ! Optical conductivity
      delt = 0.d0
      if(oc(1) .eq. oc(2)) delt = 1.d0

      sigma(:) = - zi * (eps(:) - delt) * w(:) / (4.d0*pi)

    end subroutine gensigma

end module m_gensigma
