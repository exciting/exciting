! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_genloss

  implicit none

  contains

    subroutine genloss(eps, loss, nc)
      use modxs
      use modinput
      implicit none

      ! Arguments
      complex(8), intent(in) :: eps(:, :, :)
      real(8), intent(out) :: loss(:, :, :)
      integer, intent(in) :: nc

      ! Local variables
      integer :: iw
      complex(8) :: t3(3, 3)
      character(*), parameter :: thisnam = 'genloss'


      if(any(shape(eps) .ne. shape(loss))) then
        write(unitout, '(a)') 'Error(' // thisnam // '): input and&
          & output arrays have diffenrent shape'
        call terminate
      end if

      ! Loss function
      ! stk
      if(nc.eq.1) then

        loss(1,1,:) = - aimag(1/eps(1,1,:))

      else

        do iw = 1, input%xs%energywindow%points
          call z3minv(eps(:, :, iw), t3(:, :))
          loss(:, :, iw) = -aimag(t3)
        enddo

      end if
    end subroutine genloss

end module m_genloss
