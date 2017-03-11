! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_genloss
  use modmpi
  use modxs
  use modinput

  implicit none

  contains

    subroutine genloss(eps, loss, nc)
      implicit none

      ! Arguments
      complex(8), intent(in) :: eps(:, :, :)
      real(8), intent(out) :: loss(:, :, :)
      integer, intent(in) :: nc

      ! Local variables
      integer :: iw
      complex(8) :: t3(3, 3)
      character(*), parameter :: thisnam = 'genloss'

      if(nc /= 1 .and. nc /= 3) then 
        write(unitout, '(a,i2)') 'Error(' // thisnam // '): nc invalid, nc=', nc
        call terminate
      end if

      if(any(shape(eps) .ne. shape(loss))) then
        write(unitout, '(a)') 'Error(' // thisnam // '): input and&
          & output arrays have diffenrent shape'
        call terminate
      end if

      ! Loss function
      ! stk
      if(nc .eq. 1) then

        loss(1,1,:) = - aimag(1.0d0/eps(1,1,:))

      else if(nc .eq. 3) then

        do iw = 1, input%xs%energywindow%points
          call z3minv(eps(:, :, iw), t3(:, :))
          loss(:, :, iw) = -aimag(t3)
        end do

      end if
    end subroutine genloss

end module m_genloss
