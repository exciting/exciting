! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_genwgrid
  
  implicit none

  contains

    !BOP
    ! !ROUTINE: genwgrid
    ! !INTERFACE:
    subroutine genwgrid(n, intv, timag, brd, w_real, w_cmplx)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer :: n       ! Number of grid points
    ! real(8) :: intv(2) ! Intervall limits 
    ! logical, optional :: timag ! Flag for imaginary fequencies 
    ! real(8), optional :: brd   ! Broadening parameter
    ! Out:
    ! real(8), optional :: w_real(:)     ! Real valued grid 
    ! complex(8), optional :: w_cmplx(:) ! Complex valued grid 
    !
    ! !DESCRIPTION:
    !   Generates an evenly spaced "frequency" grid.
    !   The gird is either reals OR complex valued.
    !   If timag=true then all frequencies are multiplied by i, otherwise
    !   the complex grid is mathematically identical to the real valued one.
    !   If brd is non zero the frequencies of the complex grid are shifted
    !   by i*brd.
    ! 
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      integer, intent(in) :: n
      real(8), intent(in) :: intv(2)

      ! Optional arguments
      logical, optional, intent(in) :: timag
      real(8), optional, intent(in) :: brd
      real(8), optional, intent(out) :: w_real(:)
      complex(8), optional, intent(out) :: w_cmplx(:)

      ! Local variables
      character(*), parameter :: thisnam = 'genwgrid'
      integer :: j, nerr
      real(8) :: brdt, t1
      complex(8) :: fac
      
      ! checking
      nerr = 0
      if(present(w_real) .and. present(w_cmplx)) then
        write(*, '("Error(", a, "both real and complex output grid specified")') thisnam
        nerr = nerr + 1
      end if

      if(nerr .gt. 0) stop
      
      if(present(w_real)) then

        ! Real grid
        t1 = (intv(2)-intv(1)) / dble(n)
        do j = 1, n
          w_real(j) = t1 * dble(j-1) + intv(1)
        end do

      else
        ! Complex grid

        brdt = 0.d0
        fac = (1.d0, 0.d0)
        if(present(brd)) brdt = brd

        if(present(timag)) then
          if(timag) fac = (0.d0, 1.d0)
        end if

        t1 = (intv(2)-intv(1)) / dble(n)

        do j = 1, n
          w_cmplx(j) = fac * (t1*dble(j-1)+intv(1)) + (0.d0, 1.d0) * brdt
        end do

      end if
      
    end subroutine genwgrid
    !EOC

end module m_genwgrid
