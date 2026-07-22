! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_putx0
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use precision, only: i32, long_int, dp

  implicit none

  contains

    !BOP
    ! !ROUTINE: putx0
    ! !INTERFACE:
    subroutine putx0(tp0, iq, iw, filnam, filxt, ch0, ch0wg, ch0hd)
    ! !USES:
      use modmpi
      use modmain
      use modxs
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! logical :: tp0 ! Flag if iq is the q=0 q-point
    ! integer :: iq  ! q-point index
    ! integer :: iw  ! iw'th freqency
    ! character(*) :: filnam, filxt ! File name and file extension
    ! complex(8) :: ch0(:,:)     ! Body of RPA density-density response matrix 
    ! complex(8) :: ch0wg(:,:,:) ! Wings of RPA density-density response matrix 
    ! complex(8) :: ch0hd(:,:)   ! Body of RPA density-density response matrix
    !
    ! !DESCRIPTION:
    !   Writes the V-symmetrized Kohn-Sham response function 
    !   $\tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)$ for one q-point and
    !   multiple frequencies into a direct access file. The records are
    !   numbered by the frequency index.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      logical(i32), intent(in) :: tp0
      integer(i32), intent(in) :: iq, iw
      character(*), intent(in) :: filnam, filxt
      complex(dp), intent(in) :: ch0(:, :)
      complex(dp), intent(in), optional :: ch0wg(:, :, :), ch0hd(:, :)

      ! Local variables
      character(*), parameter :: thisnam = 'putx0'
      integer(i32) :: un
      integer(long_int) :: reclen

      ! q=0 but head or wings missing
      if(tp0 .and. (.not. present(ch0wg) .or.  .not. present(ch0wg))) then
        write(*,*) 'Error(' // trim(thisnam) // '): q=0 but head or wings missing'
        call terminate
      end if

      if(tp0) then
        ! i/o record length
        call inquire_large( reclen, [ngq(iq)], vql(:, iq), ch0, ch0wg, ch0hd )
        call open_direct_unformatted_large( un, trim(filnam)//trim(filxt), "write", reclen, "unknown" )
        write(un, rec=iw) ngq(iq), vql(:, iq), ch0, ch0wg, ch0hd
      else
        ! i/o record length
        call inquire_large( reclen, [ngq(iq)], vql(:, iq), ch0 )
        call open_direct_unformatted_large( un, trim(filnam)//trim(filxt), "write", reclen, "unknown" )
        write(un, rec=iw) ngq(iq), vql(:, iq), ch0
      end if

      close(un)

    end subroutine putx0
    !EOC

end module m_putx0
