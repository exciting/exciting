! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine getscreen(iqr, ngq, scrh, scrw, scrb)
  use modinput, only: input
  use m_genfilname
  use m_getunit

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, ngq
  complex(8), intent(out) :: scrh(3, 3), scrw(ngq, 2, 3), scrb(ngq, ngq)

  ! Local variables
  character(256) :: fname
  real(8) :: rm(3, 3, 3)
  integer :: igq1, igq2, i, j, it1, it2, un, bzsampl

  ! Sampling of Brillouin zone
  bzsampl = 0
  if(input%xs%tetra%tetradf) bzsampl = 1

  ! Read in screening
  call genfilname(basename='SCREEN', iq=iqr, bzsampl=bzsampl, filnam=fname)
  call getunit(un)
  open(un, file=trim(fname), form='formatted', action='read', status='old')

  do igq1 = 1, ngq
    do igq2 = 1, ngq

      if(iqr .eq. 1) then

        ! Read head
        if((igq1 .eq. 1) .and. (igq2 .eq. 1)) then
          read(un,*) ((it1,it2,rm(1, i, j),rm(2, i, j),rm(3, i, j),j=1, 3) , i=1, 3)
          scrh(:, :) = cmplx(rm(1, :, :), rm(2, :, :), 8)
        end if

        ! Read wings
        if((igq1 .eq. 1) .and. (igq2 .ne. 1)) then
          read(un,*) (it1, it2, rm(1, 1, j), rm(2, 1, j), rm(3, 1, j), j=1, 3)
          scrw(igq2, 1, :) = cmplx(rm(1, 1, :), rm(2, 1, :), 8)
        end if
        if((igq1 .ne. 1) .and. (igq2 .eq. 1)) then
          read(un,*) (it1, it2, rm(1, 1, j), rm(2, 1, j), rm(3, 1, j), j=1, 3)
          scrw(igq1, 2, :) = cmplx(rm(1, 1, :), rm(2, 1, :), 8)
        end if

        ! Read body
        if((igq1 .ne. 1) .and. (igq2 .ne. 1)) then
          read(un,*) it1, it2, rm(1, 1, 1), rm(2, 1, 1), rm(3, 1, 1)
          scrb(igq1, igq2) = cmplx(rm(1, 1, 1), rm(2, 1, 1), 8)
        end if

      else
        
        ! Read full
        read(un,*) it1, it2, rm(1, 1, 1), rm(2, 1, 1), rm(3, 1, 1)
        scrb(igq1, igq2) = cmplx(rm(1, 1, 1), rm(2, 1, 1), 8)

      end if

    end do
  end do

  close(un)

end subroutine getscreen
