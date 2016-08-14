! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine putscreen(un, tq0, n, chi0, chi0h, chi0w)
  use mod_constants, only: krondelta

  implicit none

  ! Input parameters
  logical, intent(in) :: tq0
  integer, intent(in) :: un, n
  complex(8), intent(in) :: chi0(n, n), chi0h(3, 3), chi0w(n, 2, 3)

  ! Local variables
  integer :: ig1, ig2, i, j
  real(8) :: r1

  ! Loop over G+q and G'+q indices
  do ig1 = 1, n
    do ig2 = 1, n

      r1 = 0.d0
      if(ig1 .eq. ig2) r1 = 1.d0

      if(tq0) then

        ! Write head
        if((ig1 .eq. 1) .and. (ig2 .eq. 1)) then
          write(un, '(2i8,3g18.10)') ((-i,-j,&
            & dble(krondelta(i, j))-chi0h(i, j),&
            & abs(dble(krondelta(i, j))-chi0h(i, j)), j=1, 3),&
            & i=1, 3)
        end if

        ! Write wings
        if((ig1 .eq. 1) .and. (ig2 .ne. 1)) then
          write(un, '(2i8,3g18.10)') (-i, ig2,-chi0w(ig2, 1,&
            & i), abs(-chi0w(ig2, 1, i)), i=1, 3)
        end if
        if((ig1 .ne. 1) .and. (ig2 .eq. 1)) then
          write(un, '(2i8,3g18.10)') (ig1,-j,-chi0w(ig1, 2,&
            & j), abs(-chi0w(ig1, 2, j)), j=1, 3)
        end if

        ! Write body
        if((ig1 .ne. 1) .and. (ig2 .ne. 1)) then
          write(un, '(2i8,3g18.10)') ig1, ig2, r1 - chi0(ig1,&
            & ig2), abs(r1-chi0(ig1, ig2))
        end if

      else

        ! Write full
        write(un, '(2i8,3g18.10)') ig1, ig2, r1 - chi0(ig1,&
          & ig2), abs(r1-chi0(ig1, ig2))

      end if

    end do
  end do

end subroutine
