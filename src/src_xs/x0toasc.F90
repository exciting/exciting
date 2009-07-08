

! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine x0toasc
  use modmain
use modinput
  use modxs
  use m_getx0
  use m_putx0
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='x0toasc'
  character(256) :: filnam, filnama
  integer :: n, iq, igq, igqp, iw, oct1, oct2, noct, un
  complex(8) :: zt
  complex(8), allocatable :: chi0(:, :), chi0wg(:, :, :), chi0hd(:, :)
  logical :: tq0
  logical, external :: tqgamma
  call init0
  ! initialise universal variables
  call init1
  ! save Gamma-point variables
  call xssave0
  ! initialize q-point set
  call init2
  ! loop over q-points
  do iq=1, nqpt
     tq0=tqgamma(iq)
     ! calculate k+q and G+k+q related variables
     call init1offs(qvkloff(1, iq))
     ! size of local field effects
     n=ngq(iq)
     ! allocate
     allocate(chi0(n, n), chi0wg(n, 2, 3), chi0hd(3, 3))
     ! filenames
     call genfilname(asc = .true., basename = 'X0', bzsampl = bzsampl, acont = input%xs%tddft%acont, &
	  nar = .not.input%xs%tddft%aresdf, iqmt = iq, filnam = filnama)
     call genfilname(basename = 'X0', bzsampl = bzsampl, acont = input%xs%tddft%acont, &
	  nar = .not.input%xs%tddft%aresdf, iqmt = iq, filnam = filnam)
     ! open file to write ASCI
     call getunit(un)
     open(unit = un, file = trim(filnama), form = 'formatted', &
	  action = 'write', status = 'replace')
     noct=1
     if (tq0) noct=3
     do iw=1, nwdf
        ! read from binary file
	call getx0(tq0, iq, iw, trim(filnam), '', chi0, chi0wg, chi0hd)
        ! write to ASCII file
	do igq=1, ngq(iq)
	   do igqp=1, ngq(iq)
	      if (tq0) then
		 do oct1=1, noct
		    do oct2=1, noct
                       ! head
		       if ((igq.eq.1).and.(igqp.eq.1)) &
			    chi0(igq, igqp) = chi0hd(oct1, oct2)
                       ! wings
		       if ((n.gt.1).and.(igq.eq.1).and.(igqp.gt.1)) &
			    chi0(igq, igqp) = chi0wg(igqp, 1, oct1)
		       if ((n.gt.1).and.(igq.gt.1).and.(igqp.eq.1)) &
			    chi0(igq, igqp) = chi0wg(igq, 2, oct2)		
		       zt=chi0(igq, igqp)
		       write(un, '(6i6, 3g18.10)') iq, iw, igq, igqp, oct1, oct2, zt, &
			    abs(zt)
		    end do
		 end do
	      else
		 zt=chi0(igq, igqp)
		 write(un, '(6i6, 3g18.10)') iq, iw, igq, igqp, 0, 0, zt, abs(zt)
	      end if
	   end do
	end do
     end do
     ! close file
     close(un)
     deallocate(chi0, chi0wg, chi0hd)
     write(unitout, '(a, i8)') 'Info('//thisnam//'): Kohn Sham response &
	  &function converted to ASCII file for q - point:', iq
  end do
  call genfilname(setfilext=.true.)
end subroutine x0toasc
