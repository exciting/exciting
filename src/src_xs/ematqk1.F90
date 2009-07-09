


! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine ematqk1(iq, ik)
  use modmain
use modinput
  use modxs
  use modmpi
  use m_putemat
  implicit none
  ! arguments
  integer, intent(in) :: iq, ik
  ! set band combinations
  if (.not.(task.eq.430)) then
     call ematbdlims(2*input%xs%emattype, nst1, istl1, istu1, nst2, istl2, istu2)
     if (allocated(xiou)) deallocate(xiou)
     if (allocated(xiuo)) deallocate(xiuo)
     allocate(xiou(nst1, nst2, ngq(iq)))
     call ematqk(iq, ik)
  end if
  if (input%xs%emattype.eq.0) then
     ! all band combinations
     nst3=nstsv; nst4=nstsv
     if (.not.((task.ge.400).and.(task.le.499))) &
	  call putemat(iq, ik, .true., trim(fnemat), istl1, istu1, istl2, istu2, &
	  xiou)
  else
     ! o-u/u-o or o-o/u-u band combinations
     if (.not.(task.eq.430)) then
	allocate(xiuo(nst1, nst2, ngq(iq)))
	xiuo(:, :, :)=xiou(:, :, :)
     end if
     call ematbdlims(2*input%xs%emattype-1, nst1, istl1, istu1, nst2, istl2, istu2)
     istl3=istl2; istu3=istu2; istl4=istl1; istu4=istu1
     nst3=nst2; nst4=nst1     
     if (allocated(xiou)) deallocate(xiou)
     allocate(xiou(nst1, nst2, ngq(iq)))
     call ematqk(iq, ik)
     if (.not.tscreen) &
	  call putemat(iq, ik, .true., trim(fnemat), istl1, istu1, istl2, istu2, &
	  xiou, istl3, istu3, istl4, istu4, xiuo)
  end if
end subroutine ematqk1
