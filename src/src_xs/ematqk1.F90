! Copyright(C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine ematqk1(iq, ik)
  use mod_misc, only: task
  use mod_eigenvalue_occupancy, only: nstsv
  use modinput, only: input
  use modxs, only: nst1, nst2, nst3, nst4,&
                 & istl1, istu1, istl2, istu2,&
                 & istl3, istu3, istl4, istu4,&
                 & xiou, xiuo, ngq, fnemat,&
                 & tscreen
  use m_putemat
  
  implicit none

  ! Arguments
  integer, intent(in) :: iq, ik

  ! Set band combinations
  ! Task 430 is 'screen'
  if( .not. (task .eq. 430)) then

    call ematbdlims(2*input%xs%emattype, nst1, istl1, istu1, nst2, istl2, istu2)
    if(allocated(xiou)) deallocate(xiou)
    if(allocated(xiuo)) deallocate(xiuo)
    allocate(xiou(nst1, nst2, ngq(iq)))

    call ematqk(iq, ik)

  end if

  if(input%xs%emattype .eq. 0) then

    ! All band combinations
    nst3 = nstsv
    nst4 = nstsv
    if( .not. ((task .ge. 400) .and. (task .le. 499))) then
      call putemat(iq, ik, .true., trim(fnemat), istl1, istu1, istl2, istu2, xiou)
    end if

  else

    ! o-u/u-o or o-o/u-u band combinations
    if( .not. (task .eq. 430)) then
      allocate(xiuo(nst1, nst2, ngq(iq)))
      xiuo(:, :, :) = xiou(:, :, :)
    end if

    ! When coming from sreen->df->dfq emattype is set to 1 
    ! and the following sets up o-u combinations:
    !   occupied: nst1 = sto1-sta1+1, istl1 = sta1, istlu1 = sto1
    !   unoccupied: nst2 = sto2-sta2+1, istl2 = istunocc0+sta2-1, istlu2 = istunocc+sto2-1
    call ematbdlims(2*input%xs%emattype-1, nst1, istl1, istu1, nst2, istl2, istu2)
    ! Set index 3 to unoccupied
    istl3 = istl2
    istu3 = istu2
    nst3 = nst2
    ! Set index 4 to occupied
    istl4 = istl1
    istu4 = istu1
    nst4 = nst1

    if(allocated(xiou)) deallocate(xiou)
    allocate(xiou(nst1, nst2, ngq(iq)))

    ! Calculate plane wave elements xiou
    call ematqk(iq, ik)

    if( .not. tscreen) then
      call putemat(iq, ik, .true., trim(fnemat),&
        & istl1, istu1, istl2, istu2, xiou, istl3, istu3, istl4, istu4, xiuo)
    end if

  end if
end subroutine ematqk1
