! Copyright(C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematqk1
! !INTERFACE:
subroutine ematqk1(iq, ik)
! !USES:
  use mod_misc, only: task
  use mod_eigenvalue_occupancy, only: nstsv
  use modinput, only: input
  use modxs, only: nst1, nst2, nst3, nst4,&
                 & istl1, istu1, istl2, istu2,&
                 & istl3, istu3, istl4, istu4,&
                 & xiou, xiuo, ngq, fnemat,&
                 & tscreen
  use m_putemat
! !INPUT/OUTPUT PARAMETERS:
!   IN:
!   iq, integer : Q-point index
!   ik, integer : K-point index
!   INDIRECT I/O:
!   Manipulates modxs:xiou and modxs:xiuo
!
! !DESCRIPTION:
! This routine is a wrapper routine for {\tt ematqk}. Depending
! on the selected task (e.g. 'screen' or not 'screen') band combinations
! are chosen, {\tt ematqk} is called and its output to the modxs array
! {\tt xiou} is managed.
!
! !REVISION HISTORY:
! Added to documentation scheme. (Aurich)
!
!EOP
!BOC
  
  implicit none

  ! Arguments
  integer, intent(in) :: iq, ik

  ! Set band combinations
  ! Task 430 is 'screen', task 440 is 'scrcoulint'
  ! When coming from scrcoulint (task 440) emattype is set to 2
  if( .not. (task .eq. 430)) then

    ! For 'scrcoulint' this selects 12=uu combinations
    call ematbdlims(2*input%xs%emattype, nst1, istl1, istu1, nst2, istl2, istu2)
    ! Clear plane wave arrays
    if(allocated(xiou)) deallocate(xiou)
    if(allocated(xiuo)) deallocate(xiuo)
    ! ematqk only operates on xiou, regardless of the chosen bands
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
      ! In the case of 'scrcoulint' this copies the u-u elements calculated
      ! in the previous call to ematqk to the modxs array xiuo.
      allocate(xiuo(nst1, nst2, ngq(iq)))
      xiuo(:, :, :) = xiou(:, :, :)
    end if

    ! When coming from sreen->df->dfq ('screen') emattype is set to 1 
    ! and the following sets up o-u combinations:
    !   occupied: nst1 = sto1-sta1+1, istl1 = sta1, istlu1 = sto1
    !   unoccupied: nst2 = sto2-sta2+1, istl2 = istunocc0+sta2-1, istlu2 = istunocc+sto2-1
    ! In the case of 'scrcoulint' emattype is set to 2, so that the following
    ! sets the band combinations to 12=oo
    call ematbdlims(2*input%xs%emattype-1, nst1, istl1, istu1, nst2, istl2, istu2)
    ! 'screen': Set index 3 to unoccupied
    ! 'scrcoulint': Set index 3 to occupied
    istl3 = istl2
    istu3 = istu2
    nst3 = nst2
    ! 'screen': Set index 4 to occupied
    ! 'scrcoulint': Set index 4 to occupied
    istl4 = istl1
    istu4 = istu1
    nst4 = nst1

    ! 'scrcoulint': Set up xiou for oo band ranges
    if(allocated(xiou)) deallocate(xiou)
    allocate(xiou(nst1, nst2, ngq(iq)))

    ! Calculate plane wave elements xiou
    ! 'scrcoulint': Calculate the oo plane wave elements.
    call ematqk(iq, ik)
    ! To summarize:
    !   'screen': only o-u elements are calculated and stored in xiou
    !   'scrcoulint': u-u and o-o elements are calculated and stored in 
    !                 xiuo and xiou respectively

    if( .not. tscreen) then
      call putemat(iq, ik, .true., trim(fnemat),&
        & istl1, istu1, istl2, istu2, xiou, istl3, istu3, istl4, istu4, xiuo)
    end if

  end if
end subroutine ematqk1
!EOC
