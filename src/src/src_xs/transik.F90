
! Copyright (C) 2007-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: transik
! !INTERFACE:
logical function transik(ik)
! !USES:
  use modinput
! !DESCRIPTION:
!   This function returns "true" if the transition specified by the input k-point
!   initial and final state matches with the ones specified in the {\tt transitions}
!   element. Wether a state is included or excluded depends on the sequence
!   of the transitions specified (using the "include" and "exclude" attribute).
!
! !REVISION HISTORY:
!   Created July 2010 (Sagmeister)
!EOP
!BOC
  implicit none
  integer, intent(in) :: ik
  integer :: itrans,istat,jstat,irange,jrange
  integer :: ikr,jkr
  logical :: tex,texi,texj
  transik=.true.
  if (associated(input%xs%transitions)) then
    ! start from scratch (no transitions)
    transik=.false.
    ! individual transitions
    if (associated(input%xs%transitions%individual)) then
      do itrans=1,size(input%xs%transitions%individual%transarray)
        ! get include/exclude attribute
        tex=input%xs%transitions%individual%transarray(itrans)%trans%action .eq. "exclude"
        ! if in exclude mode return drop the exception, since not all transitions
        ! might be excluded for the current k-point
        if (tex) cycle
        ikr=input%xs%transitions%individual%transarray(itrans)%trans%kpointnumber
        if ((ikr .eq. 0).or.(ikr .eq. ik)) then
          transik=.not. tex
        end if
      end do
    end if
    ! lists
    if (associated(input%xs%transitions%lists)) then
      do istat=1,size(input%xs%transitions%lists%istatearray)
        if (input%xs%transitions%lists%istatearray(istat)%istate%statestype .eq. "initialstates") then
          texi=input%xs%transitions%lists%istatearray(istat)%istate%action .eq. "exclude"
          ikr=input%xs%transitions%lists%istatearray(istat)%istate%kpointnumber
          do jstat=1,size(input%xs%transitions%lists%istatearray)
            if (input%xs%transitions%lists%istatearray(jstat)%istate%statestype .eq. "finalstates") then
              texj=input%xs%transitions%lists%istatearray(jstat)%istate%action .eq. "exclude"
              ! if in exclude mode return drop the exception, since not all transitions
              ! might be excluded for the current k-point
              if (texi .or. texj) cycle
              jkr=input%xs%transitions%lists%istatearray(jstat)%istate%kpointnumber
              ! vertical transitions between k-point and (k+q)-point sets
              if (((ik .eq. ikr).or.(ikr .eq. 0)).and.((ik .eq. jkr).or.(jkr .eq. 0))) then
                transik=(.not. texi) .and. (.not. texj)
              end if
            end if
          end do
        end if
      end do
    end if
    ! ranges
    if (associated(input%xs%transitions%ranges)) then
      do irange=1,size(input%xs%transitions%ranges%rangearray)
        if (input%xs%transitions%ranges%rangearray(irange)%range%statestype .eq. "initialstates") then
          texi=input%xs%transitions%ranges%rangearray(irange)%range%action .eq. "exclude"
          ikr=input%xs%transitions%ranges%rangearray(irange)%range%kpointnumber
          do jrange=1,size(input%xs%transitions%ranges%rangearray)
            if (input%xs%transitions%ranges%rangearray(jrange)%range%statestype .eq. "finalstates") then
              texj=input%xs%transitions%ranges%rangearray(jrange)%range%action .eq. "exclude"
              ! if in exclude mode return drop the exception, since not all transitions
              ! might be excluded for the current k-point
              if (texi .or. texj) cycle
              jkr=input%xs%transitions%ranges%rangearray(jrange)%range%kpointnumber
              ! vertical transitions between k-point and (k+q)-point sets
              if (((ik .eq. ikr).or.(ikr .eq. 0)).and.((ik .eq. jkr).or.(jkr .eq. 0))) then
                transik=(.not. texi) .and. (.not. texj)
              end if
            end if
          end do
        end if
      end do
    end if
  end if
end function
