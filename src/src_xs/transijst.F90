
! Copyright (C) 2007-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: transijst
! !INTERFACE:
logical function transijst(ik,ist,jst)
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
  integer, intent(in) :: ik,ist,jst
  integer, parameter :: infinity=1000000000
  integer :: itrans,istat,jstat,irange,jrange
  integer :: ikr,jkr,istr,jstr
  integer :: istar,istor,jstar,jstor
  logical :: tex,texi,texj
  transijst=.true.
  if (associated(input%xs%transitions)) then
    ! start from scratch (no transitions)
    transijst=.false.
    ! individual transitions
    if (associated(input%xs%transitions%individual)) then
      do itrans=1,size(input%xs%transitions%individual%transarray)
        ! get include/exclude attribute
        tex=input%xs%transitions%individual%transarray(itrans)%trans%action .eq. "exclude"
        ikr=input%xs%transitions%individual%transarray(itrans)%trans%kpointnumber
        istr=input%xs%transitions%individual%transarray(itrans)%trans%initial
        jstr=input%xs%transitions%individual%transarray(itrans)%trans%final
        if ((ikr .eq. 0).or.(ikr .eq. ik)) then
          if (((istr .eq. 0) .or. (istr .eq. ist)) .and. ((jstr .eq. 0) .or. (jstr .eq. jst))) then
            transijst=.not. tex
          end if
        end if
      end do
    end if
    ! lists
    if (associated(input%xs%transitions%lists)) then
      do istat=1,size(input%xs%transitions%lists%istatearray)
        if (input%xs%transitions%lists%istatearray(istat)%istate%statestype .eq. "initialstates") then
          texi=input%xs%transitions%lists%istatearray(istat)%istate%action .eq. "exclude"
          ikr=input%xs%transitions%lists%istatearray(istat)%istate%kpointnumber
          istr=input%xs%transitions%lists%istatearray(istat)%istate%state
          do jstat=1,size(input%xs%transitions%lists%istatearray)
            if (input%xs%transitions%lists%istatearray(jstat)%istate%statestype .eq. "finalstates") then
              texj=input%xs%transitions%lists%istatearray(jstat)%istate%action .eq. "exclude"
              jkr=input%xs%transitions%lists%istatearray(jstat)%istate%kpointnumber
              jstr=input%xs%transitions%lists%istatearray(jstat)%istate%state
              ! vertical transitions between k-point and (k+q)-point sets
              if (((ik .eq. ikr).or.(ikr .eq. 0)).and.((ik .eq. jkr).or.(jkr .eq. 0))) then
                if (((istr .eq. 0) .or. (istr .eq. ist)) .and. ((jstr .eq. 0) .or. (jstr .eq. jst))) then
                  transijst=(.not. texi) .and. (.not. texj)
                end if
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
          istar=input%xs%transitions%ranges%rangearray(irange)%range%start
          istor=input%xs%transitions%ranges%rangearray(irange)%range%stop
          ! no upper limit for final states in this case of considering all final states
          if (istor .eq. 0) istor=infinity
          do jrange=1,size(input%xs%transitions%ranges%rangearray)
            if (input%xs%transitions%ranges%rangearray(jrange)%range%statestype .eq. "finalstates") then
              texj=input%xs%transitions%ranges%rangearray(jrange)%range%action .eq. "exclude"
              jkr=input%xs%transitions%ranges%rangearray(jrange)%range%kpointnumber
              jstar=input%xs%transitions%ranges%rangearray(jrange)%range%start
              jstor=input%xs%transitions%ranges%rangearray(jrange)%range%stop
              ! no upper limit for final states in this case of considering all final states
              if (jstor .eq. 0) jstor=infinity
              ! vertical transitions between k-point and (k+q)-point sets
              if (((ik .eq. ikr).or.(ikr .eq. 0)).and.((ik .eq. jkr).or.(jkr .eq. 0))) then
                if (((ist .ge. istar).and.(ist .le. istor)).and.((jst .ge. jstar).and.(jst .le. jstor))) then
                  transijst=(.not. texi) .and. (.not. texj)
                end if
              end if
            end if
          end do
        end if
      end do
    end if
  end if
end function
