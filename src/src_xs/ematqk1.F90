
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqk1(iq,ik)
  use modmain
  use modxs
  use modmpi
  use m_putemat
  implicit none
  ! arguments
  integer, intent(in) :: iq,ik
  ! set band combinations
  if (.not.(task.eq.430)) then
     call ematbdlims(2*emattype,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     if (allocated(xiou)) deallocate(xiou)
     if (allocated(xiuo)) deallocate(xiuo)
     allocate(xiou(nst1,nst2,ngq(iq)))
     call ematqk(iq,ik)
  end if
  if (emattype.eq.0) then
     ! all band combinations
     nst3=nstsv; nst4=nstsv
     if (.not.((task.ge.400).and.(task.le.499))) &
          call putemat(iq,ik,.true.,trim(fnemat),istlo1,isthi1,istlo2,isthi2, &
          xiou)
  else
     ! o-u/u-o or o-o/u-u band combinations
     if (.not.(task.eq.430)) then
        allocate(xiuo(nst1,nst2,ngq(iq)))
        xiuo(:,:,:)=xiou(:,:,:)
     end if
     call ematbdlims(2*emattype-1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     istlo3=istlo2; isthi3=isthi2; istlo4=istlo1; isthi4=isthi1
     nst3=nst2; nst4=nst1     
     if (allocated(xiou)) deallocate(xiou)
     allocate(xiou(nst1,nst2,ngq(iq)))
     call ematqk(iq,ik)
     if (.not.tscreen) &
          call putemat(iq,ik,.true.,trim(fnemat),istlo1,isthi1,istlo2,isthi2, &
          xiou,istlo3,isthi3,istlo4,isthi4,xiuo)
  end if
end subroutine ematqk1
