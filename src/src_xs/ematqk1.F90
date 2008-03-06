
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
  call ematbdlims(2*emattype,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
  if (allocated(xiou)) deallocate(xiou)
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiou(nst1,nst2,ngq(iq)))
  call ematqk(iq,ik)
  if (emattype.eq.0) then
     ! all band combinations
     nst3=nstsv; nst4=nstsv
     if (.not.((task.ge.400).and.(task.le.499))) &
          call putemat(iq,ik,.false.,trim(fnemat_t),x1=xiou)
  else
     ! v-c/c-v or v-v/c-c band combinations
     allocate(xiuo(nst1,nst2,ngq(iq)))
     xiuo(:,:,:)=xiou(:,:,:)
     deallocate(xiou)
     call ematbdlims(2*emattype-1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
     nst3=nst2; nst4=nst1
     allocate(xiou(nst1,nst2,ngq(iq)))
     call ematqk(iq,ik)
     if (.not.((task.ge.400).and.(task.le.499))) &
          call putemat(iq,ik,.false.,trim(fnemat_t),x1=xiou,x2=xiuo)
  end if
!  deallocate(xiou)
!  if (allocated(xiuo)) deallocate(xiuo)
end subroutine ematqk1
