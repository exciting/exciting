
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findkmapkq(iq,vq,voff,map)
  use modmain
  use modxs
  use modmpi
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  real(8), intent(in) :: vq(3),voff(3)
  integer, intent(out) :: map(nkpt)
  ! local variables
  integer :: ik,un,iv(3),ivt(3)
  real(8) :: vkq(3)
  real(8), external :: r3taxi
  character(256) :: filnam

!!$    if (task.ge.400) then
  write(*,*) '............... new method in findkmapkq'
  do ik=1,nkpt
     vkq(:)=vkl(:,ik)+vq(:)
     call r3frac(epslat,vkq,ivt)
     iv(:)=nint(vkq(:)*ngridk(:)-voff(:))
     map(ik)=ikmap(iv(1),iv(2),iv(3))
  end do
!!$
!!$    else
!!$
!!$       allocate(done(nkptnr))
!!$       done(:)=.false.
!!$       map(:)=0
!!$
!!$       vqt=vq
!!$       vofft=voff-vkloff
!!$       call mapkto01(vqt)
!!$       do ik=1,nkpt
!!$          ! k+q point from k-point
!!$          vkq(:)=vkl(:,ik)+vqt(:)
!!$          call mapkto01(vkq)
!!$          do ikt=1,nkptnr
!!$             if (.not.done(ikt)) then
!!$                vkqt(:)=vklnr(:,ikt)+vofft(:)/dble(ngridk)
!!$                if (r3taxi(vkq,vkqt).lt.epslat) then
!!$                   iktr=ikmap(ivknr(1,ikt),ivknr(2,ikt),ivknr(3,ikt))
!!$                   done(ikt)=.true.
!!$                   map(ik)=iktr
!!$                   exit
!!$                end if
!!$             end if
!!$          end do
!!$          if (map(ik).eq.0) then
!!$             write(*,*) 'Error(findkmapkq): mapping between k and k+q point set &
!!$                  &failed for k-point:',ik
!!$             call terminate
!!$          end if
!!$       end do
!!$       deallocate(done)
!!$
!!$    end if

  if (rank.eq.0) then
     call getunit(un)
     call genfilname(basename='KMAPKQ',iq=iq,filnam=filnam)
     open(un,file=trim(filnam),form='formatted',action='write', &
          status='replace')
     write(un,'(i9,a)') nkpt, ' : nkpt; k-point, ikmapikq below'
     do ik=1,nkpt
        write(un,'(2i9)') ik,map(ik)
     end do
     close(un)
  end if

end subroutine findkmapkq

subroutine mapkto01(v)
  implicit none
  ! arguments
  real(8), intent(inout) :: v(3)
  ! local variables
  !integer :: id(3)
  integer(8) :: v2(3),v3(3)
  real(8),parameter :: fac=1.d15

  !call r3frac(epslat,v,id)

  v2=dint(v)
  v3=dint(fac*v)
  v3=v3-v2*dint(fac)
  v=dble(v3/dint(fac))
  where(v.lt.0.d0) v=v+1.d0

end subroutine mapkto01

