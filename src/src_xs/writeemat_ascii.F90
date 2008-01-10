
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeemat_ascii
  use modmain
  use modxs
  use m_getunit
  use m_getemat2
  use m_genfilname
  implicit none
  character(256) :: filnam
  integer :: un,iq,ik,i,j,ib,jb,igq
  complex(8) :: zt
  real(8) :: vkloff_save(3)
  ! save k-point offset
  vkloff_save = vkloff
  call init0
  call init1
  call tdsave0
  call init2xs
  call getunit(un)
  ! loop over q-points
  do iq=1,nqpt
     ! find highest (partially) occupied and lowest (partially) unoccupied
     ! states
     call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0, &
          istu)
     ! set limits for band combinations
     call ematbdcmbs(emattype)
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1xs
     if (allocated(xiou)) deallocate(xiou)
     allocate(xiou(nst1,nst2,ngq(iq)))
     if (emattype.ne.0) then
        if (allocated(xiuo)) deallocate(xiuo)
        allocate(xiuo(nst3,nst4,ngq(iq)))
     end if
     ! filename for matrix elements file
     call genfilname(basename='EMAT',asc=.true.,iq=iq,etype=emattype, &
          filnam=filnam)
     open(un,file=trim(filnam),action='write')
     ! read matrix elements of exponential expression
     call genfilname(basename='EMAT',iq=iq,etype=emattype,filnam=fnemat)
     write(un,'(a)') 'iq,ik,igq,i1,i2,emat,|emat|^2, below'
     ! loop over k-points
     do ik=1,nkpt
        if (emattype.eq.0) then
           call getemat2(iq,ik,trim(fnemat),x1=xiou)
        else
           call getemat2(iq,ik,trim(fnemat),x1=xiou,x2=xiuo)
        end if
        do igq=1,ngq(iq)
           do i=1,nst1
              ib=i+istlo1-1
              do j=1,nst2
                 jb=j+istlo2-1
                 zt=xiou(i,j,igq)
                 write(un,'(5i8,3g18.10)') iq,ik,igq,ib,jb,zt,abs(zt)**2
              end do
           end do
        end do
        do igq=1,ngq(iq)
           do i=1,nst3
              ib=i+istlo3-1
              do j=1,nst4
                 jb=j+istlo4-1
                 zt=xiuo(i,j,igq)
                 write(un,'(5i8,3g18.10)') iq,ik,igq,ib,jb,zt,abs(zt)**2
              end do
           end do
        end do
     end do ! ik
     close(un)
     deallocate(xiou)
     if (emattype.ne.0) deallocate(xiuo)
  end do ! iq
  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)
end subroutine writeemat_ascii
