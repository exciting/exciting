
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematgather2
  use modmain
  use modxs
  use modmpi
  use m_filedel
  use m_getemat2
  use m_putemat2
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='ematgather2'
  integer :: iq,ik,iproc
  ! allocate matrix elements array
  if (allocated(xiou)) deallocate(xiou)
  if (emattype.ne.0) then
     if (allocated(xiuo)) deallocate(xiuo)
  end if
  ! loop over q-points
  do iq=1,nqpt
     ! find highest (partially) occupied and lowest (partially) unoccupied
     ! states
     call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0, &
          istu)
     ! set limits for band combinations
     call ematbdcmbs(emattype)
     ! calculate k+q and G+k+q related variables
     call init1xs(qvkloff(1,iq))
     allocate(xiou(nst1,nst2,ngq(iq)))
     if (emattype.ne.0) allocate(xiuo(nst3,nst4,ngq(iq)))
     ! file extension for q-point
     do iproc=0,procs-1
        call genfilname(basename='EMAT',iq=iq,procs=procs,rank=iproc,&
             filnam=fnemat_t)
        kpari=firstofset(iproc,nkpt)
        kparf=lastofset(iproc,nkpt)
        do ik=kpari,kparf
           if (emattype.ne.0) then
              call getemat2(iq,ik,.false.,trim(fnemat_t),x1=xiou,x2=xiuo)
              call putemat2(iq,ik,.true.,trim(fnemat),x1=xiou,x2=xiuo)
           else
              call getemat2(iq,ik,.false.,trim(fnemat_t),x1=xiou)
              call putemat2(iq,ik,.true.,trim(fnemat),x1=xiou)
           end if
        end do
     end do
     do iproc=0,procs-1
        call genfilname(basename='EMAT',iq=iq,procs=procs,rank=iproc,&
             filnam=fnemat_t)
        call filedel(trim(fnemat_t))
     end do
     deallocate(xiou)
     if (emattype.ne.0) deallocate(xiuo)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Matrix elements of &
          &exponential factor gathered for q-point:',iq
  end do
end subroutine ematgather2
