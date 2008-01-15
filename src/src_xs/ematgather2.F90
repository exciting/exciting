
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
  character(*), parameter :: thisnam = 'ematgather2'
  integer :: iq,ik,iproc
  real(8) :: vkloff_save(3)

  ! save k-point offset
  vkloff_save(:)=vkloff(:)

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
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1xs
     allocate(xiou(nstval,nstcon,ngq(iq)))
     if (emattype.ne.0) allocate(xiuo(nstcon,nstval,ngq(iq)))
     ! file extension for q-point
     do iproc=0,procs-1
        call genfilname(basename='EMAT',iq=iq,procs=procs,rank=iproc,&
             filnam=fnemat_t)
        kpari=firstofset(iproc,nkpt)
        kparf=lastofset(iproc,nkpt)
        do ik=kpari,kparf
           ! exponential factor matrix elements
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

  ! restore offset
  vkloff(:)=vkloff_save(:)
  call genfilname(setfilext=.true.)

end subroutine ematgather2
