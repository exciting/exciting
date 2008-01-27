
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine devalsvgather2
  use modmain
  use modxs
  use modmpi
  use m_filedel
  use m_genfilname
  use m_getdevalsv2
  use m_putdevalsv2
  implicit none
  ! local variables
  character(*), parameter :: thisnam='devalsvgather2'
  integer :: iq,ik,iproc
  ! allocate matrix elements array
  if (allocated(deou)) deallocate(deou)
  if (allocated(docc12)) deallocate(docc12)
  if (emattype.ne.0) then
     if (allocated(deuo)) deallocate(deuo)
     if (allocated(docc21)) deallocate(docc21)
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
     allocate(deou(nst1,nst2))
     allocate(docc12(nst1,nst2))
     if (emattype.ne.0) then
        allocate(deuo(nst3,nst4))
        allocate(docc21(nst3,nst4))
     end if
     ! file extension for q-point
     do iproc=0,procs-1
        call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=iproc,&
             filnam=fndevalsv_t)
        kpari=firstofset(iproc,nkpt)
        kparf=lastofset(iproc,nkpt)
        do ik=kpari,kparf
           if (emattype.ne.0) then
              call getdevalsv2(iq,ik,.false.,trim(fndevalsv_t), &
                   deou,docc12,deuo,docc21)
              call putdevalsv2(iq,ik,.true.,trim(fndevalsv), &
                   deou,docc12,deuo,docc21)
           else
              call getdevalsv2(iq,ik,.false.,trim(fndevalsv_t), &
                   deou,docc12)
              call putdevalsv2(iq,ik,.true.,trim(fndevalsv), &
                   deou,docc12)
           end if
        end do
     end do
     do iproc=0,procs-1
        call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=iproc,&
             filnam=fndevalsv_t)
        call filedel(trim(fndevalsv_t))
     end do
     deallocate(deou,docc12)
     if (emattype.ne.0) deallocate(deuo,docc21)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn-Sham eigenvalue &
          &differences gathered for q-point:',iq
  end do
end subroutine devalsvgather2
