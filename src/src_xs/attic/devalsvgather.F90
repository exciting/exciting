
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine devalsvgather
  use modmain
  use modxs
  use modmpi
  use m_filedel
  use m_getdevalsv
  use m_putdevalsv
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'devalsvgather'
  integer :: iq,ik,iproc
  real(8) :: vkloff_save(3)

  ! save k-point offset
  vkloff_save = vkloff

  ! allocate matrix elements array
  if (allocated(deou)) deallocate(deou)
  if (allocated(deuo)) deallocate(deuo)
  if (allocated(docc12)) deallocate(docc12)
  if (allocated(docc21)) deallocate(docc21)

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
     allocate(deou(nstval,nstcon))
     allocate(deuo(nstcon,nstval))
     allocate(docc12(nst1,nst2))
     allocate(docc21(nst2,nst1))
     ! file extension for q-point
     do iproc=0,procs-1
        call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=iproc,&
             filnam=fndevalsv_t)
        kpari=firstofset(iproc,nkpt)
        kparf=lastofset(iproc,nkpt)
        do ik=kpari,kparf
           ! exponential factor matrix elements
           call getdevalsv(iq,ik,.false.,trim(fndevalsv_t), &
                deou,docc12,deuo,docc21)
           call putdevalsv(iq,ik,.true.,trim(fndevalsv),deou,docc12,deuo,docc21)
        end do
     end do
     do iproc=0,procs-1
        call genfilname(basename='DEVALSV',iq=iq,procs=procs,rank=iproc,&
             filnam=fndevalsv_t)
        call filedel(trim(fndevalsv_t))
     end do

     deallocate(deou,deuo)
     deallocate(docc12,docc21)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn-Sham eigenvalue &
          &differences gathered for q-point:',iq
  end do

  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)

end subroutine devalsvgather
