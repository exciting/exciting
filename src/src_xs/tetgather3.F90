
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tetgather3
  use modmain
  use modxs
  use modmpi
  use m_gettetcw3
  use m_puttetcw3
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='tetgather3'
  character(256) :: filnam,filnam_t
  integer :: iq,iproc,ik,i1,i2,nwdfp
  real(8), allocatable :: cw(:),cwa(:),cwsurf(:)
  real(8), allocatable :: cwp(:),cwap(:),cwsurfp(:)
  allocate(cw(nwdf),cwa(nwdf),cwsurf(nwdf))
  ! loop over q-points
  do iq=1,nqpt
     ! calculate k+q and G+k+q related variables
     call init1xs(qvkloff(1,iq))
     ! file name for output file
     call genfilname(basename='TETW',iqfmt=iq,filnam=filnam)
     do ik=1,nkpt
        do i1=1,nst1
           do i2=1,nst2
              ! collect weights from processes
              do iproc=0,procs-1
                 ! filename for input file
                 call genfilname(basename='TETW',iqfmt=iq,rank=iproc,&
                      procs=procs,filnam=filnam_t)                 
                 wpari=firstofset(iproc,nwdf)
                 wparf=lastofset(iproc,nwdf)
                 nwdfp=wparf-wpari+1
                 allocate(cwp(nwdfp),cwap(nwdfp),cwsurfp(nwdfp))
                 call gettetcw3(iq,ik,i1,i2,nwdfp,trim(filnam_t),&
                      cwp,cwap,cwsurfp)
                 cw(wpari:wparf)=cwp(:)
                 cwa(wpari:wparf)=cwap(:)
                 cwsurf(wpari:wparf)=cwsurfp(:)
                 deallocate(cwp,cwap,cwsurfp)
              end do ! iproc
              ! write weights
              call puttetcw3(iq,ik,i1,i2,trim(filnam),cw,cwa,cwsurf)
           end do
        end do
        ! end loop over k-points
     end do
     do iproc=0,procs-1
        call genfilname(basename='TETW',iqfmt=iq,rank=rank,procs=procs,&
             filnam=filnam_t)
        call filedel(trim(filnam_t))
     end do
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method gathered for q-point:',iq
     ! end loop over q-points
  end do
  deallocate(cw,cwa,cwsurf)
end subroutine tetgather3
