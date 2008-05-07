
! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematq(iq)
  use modmain
  use modxs
  use modmpi
  use m_writegqpts
  use m_filedel
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='ematq'
  integer :: ik
  ! filenames
  call genfilname(basename='EMAT',iqmt=iq,etype=emattype,filnam=fnemat)
  call genfilname(basename='EMAT',iqmt=iq,etype=emattype,procs=procs,rank=rank,&
       filnam=fnemat_t)
  call genfilname(nodotpar=.true.,basename='EMAT_TIMING',iqmt=iq,&
       etype=emattype,procs=procs,rank=rank,filnam=fnetim)
  ! file extension for q-point
  call genfilname(iqmt=iq,setfilext=.true.)
  ! calculate k+q and G+k+q related variables
  call init1xs(qvkloff(1,iq))
  ! write G+q-vectors
  if (rank.eq.0) then
     call writegqpts(iq)
     call writekmapkq(iq)
  end if
  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(iq,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  call ematbdlims(1,nst1,istlo1,isthi1,nst2,istlo2,isthi2)
  ! generate radial integrals wrt. sph. Bessel functions
  call ematrad(iq)
  ! delete timing information of previous runs
  call filedel(trim(fnetim))
  ! write information
  write(unitout,'(a,i6)') 'Info('//thisnam//'): number of G+q vectors:', &
       ngq(iq)
  call ematqalloc
  ! loop over k-points
  do ik=kpari,kparf
     call ematqk1(iq,ik)
     ! synchronize for common number of k-points to all processes
     if ((partype.eq.'k').and.(ik-kpari+1 <= nkpt/procs)) call barrier
     ! end loop over k-points
  end do
  call ematqdealloc
end subroutine ematq
