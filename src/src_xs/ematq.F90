

! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine ematq(iq)
  use modmain
use modinput
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
  call genfilname(basename='EMAT', iqmt=iq, etype=input%xs%emattype, filnam=fnemat)
  call genfilname(basename = 'EMAT', iqmt = iq, etype = input%xs%emattype, procs = procs, rank = rank, &
       filnam = fnemat_t)
  call genfilname(nodotpar = .true., basename = 'EMAT_TIMING', iqmt = iq, &
       etype = input%xs%emattype, procs = procs, rank = rank, filnam = fnetim)
  ! file extension for q-point
  call genfilname(iqmt=iq, setfilext=.true.)
  ! calculate k+q and G+k+q related variables
  call init1offs(qvkloff(1, iq))
  ! write G+q-vectors
  if (rank.eq.0) then
     call writegqpts(iq, filext)
     call writekmapkq(iq)
  end if
  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(iq, istocc0, istocc, istunocc0, istunocc, isto0, isto, istu0, istu)
  call ematbdlims(1, nst1, istl1, istu1, nst2, istl2, istu2)
  ! generate radial integrals wrt. sph. Bessel functions
  call ematrad(iq)
  ! delete timing information of previous runs
  call filedel(trim(fnetim))
  ! write information
  write(unitout, '(a, i6)') 'Info('//thisnam//'): number of G+q vectors:', &
       ngq(iq)
  call ematqalloc
  ! loop over k-points
  do ik=kpari, kparf
     call ematqk1(iq, ik)
  end do
  call ematqdealloc
  call barrier
end subroutine ematq
