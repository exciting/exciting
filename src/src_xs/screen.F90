
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine screen
  use modmain
  use modmpi
  use modxs
  use m_genfilname
  use m_writegqpts
  implicit none
  ! local variables
  character(*), parameter :: thisnam='screen'
  real(8) :: vklofft(3),rgkmaxt
  integer :: ngridkt(3),nemptyt,nwdft
  logical :: nosymt,reducekt
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  rgkmaxt=rgkmax
  nemptyt=nempty
  nwdft=nwdf
  ! map variables for screening
  call initscr
  nosym=nosymscr
  ! no symmetries implemented for screening
  reducek=.false.
  ngridk(:)=ngridkscr(:)
  vkloff(:)=vkloffscr(:)
  rgkmax=rgkmaxscr
  nempty=nemptyscr
  ! only one frequency w=0
  nwdf=1
  emattype=1
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  ! call dielectric function with only one frequency point
  call df
  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  rgkmax=rgkmaxt
  nempty=nemptyt
  nwdf=nwdft
  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"
end subroutine screen
