
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine scrgeneigvec
  use modmain
  use modmpi
  use modxs
  use m_filedel
  use m_genfilname
  use m_writegqpts
  implicit none

!!$  ! local variables
!!$  character(*), parameter :: thisnam='scrgeneigvec'
!!$  real(8) :: vklofft(3),rgkmaxt
!!$  integer :: ngridkt(3),nemptyt
!!$  logical :: nosymt,reducekt

!!$  call init0
!!$  call init1
!!$  ! save global variables
!!$  nosymt=nosym
!!$  reducekt=reducek
!!$  ngridkt(:)=ngridk(:)
!!$  vklofft(:)=vkloff(:)
!!$  rgkmaxt=rgkmax
!!$  nemptyt=nempty
!!$  nosym=nosymscr
!!$  ! no symmetries implemented for screening
!!$  reducek=.false.
!!$  ngridk(:)=ngridkscr(:)
!!$  vkloff(:)=vkloffscr(:)
!!$  rgkmax=rgkmaxscr
!!$  nempty=nemptyscr

  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  ! generate eigenvectors, eigenvalues, occupancies and APW MT coefficients
  call xsgeneigvec
  write(unitout,'("Info(scrgeneigvec): eigenvectors for screening finished")')

!!$  ! restore global variables
!!$  nosym=nosymt
!!$  reducek=reducekt
!!$  ngridk(:)=ngridkt(:)
!!$  vkloff(:)=vklofft(:)
!!$  rgkmax=rgkmaxt
!!$  nempty=nemptyt
end subroutine scrgeneigvec
