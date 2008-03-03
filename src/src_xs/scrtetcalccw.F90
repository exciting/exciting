subroutine scrtetcalccw
  use modmain
  use modxs
  use m_genfilname
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
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  ! calculate tetrahedron weights with only one frequency point
  call tetcalccw
  ! restore global variables
  nosym=nosymt
  reducek=reducekt
  ngridk(:)=ngridkt(:)
  vkloff(:)=vklofft(:)
  rgkmax=rgkmaxt
  nempty=nemptyt
  nwdf=nwdft
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
end subroutine scrtetcalccw
