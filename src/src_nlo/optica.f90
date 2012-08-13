module optica

  real*8  :: idel1,dw,sc,tol
  integer :: banmin,banmax,emesh

  complex*8, allocatable :: r(:,:,:,:)          ! matrix elements all
  integer,   allocatable :: noval(:)
  integer,   allocatable :: nocond(:)
  integer,   allocatable :: sym(:,:,:)
  integer,   allocatable :: tot_ban(:)
  
  integer v1,v2,v3                              ! component calculated

! for error messages
  character*80 title
  character*4  lattic
  character*6  modus

end module optica
