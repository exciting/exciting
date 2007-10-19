

!BOP
!
! !ROUTINE: tetcw 
!
! !INTERFACE:
subroutine tetcw(nik,nt,nb,ebd,tetc,linkt,v,efer,omeg,sigfreq,rtyp,cw)
  !
  ! !DESCRIPTION: 
  ! 
  ! This subroutine caluclates the value of the weights for convolution
  ! using the improved tetrahedron method. cw(ik,ib,jb) tells us for a 
  ! given q, the weight on one k point for two bands ib and jb. 
  !  

  ! !USES:

  use tetra_internal

  implicit none      

  ! !INPUT PARAMETERS:

  integer(4), intent(in) :: nik        ! Number of k-points

  integer(4), intent(in) :: nt         ! Number of tetrahedra

  integer(4), intent(in) :: nb         ! Number of bands

!!$  real(8), intent(in) :: ebd(nik,nb)   ! Band energies
!!$
!!$  integer(4), intent(in) :: tetc(4,*) ! id. numbers of the
!!$  !                                                    corners of the 
!!$  !                                                    tetrahedra
!!$
!!$  integer(4), intent(in) :: linkt(*)  ! index to tell with which
!!$  !                                              tetrahedron in k-q mesh is the
!!$  !                                              tetrahedron is k mesh linked by 
!!$  !                                              this q vector 

  real(8), target, intent(in) :: ebd(nik,nb)   ! Band energies

  integer(4), target, intent(in) :: tetc(4,*) ! id. numbers of the
  !                                                    corners of the 
  !                                                    tetrahedra

  integer(4), target, intent(in) :: linkt(*)  ! index to tell with which
  !                                              tetrahedron in k-q mesh is the
  !                                              tetrahedron is k mesh linked by 
  !                                              this q vector 

  real(8), intent(in)    :: v        ! the volume of the tetrahedra

  real(8), intent(in)    :: efer     ! fermi energy

  real(8), intent(in)    :: omeg     ! the fequency for which the
  !                                           weights are calculated, 
  !                                           only used ofr sigfreq =2,3 and
  !                                           for sigfreq=4 which means the 
  !                                           surface integration

  integer(4), intent(in) :: sigfreq  ! Select the kind of bulk convolution
  !                                           weights when it equals 1,2,3. And
  !                                           surface integration for 4.
  ! additional interface parameter resonance type
  integer, intent(in) :: rtyp
  ! !OUTPUT PARAMETERS:

  real(8), intent(out)    :: cw(nik,nik,nb,nb) ! the values of
  !                                                    convolution weights


  ! !LOCAL VARIABLES:
  ! <contribution action="changes">
  real(8), target, allocatable :: target_ebd(:,:)
  integer(4), target, allocatable:: target_tetc(:,:)
  integer(4), target, allocatable :: target_linkt(:)
  ! </contribution>

  real(8), allocatable :: cweit(:,:,:)

  ! !SYSTEM SUBROUTINES:


  external convw

  ! !REVISION HISTORY:
  !
  !   Created: 4th. March 2004 by RGA
  !
  !EOP
  !BOC 

  nirkp = nik
  ntet  = nt
  nband=nb
  vt = v
  allocate(cweit(nirkp,nband,nband))

  ! <contribution action="changes">
  restype=rtyp

  ! additional targets to get around with core dumps in combination with
  ! the Portland compiler
  allocate(target_tetc(1:4,1:ntet))
  allocate(target_linkt(1:ntet))
  allocate(target_ebd(1:nik,1:nb))
  target_tetc(1:4,1:ntet)=tetc(1:4,1:ntet)
  target_linkt(1:ntet)=linkt(1:ntet)
  target_ebd(1:nik,1:nb)=ebd(1:nik,1:nb)

  tetcorn => target_tetc(1:4,1:ntet)
  tetln => target_linkt(1:ntet)
  eband   => target_ebd(1:nik,1:nb)
  ! </contribution>

  !
  !     convw is called to calculate the weights.
  !
  call convw(efer,omeg,sigfreq,cw)

  ! <contribution>
  deallocate(target_tetc,target_linkt,target_ebd)
  ! </contribution>


end subroutine tetcw
      
!EOC
