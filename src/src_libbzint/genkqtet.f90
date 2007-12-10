
! Copyright (C) 2007 S. Sagmeister, R. Gomez-Abal and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkqtet(nsymcr,symcr,omega,bvec,vql,ngridk,nkptnr,nkpt,vkloff,ivk,&
     ivknr,ikmap,ikmapnr,wkpt,tvol,tnodes,link)
  use kgen_internals
  implicit none
  ! arguments
  integer, intent(in) :: nsymcr
  integer, intent(in) :: symcr(3,3,nsymcr)
  real(8), intent(in) :: omega
  real(8), intent(in) :: bvec(3,3)
  real(8), intent(in) :: vql(3)
  integer, intent(in) :: ngridk(3)
  integer, intent(in) :: nkptnr
  integer, intent(in) :: nkpt
  real(8), intent(in) :: vkloff(3)
  integer, intent(in) :: ivk(3,ngridk(1)*ngridk(2)*ngridk(3))
  integer, intent(in) :: ivknr(3,ngridk(1)*ngridk(2)*ngridk(3))
  integer, intent(in) :: ikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
  integer, intent(in) :: ikmapnr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
  real(8), intent(in) :: wkpt(nkpt)
  real(8), intent(out) :: tvol
  integer, intent(out) :: tnodes(4,6*nkptnr)
  integer, intent(out) :: link(6*nkptnr)
  ! local variables
  integer, allocatable :: map(:)
  integer :: ik,ntet,divsh(3),j,im,jm,tet(3,4,6)
  real(8), parameter :: epslat=1.d-6
  logical :: done(ngridk(1)*ngridk(2)*ngridk(3))

  ! assign symmetries [iio]
  if (allocated(iio)) deallocate(iio)
  allocate(iio(3,3,nsymcr))
  iio(:,:,:)=symcr(:,:,1:nsymcr)
  
  ! assign lattice vectors [gbas]
  gbas(:,:)=bvec(:,:)

  ! number of symmetry-reduced k-points [nirkp]
  nirkp=nkpt
  ! assign grid division [div]
  div(:)=ngridk(:)
  ! determine integer representation of submesh offset [shift],(divsh)
  call factorize(3,vkloff,shift,divsh)
  ! check offset factorization
  if (any(abs(dble(shift)/dble(divsh)-vkloff).gt.epslat)) then
     write(*,*)
     write(*,'("Error(libbzint:genkqtet):")')
     write(*,'(" factorization of k-point offest failed")')
     write(*,'(" offset                   :",3g18.10)') vkloff
     write(*,'(" offset from factorization:",3g18.10)') &
          dble(shift)/dble(divsh)
     write(*,*)
     stop
  end if

  ! allocate arrays for k-points
  if(allocated(ikpid)) deallocate(ikpid)
  allocate(ikpid(nkptnr))
  if(allocated(redkp)) deallocate(redkp)
  allocate(redkp(nkptnr))

  ! [ikpid], [redkp]
write(*,*) 'genkqtet:' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(map(nkptnr))
  done(:)=.false.
  j=0
  do ik=1,nkptnr
     im=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
     jm=im
     if (done(jm)) jm=0
     if (jm.ne.0) then
        done(jm)=.true.
        j=j+1
        map(j)=ik
     end if
     ikpid(ik)=jm
     redkp(ik)=map(im)
write(*,*) ik,ikpid(ik),redkp(ik) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end do
  deallocate(map)

  ! number of tetrahedra
  ntet=6*nkpt
  ! tetrahedra volume [vt]
  vt=omega/6.d0
  tvol=vt
  ! allocate arrays for tetrahedra
  if (allocated(redtet)) deallocate(redtet)
  allocate(redtet(ntet))
        
  ! [mndg] from "tetinit"
  call tetinit(tet)

  ! [redtet]
  redtet(:)=0

  ! (tnodes)
  tnodes(:,:)=0

  ! (link) generate link array
  link(:)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***!!$      integer(4) :: nirkp            ! Number of irreducible k-points
!!$
!***!!$      integer(4), dimension(3) :: div      ! Number of subdivisions of 
!!$!                                            the BZ in each direction
!!$
!***!!$      integer(4), allocatable  :: ikpid(:) ! Order number of the 
!!$!                                            irreducible k-points
!!$                                            
!***!!$      integer(4), allocatable  :: redkp(:) ! Identification number of 
!!$!                                            the irreducible k-point
!!$!                                            associated to the general
!!$!                                            k-point i.
!!$      integer(4), allocatable :: redtet(:)
!!$                                            
!***!!$      integer(4) :: mndg                                      
!!$                                            
!***!!$      integer(4), dimension(3) :: shift    ! Shift of the sublattice from
!!$!                                            the origin
!!$                                            
!***!!$      integer(4), allocatable :: iio(:,:,:) ! The symmetry operations
!!$!                                            matrices. 
!!$!                                            Dimension >= nsymt
!!$
!***!!$      real(8) :: vt                    ! Volume of the tetrahedra
!!$      
!***!!$      real(8), dimension(3,3) :: gbas  ! Basis vectors of the
!!$!                                        reciprocal lattice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine genkqtet
