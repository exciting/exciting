!BOP
! 
! !ROUTINE: kgen
!
! !INTERFACE:
      subroutine kgen(bvec,nsym,symop,ndiv,nshift,divsh,nik,klist,idiv,&
     &                ikpredin,wk,ntet,tetk,wtet,vtet,mnd)
!
! !DESCRIPTION:
!
! This subroutine generates the set of irreducible k-points for tetrahedron
! integration in the Brillouin Zone. The list of k-points are given in
! integer internal coordinates of the reciprocal lattice. The symmetry
! operations must also be given to the subroutine in internal coordinates.
! If the calling program is using cartesian coordinates (as is the case in
! the WIEN2k package, when \verb"ortho" = .\verb"true".) the subroutine
! \verb"sym2int" (also contained in this library) has to be called before
!to transform the symmetry operations to internal coordinates and the
!subroutine \verb"cartezian" has to be called afterwards to transform the
!k-points coordinates back to cartesian ones.
!
! Since the dimension of most of the output parameters are determined
! within the subroutine, these dummy array arguments are allocatable, requiring an
! explicit interface to be written in the calling program.
!

! !USES:
      
      use tetra_internal, only: redtet,mndg,fout
      use kgen_internals

      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: nsym                  ! Number of symmetry operations
      integer, dimension(3), intent(in) :: ndiv    ! Number of k-points in each direction
      integer, intent(in) :: symop(1:3,1:3,nsym)   ! Symmetry operations in internal coordinates 
      integer, dimension(3), intent(in) :: nshift  ! Shift vector for the k-points mesh
      integer, intent(in) :: divsh                 ! common divisor of the shift vector
      real(8), intent(in) :: bvec(3,3)             ! recipirical lattice vectors 
      
! !OUTPUT PARAMETERS:

      integer, intent(out) :: idiv   ! Minimum common divisor of k-points
      integer, intent(out) :: nik    ! Number of irreducible  k-points
      integer, intent(out) :: ntet   ! number of irreducible tetrahedra
      integer, intent(out) :: mnd    ! an integer to tell how to devide the cube cell into 6 tetrahedron
      integer, intent(out) :: wk(*)  ! The weight of each k-point
      integer, intent(out) :: ikpredin(*) 
      integer, intent(out) :: klist(3,*) ! The list of irreducible k-points
      integer, intent(out) :: tetk(4,*)  ! index of the k-points corresponding to the nodes of each tetrahedron.
      integer, intent(out) :: wtet(*)    ! index of the tetrahedron linked to the one in the index by the vector q
      real(8),    intent(out) :: vtet       ! volume of each tetrahedron

! !LOCAL VARIABLES:
      integer:: ierr
      integer :: i,j,nsr                   ! Just some counters
      integer :: nkpt                      ! Total number of k-points
      integer, dimension(3) :: onek        ! Temporary storage of one k-point coordinates
      integer, allocatable  :: ikp(:,:)    ! Coordinates of the irreducible k-points
      integer, allocatable  :: ikpint(:,:) ! Coordinates of the irreducible k-points in internal coordinates
      integer, allocatable  :: ktet(:,:)   ! Temporary storage of tetk
      integer, allocatable  :: tetw(:)     ! Temporary storage of wtet
      integer, allocatable  :: wt(:)       ! Temporary storage of wk
      
!EOP
!BOC
      div=ndiv
      shift=nshift
      gbas=bvec
!
!     Set the total number of k-points
!
      nkpt=ndiv(1)*ndiv(2)*ndiv(3)
!
!     Allocate arrays depending on nkpt
!
      allocate(ikpid(nkpt),redkp(nkpt),wt(nkpt),stat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "ERROR in kgen: fail to allocate ikpid,redkp,wt"
        stop "ERROR in kgen"
      endif 
 
      call redusym(nsym,symop,divsh,nsr)
      call reduk(nsr,divsh,wt)
!
!     Allocate arrays depending on nirkp
!
      allocate(ikp(3,nirkp),ikpint(3,nirkp),stat=ierr)
      if(ierr.ne.0) then
        write(6,*) "ERROR in kgen: fail to allocate ikp,ikpint"
        stop "ERROR in kgen"
      endif

      wk(1:nirkp)=wt(1:nirkp)
      wk(nirkp+1:nkpt)=0.0d0
!
!     Calculate the sublattice coordinates of the irreducible k-points
!
      do i=1,nkpt
        ikpredin(i)=ikpid(redkp(i))
      enddo

      do i=1,nkpt
        if(ikpid(i).ne.0)then
          call coorskp(i,onek)
          ikp(1:3,ikpid(i))=onek(1:3)
        endif
      enddo
!
!     Transform the irreducible k-points to internal coordinates
!
      call intern(nirkp,ikp,ndiv,nshift,divsh,ikpint,idiv)
      klist(1:3,1:nirkp)=ikpint(1:3,1:nirkp)
!
!     Generate the tetrahedra
!
      allocate(tetw(6*nkpt))
      allocate(ktet(4,6*nkpt))
      call tgen(ntet,ktet,tetw)
      wtet(1:ntet)=tetw(1:ntet)
      tetk(1:4,1:ntet)=ktet(1:4,1:ntet)
      nik=nirkp
      vtet=vt
      mnd=mndg
     
      deallocate(ikpid,ikpint,redkp,ikp,wt,tetw,ktet,redtet)

      return

      end subroutine kgen
!EOC
