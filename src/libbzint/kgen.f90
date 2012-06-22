
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
! 
! !ROUTINE: kgen
!
! !INTERFACE:
      subroutine kgen(bvec,nsymt,symop,ndiv,nshift,divsh,nik,klist,idiv,&
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
      
      use kgen_internals

      implicit none

! !INPUT PARAMETERS:

      
      integer(4), intent(in) :: nsymt              ! Number of symmetry 
!                                                    operations

      integer(4), dimension(3), intent(in) :: ndiv ! Number of divisions
!                                                    of the BZ in each
!                                                    direction

      integer(4), intent(in), target :: symop(1:3,1:3,*)   ! Symmetry 
!                                                            operations 
!                                                            in internal
!                                                            coordinates 
      
      integer(4), dimension(3), intent(in) :: nshift  ! Shift vector for the
!                                                    k-points mesh
      integer(4), intent(in) :: divsh ! common divisor of the shift vector
      
      real(8), intent(in) :: bvec(3,3)
      
! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: idiv       ! Minimum common divisor or
!                                             the integer sublattice 
!                                             coordinates of the
!                                             irreducible k-points

      integer(4), intent(out) :: nik     ! Number of irreducible
!                                             k-points

      integer(4), intent(out) :: ntet       ! number of irreducible
!                                             tetrahedra

      integer(4), intent(out) :: mnd        ! an integer to tell how to 
!                                             devide the cube cell into
!                                             6 tetrahedron

      integer(4), intent(out) :: wk(*)      ! The weight of 
!                                                          each k-point
      integer(4), intent(out) :: ikpredin(*)
      integer(4), intent(out) :: klist(3,*) ! The list of 
!                                                          irreducible 
!                                                          k-points

      integer(4), intent(out) :: tetk(4,*)  ! id. number of
!                                                          the k-points
!                                                          corresponding 
!                                                          to the nodes  
!                                                          of each 
!                                                          tetrahedron.

      integer(4), intent(out) :: wtet(*)    ! index of the
!                                                          tetrahedron 
!                                                          linked to the 
!                                                          one in the 
!                                                          index by the 
!                                                          vector q

!      integer(4), intent(out) :: ikpid1(*)
!      integer(4), intent(out) :: redk(*)    ! id. number of 
!                                                          the irreducible
!                                                          k-point associ-
!                                                          ated to the 
!                                                          general k-point
!    
      real(8),    intent(out) :: vtet       ! volume of each tetrahedron


! !LOCAL VARIABLES:

      integer(4) :: i,j,nsr                   ! Just some counters
      integer(4) :: nkpt                    ! Total number of k-points
      integer(4), dimension(3) :: onek       ! Temporary storage of one
!                                             k-point coordinates
      integer(4), allocatable  :: ikp(:,:) ! Coordinates of the
!                                            irreducible k-points
      integer(4), allocatable  :: ikpint(:,:) ! Coordinates of the
!                                            irreducible k-points
      integer(4), allocatable  :: ktet(:,:) ! Temporary storage of tetk
      integer(4), allocatable  :: tetw(:)   ! Temporary storage of wtet
      integer(4), allocatable  :: wt(:)     ! Temporary storage of wk
      
!EOP
!BOC
!      iio=>symop(1:3,1:3,1:nsymt)
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
      if(allocated(ikpid)) deallocate(ikpid)
      allocate(ikpid(nkpt))
      if(allocated(redkp)) deallocate(redkp)
      allocate(redkp(nkpt))
      if(allocated(wt))deallocate(wt)
      allocate(wt(nkpt))
!      allocate(redk(nkpt))
!      allocate(ikpid1(nkpt))
      call redusym(nsymt,symop,divsh,nsr)
      call reduk(nsr,divsh,wt)
!
!     Allocate arrays depending on nirkp
!      
      if(allocated(ikp)) deallocate(ikp)
      allocate(ikp(3,nirkp))
      if(allocated(ikpint)) deallocate(ikpint)
      allocate(ikpint(3,nirkp))
      wk(1:nirkp)=wt(1:nirkp)
      wk(nirkp+1:nkpt)=0.0d0
      deallocate(wt)
!
!     Calculate the sublattice coordinates of the irreducible k-points
!      
!      ikpid1(1:nkpt)=ikpid(1:nkpt)
!      redk(1:nkpt)=redkp(1:nkpt)
      do i=1,nkpt
        ikpredin(i)=ikpid(redkp(i))
!        write(24,'(4i4)')i,ikpid(i),redkp(i),wk(i)
      enddo
      do i=1,nkpt
        if(ikpid(i).ne.0)then
          call coorskp(i,onek)
          do j=1,3
            ikp(j,ikpid(i))=onek(j)
          enddo
        endif
      enddo
!
!     Transform the irreducible k-points to internal coordinates
!     
      call intern(nirkp,ikp,ndiv,nshift,divsh,ikpint,idiv)
      klist(1:3,1:nirkp)=ikpint(1:3,1:nirkp)
      deallocate(ikpint)
!
!     Generate the tetrahedra
!       
      allocate(tetw(6*nkpt))
      allocate(ktet(4,6*nkpt))
      call tgen(ntet,ktet,tetw)
      wtet(1:ntet)=tetw(1:ntet)
      tetk(1:4,1:ntet)=ktet(1:4,1:ntet)
      deallocate(tetw)
      deallocate(ktet)
      deallocate(redtet)
      nik=nirkp
      vtet=vt
     
    
      mnd=mndg
     
      deallocate(ikpid)
      deallocate(redkp)
      deallocate(ikp)

      return

      end subroutine kgen
!EOC
