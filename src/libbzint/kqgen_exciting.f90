
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
! 
! !ROUTINE: kqgen_exciting
!
! !INTERFACE:
      subroutine kqgen_exciting(bvec,ndiv,nshift,divsh,nkpt,iqnr,klist,qlist, &
           idivk,idivq,ntet,tetk,wtet,linkt,vtet)
!-----------------------------------------------------------------------------
! this interface is changed. Compared with the original version
! 
! kqid(nkpt,nkpt) is added to related the ID number of k point and k-q point
! 
! wkkq(nkpt,nkpt,nkpt) is deleted for fear of overflowing in high ndiv case
! 
! by XZL on July 21th
!--------------------------------------------------------------------
      
! !DESCRIPTION:
!
! This subroutine generates the set of irreducible k- and q-points for
!double convolutions of the form:

! \begin{equation}
! \bar{X}(\vec{q})=\int\limits_{BZ}\langle\Phi_n(\vec{k})|X(\vec{q})|%
! \Phi_{n'}(\vec{k}-\vec{q})\rangle \Theta(\epsilon_F-\epsilon_n(\vec{k}))%
! \Theta(\epsilon_{n'}(\vec{k}-\vec{q})-\epsilon_F) d^3\vec{k}
! \end{equation}
!
! and 
!
! \begin{equation}
! S(\vec{k})=\int\limits_{BZ}\langle\Phi_n(\vec{k})|\bar{X}(\vec{q})|%
! \Phi_{n'}(\vec{k}-\vec{q})\rangle  \Theta(\epsilon_F-\epsilon_n(\vec{k}))%
! \Theta(\epsilon_{n'}(\vec{k}-\vec{q})-\epsilon_F) d^3\vec{q}
! \end{equation}
!
! using the tetrahedron method. The k and k'(k-q) grids are identical. The
! tetrahedra can not be reduced by symmetry, and, for a given q, each tetrahedra with vertix
! k is linked to a different one by q. The q grid is similar to the k one
! except that it can not contain a shift. The geometrical weights of the q
! points depend on k and k'.
!
! Since the dimension of most of the output parameters are determined
! within the subroutine, these dummy array arguments are allocatable, requiring an
! explicit interface which is included in the module \verb"bz_interfaces"
! that has to be used by the calling program.
    
! !USES:
     
      use kgen_internals
!<sag>
      use control, only: tetradbglv,kplusq
!</sag>
      
      implicit none

! !INPUT PARAMETERS:
     
      integer(4), dimension(3), intent(in) :: ndiv ! Number of divisions
!                                                    of the BZ in each 
!                                                    direction
      integer(4), dimension(3), intent(in) :: nshift ! Shift vector for 
!                                                      the k-points mesh
      integer(4), intent(in) :: divsh
      integer(4), intent(in) :: nkpt
      ! q-point index in the non-reduced set
      integer, intent(in) :: iqnr
      real(8), intent(in) :: bvec(3,3)

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: klist(3,nkpt) ! The list of
!                                                          irreducible
!                                                          k-points

      integer(4), intent(out) :: qlist(3,nkpt) ! The list of
!                                                          irreducible
!                                                          q-points

      integer(4), intent(out) :: idivk  ! Minimum common divisor or
!                                        the integer sublattice 
!                                        coordinates of the
!                                        irreducible k-points

      integer(4), intent(out) :: idivq  ! Minimum common divisor or
!                                        the integer sublattice 
!                                        coordinates of the
!                                        irreducible q-points

      integer(4), intent(out) :: ntet       ! Number of tetrahedra

      integer(4), intent(out) :: tetk(4,6*nkpt)  ! id. number of
!                                                          the k-points
!                                                          corresponding
!                                                          to the nodes 
!                                                          of each 
!                                                          tetrahedron.
      integer(4), intent(out) :: wtet(*)    ! index of the

      integer(4), intent(out) :: linkt(6*nkpt) ! index of the
!                                                          tetrahedron 
!                                                          linked to the 
!                                                          one in the
!                                                          index by the 
!                                                          vector q
      
      

      real(8),    intent(out) :: vtet    ! volume of each tetrahedron
!
! !LOCAL VARIABLES:

      integer(4) :: i,j,itet                 ! Just some counters
      integer(4), dimension(3) :: onek       ! Temporary storage of one
!                                             k-point coordinates
      integer(4), dimension(3) :: q          ! Vector q for which the set
!                                             of(k,k') points are 
!                                             calculated
      integer(4) :: kpt(3,nkpt) ! Coordinates of the k-points
      integer(4) :: link1(6*nkpt)  ! Temporary storage for the
!                                             index of the linked 
!                                             tetrahedra corresponding to
!                                             one q.


      integer(4), external :: idkp     

! ! SYSTEM ROUTINES:
    
      intrinsic modulo

!EOP
!     
!      interface
!        subroutine asockkp(wqk)
!          integer(4), intent(out) :: wqk(:,:,:) 
!        end subroutine
!      end interface
!      
   
!BOC
!      iio=>symop(1:3,1:3,1:nsymt)
     div=ndiv
      gbas=bvec
      if (tetradbglv > 0) then     ! sag
         write(*,*)'div =',div
         write(*,*)'nkpt =',nkpt
      end if                       ! sag
      linkt=0
!
!     Set the total number of k-points
!
      if(nkpt.lt.ndiv(1)*ndiv(2)*ndiv(3))goto 999
!
!     Allocate arrays depending on nkpt
!      
!      allocate(link1(6*nkpt))
!      allocate(kpt(3,nkpt))
!
!     Generate the q-grid.
!      
!
!     Set the sublattice coordinates of all the k-points (kpt) and only
!     the irreducible ones (ikp)
!
      do i=1,nkpt
        call coorskp(i,onek)
        do j=1,3
          kpt(j,i)=onek(j)
        enddo
      enddo
!
!     Change  the coordinates of the k-points to internal reciprocal
!     lattice ones
!     
      call intern(nkpt,kpt,ndiv,nshift,divsh,klist,idivk)
      
!
!     Change  the coordinates of the q-points to internal reciprocal
!     lattice ones
!     
      if((nshift(1).eq.0).and.(nshift(2).eq.0).and.&
     &   (nshift(3).eq.0))then
        qlist(1:3,1:nkpt)=klist(1:3,1:nkpt)
        idivq=idivk
      else
        call intern(nkpt,kpt,ndiv,(/0,0,0/),1,qlist,idivq)
      endif
!
!     Calculate the k-dependent geometrical weights of the (q,k') pairs
!
!       call asockkp(wkkq)
      wtet(1:ntet)=1
 
!!$      do i=1,nkpt
        q(1:3)=kpt(1:3,iqnr) !+-+-
        call tgenq(q,ntet,tetk,link1)
        do itet=1,ntet
          linkt(itet)=link1(itet)       !+-+-
        enddo
!!$      enddo
!!$!----------------------------------------------------------------------
!!$! to calculate kqid(k,q)=ID(k-q)  added by XZL on July 21th
!!$!----------------------------------------------------------------------    
!!$      do i=1,nkpt
!!$       call coorskp(i,q)
!!$       do j=1,nkpt
!!$         call coorskp(j,k)
!!$         do ik=1,3
!!$            kq(ik)=modulo(k(ik)-q(ik),ndiv(ik))
!!$         enddo
!!$         kqid(j,i)=idkp(kq)
!!$       enddo
!!$      enddo  
!!$!      write(*,*)allocated(link1),sizeof(link1)
!!$!      write(*,*)allocated(kpt),sizeof(kpt)
!!$!----------------------------------------------------------------------
!!$! ended
!!$!----------------------------------------------------------------------
      vtet=vt
      
!     deallocate(link1)
!      deallocate(kpt)
      return
      
  999 write(*,*)'nkpt too small, should be',ndiv(1)*ndiv(2)*ndiv(3)
      stop 'error in kqgen'     

    end subroutine kqgen_exciting
!EOC
