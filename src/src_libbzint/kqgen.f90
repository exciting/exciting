!BOP
! 
! !ROUTINE: kqgen
!
! !INTERFACE:
      subroutine kqgen(bvec, ndiv, nshift, divsh,&
                     & nkpt, klist, qlist, idivk,&
                     & idivq, kqid, ntet, tetk, wtet, linkt, vtet)
!-----------------------------------------------------------------------------
! this interface is changed. Compared with the original version
!  kqid(nkpt, nkpt) is added to related the ID number of k point and k-q point
!  wkkq(nkpt, nkpt, nkpt) is deleted for fear of overflowing in high ndiv case
! 
! by XZL on July 21th
!--------------------------------------------------------------------
      
! !DESCRIPTION:
!
! This subroutine generates the set of irreducible k- and q-points for
! double convolutions of the form:

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
      
      implicit none

! !INPUT PARAMETERS:
     
      integer(4), dimension(3), intent(in) :: ndiv   ! Number of k-points in each direction
      integer(4), dimension(3), intent(in) :: nshift ! Shift vector for the k-points mesh
      integer(4), intent(in) :: divsh
      integer(4), intent(in) :: nkpt
      real(8), intent(in) :: bvec(3,3)

! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: klist(3, nkpt) ! The list of irreducible k-points
      integer(4), intent(out) :: qlist(3, nkpt) ! The list of irreducible q-points
      integer(4), intent(out) :: kqid(nkpt, nkpt) ! ID number of k-q point for k(ID) to a  specific q(ID)

      integer(4), intent(out) :: idivk  ! Minimum common divisor or the integer sublattice 
                                        ! coordinates of the irreducible k-points
      integer(4), intent(out) :: idivq  ! Minimum common divisor or the integer sublattice 
                                        ! coordinates of the irreducible q-points
      integer(4), intent(out) :: ntet   ! Number of tetrahedra
      integer(4), intent(out) :: tetk(4,6*nkpt)  ! id. number of the k-points corresponding
                                                 ! to the nodes of each tetrahedron.
      integer(4), intent(out) :: wtet(*)    ! index of the
      integer(4), intent(out) :: linkt(6*nkpt, nkpt) ! index of the tetrahedron linked to the one in the
                                                    ! index by the vector q
      real(8),    intent(out) :: vtet    ! volume of each tetrahedron
!
! !LOCAL VARIABLES:

      integer(4) :: iq, ik           ! Just some counters
      integer(4), dimension(3) :: onek    ! Temporary storage of one k-point coordinates
      integer(4), dimension(3) :: q, k, kq  ! Vector q for which the set of(k, k') points are calculated
      integer(4) :: kpt(3, nkpt)           ! Coordinates of the k-points
                                          ! tetrahedra corresponding to one q.
      integer(4), external :: idkp     

! ! SYSTEM ROUTINES:
    
      intrinsic modulo

!EOP
!BOC

      div=ndiv
      gbas=bvec
      linkt=0
!
!     Set the total number of k-points
!
      if(nkpt.lt.ndiv(1)*ndiv(2)*ndiv(3)) then 
        write(6,*)'ERROR in kggen: nkpt too small, should be', &
     &         ndiv(1)*ndiv(2)*ndiv(3)
        stop "ERROR in kggen"
      endif 

!
!     Set the sublattice coordinates of all the k-points (kpt) and only
!     the irreducible ones (ikp)
!
      do ik=1, nkpt
        call coorskp(ik, onek)
        kpt(1:3, ik) = onek(1:3) 
      enddo
!
!     Change  the coordinates of the k-points to internal reciprocal
!     lattice ones
!
      call intern(nkpt, kpt, ndiv, nshift, divsh, klist, idivk)

!
!     Change  the coordinates of the q-points to internal reciprocal lattice ones
!
      if( all(nshift .eq. 0) ) then 
        qlist(1:3,1:nkpt)=klist(1:3,1:nkpt)
        idivq=idivk
      else
        call intern(nkpt, kpt, ndiv,(/0, 0, 0/),1, qlist, idivq)
      endif
!
!     Calculate the k-dependent geometrical weights of the (q, k') pairs
!
      wtet(1:ntet)=1
 
      do iq=1, nkpt
        q(1:3)=kpt(1:3, iq)
        call tgenq(q, ntet, tetk, linkt(:, iq))
      enddo

!----------------------------------------------------------------------
! to calculate kqid(k, q)=ID(k-q)  added by XZL on July 21th
!----------------------------------------------------------------------    
      do iq=1, nkpt
        call coorskp(iq, q)
        do ik=1, nkpt
          call coorskp(ik, k)
          kq(:)= modulo(k(:)-q(:), ndiv(:))
          kqid(ik, iq)=idkp(kq)
        enddo
      enddo  
      vtet=vt
      
      return
      end subroutine kqgen
!EOC
