
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: tetcorecw 
!
! !INTERFACE:
       subroutine tetcorecw(nik,nt,nb,nc,ebd,ec,tetc,linkt,v,efer,omeg,sigfreq,cw)
!     
! !DESCRIPTION: 
! 
! This subroutine calculates the values of the weights for convolutions for
! the case when the occupied state is a core state and the unoccupied state is
! a conduction state using the improved tetrahedron method


! !USES:
 
       use tetra_internal
       
       implicit none      

! !INPUT PARAMETERS:

       integer(4), intent(in) :: nik        ! Number of irreducible 
!                                             k-points
       
       integer(4), intent(in) :: nt         ! Number of tetrahedra
       
       integer(4), intent(in) :: nb         ! Number of bands
       
       integer(4), intent(in) :: nc         ! Number of core states involved

       real(8), target, intent(in) :: ebd(nb,nik)   ! Band energies

       real(8), target, intent(in) :: ec(*)     ! Core energies
       
       integer(4), target, intent(in) :: tetc(4,*) ! id. numbers of the
!                                                    corners of the 
!                                                    tetrahedra
   
       integer(4), target, intent(in) :: linkt(*)  ! the index which links the 
!                                                   tetrahedron in the k mesh to
!                                                   the tetrahedron in the k-q
!                                                   mesh
       
       real(8), intent(in)    :: v        ! the volume of the tetrahedra
 
       real(8), intent(in)    :: efer     ! fermi energy
       
       real(8), intent(in)    :: omeg     ! the fequency for which the
!                                           weights are calculated, 
!                                           only used ofr sigfreq =2,3

       integer(4), intent(in) :: sigfreq  ! Select the kind of convolution
!                                           weights

! !OUTPUT PARAMETERS:
       
       real(8), intent(out)    :: cw(nik,nb,nc) ! the values of
!                                                    convolution weights
       

! !EXTERNAL ROUTINES:
             
       external convcorew

! !REVISION HISTORY:
!
!   Created: 17th Jan, 2005 by XZL
!
!EOP
!BOC 
   
      nirkp = nik
      ntet  = nt
      nband=nb
      ncore=nc
      vt = v
      tetcorn => tetc(1:4,1:ntet)
      tetln => linkt(1:ntet)
      ecore => ec(1:ncore)
      eband   => ebd(1:nband,1:nik)
      call convcorew(efer,omeg,sigfreq,cw)

      end subroutine tetcorecw
      
!EOC
