
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: tetcw_1k 
!
! !INTERFACE:
subroutine tetcw_1k(iklib,nik,nt,nb,wt,ebd,tetc,linkt,v,efer,omeg, &
     sigfreq,cw)
!
! !DESCRIPTION: 
! 
! This subroutine caluclates the value of the weights for convolution
! using the improved tetrahedron method. cw(ib,jb) contains the weight 
! on ik point for the bands ib and the band jb at ik+q. 
!  
! !USES:
       use tetra_internal
       !<sag>
       use control, only: pointerhandling
       !</sag>
       implicit none      
! !INPUT PARAMETERS:
       integer, intent(in) :: iklib ! library k-point number
       integer(4), intent(in) :: nik        ! Number of k-points
       integer(4), intent(in) :: nt         ! Number of tetrahedra
       integer(4), intent(in) :: nb         ! Number of bands
       integer(4), target, intent(in) :: wt(*)
       real(8), target, intent(in) :: ebd(nb,nik)   ! Band energies
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
! !OUTPUT PARAMETERS:
       real(8), intent(out)    :: cw(nb,nb) ! the values of
  !<sag>
  ! local variables
  real(8), target, allocatable :: target_ebd(:,:)
  integer(4), target, allocatable:: target_tetc(:,:)
  integer(4), target, allocatable :: target_linkt(:)
  integer(4), target, allocatable :: target_wt(:)
  !</sag>
! !SYSTEM SUBROUTINES:
  external convw
! !REVISION HISTORY:
!
!   Created: August 2008 by S. Sagmeister
!   Based upon the routine {\tt tetcw.f90}
!
!EOP
!BOC 
!
  ! Asign the input variables to the corresponding variable, or pointer in the
  ! module
  nirkp = nik
  ntet  = nt
  nband=nb
  vt = v
  !<sag>
  if (pointerhandling.eq.0) then
     ! default treatment
     !</sag>
     tetcorn => tetc(1:4,1:ntet)
     tetln => linkt(1:ntet)
     eband   => ebd(1:nb,1:nik)
     tetweig => wt(1:ntet)
     !<sag>
  else if (pointerhandling.eq.1) then
     ! additional target arrays to get around with core dumps in
     ! combination with the Portland compiler
     allocate(target_tetc(1:4,1:ntet))
     allocate(target_linkt(1:ntet))
     allocate(target_ebd(1:nb,1:nik))
     allocate(target_wt(1:ntet))
     ! store copy of input parameters locally for this routine
     target_tetc(:,:)=tetc(1:4,1:ntet)
     target_linkt(:)=linkt(1:ntet)
     target_ebd(:,:)=ebd(1:nb,1:nik)
     target_wt(:)=wt(1:ntet)
     ! assign pointers to local objects
     tetcorn=>target_tetc(:,:)
     tetln=>target_linkt(:)
     eband=>target_ebd(:,:)
     tetweig=>target_wt(:)
  end if
  !</sag>
  !
  !     Calculate the weights.
  !
  call convw_1k(iklib,efer,omeg,sigfreq,cw)
  !<sag>
  if (pointerhandling.eq.1) then
     deallocate(target_tetc,target_linkt,target_ebd,target_wt)
  end if
  !</sag>  
end subroutine tetcw_1k
!EOC
