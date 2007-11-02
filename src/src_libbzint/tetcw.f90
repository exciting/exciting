!BOP
!
! !ROUTINE: tetcw 
!
! !INTERFACE:
       subroutine tetcw(nik,nt,nb,wt,ebd,tetc,linkt,v,efer,omeg,sigfreq,cw)
!
! !DESCRIPTION: 
! 
! This subroutine caluclates the value of the weights for convolution
! using the improved tetrahedron method. cw(ib,jb,ik) contains the weight 
! on ik point for the bands ib and the band jb at ik+q. 
!  

! !USES:
 
       use tetra_internal
       
       implicit none      

! !INPUT PARAMETERS:
 
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
       
       real(8), intent(out)    :: cw(nb,nb,nik) ! the values of
       
   
! !SYSTEM SUBROUTINES:
       
      
       external convw

! !REVISION HISTORY:
!
!   Created: 4th. March 2004 by RGA
!
!EOP
!BOC 
!
!      Asign the input variables to the corresponding variable, or pointer in the module
!
      nirkp = nik
      ntet  = nt
      nband=nb
      vt = v
      tetcorn => tetc(1:4,1:ntet)
      tetln => linkt(1:ntet)
      eband   => ebd(1:nb,1:nik)
      tetweig => wt(1:ntet)
!
!     Calculate the weights.
!
      call convw(efer,omeg,sigfreq,cw)

      end subroutine tetcw
      
!EOC
