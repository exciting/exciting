!BOP
!
! !ROUTINE: tetcw 
!
! !INTERFACE:
       subroutine tetcw(nik,nt,nb,wt,ebd,tetc,linkt,linkq,v,efer,omeg,sigfreq,cw)
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
 
       integer, intent(in) :: nik                  ! Number of k-points
       integer, intent(in) :: nt                   ! Number of tetrahedra
       integer, intent(in) :: nb                   ! Number of bands
       integer, target, intent(in) :: wt(*)        ! weights of each tetrahedron
       real(8), target, intent(in) :: ebd(nb,nik)  ! Band energies
       integer, target, intent(in) :: tetc(4,*)    ! id. numbers of the corners of the tetrahedra
       integer, target, intent(in) :: linkt(*)     ! index to tell with which tetrahedron in k-q mesh is the
                                                   ! tetrahedron in k mesh linked by this q vector
       integer, target, intent(in) :: linkq(*)     ! index to tell with which kpoint corresponds to k-q
       real(8), intent(in)    :: v                 ! the volume of the tetrahedra
       real(8), intent(in)    :: efer              ! fermi energy
       real(8), intent(in)    :: omeg              ! the fequency for which the weights are calculated, 
                                                   ! only used for sigfreq =2,3 and 4
       integer, intent(in) :: sigfreq              ! Select the kind of bulk convolution weights when it equals 1,2,3 or 4

! !OUTPUT PARAMETERS:
       real(8), intent(out)    :: cw(nb,nb,nik)    ! the values of weights 
   
! !SYSTEM SUBROUTINES:
       external convw
       external average_degen_convw

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
      klinkq => linkq(1:nik)
      eband   => ebd(1:nb,1:nik)
      tetweig => wt(1:ntet)

      if(iop_integ.eq.1) then 
        if(.not.allocated(x_gauq)) then 
          allocate(x_gauq(n_gauq),w_gauq(n_gauq))
        endif 
        call gauleg(0.d0,1.d0,x_gauq,w_gauq,n_gauq)
      endif 
!
!     Calculate the weights.
!
      call convw(efer,omeg,sigfreq,cw)

!      call average_degen_convw(cw)

      end subroutine tetcw
      
!EOC
