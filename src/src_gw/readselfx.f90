!BOP
!
! !ROUTINE: readselfx
!
! !INTERFACE:
      subroutine readselfx

! !DESCRIPTION:
! 
! This subroutine writes the self-energy to file
!
! !USES:
!
      use modmain
      use modgw     
       
! !LOCAL VARIABLES:
      
      integer(4) :: ie,i     !(Counter) Runs over bands
      integer(4) :: ikp, ik  !(Counter) Runs over k-points
      real(8)    :: kvec(3)
      real(8)    :: sxr, sxi
      integer(4) :: ib, nb
      
! !REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
      open(92,file='SELFX.OUT',action='READ',form='FORMATTED')
      read(92,*)
      do ikp = 1, nkpt
        read(92,*) kvec, ik, ib, nb
        if(ik.ne.ikp)then
          write(6,*)'ERROR(readselfx): Inconsistent input parameters'
          write(6,*)'ik=',ik,'ikp=',ikp
          stop
        end if
        if ((ib.ne.ibgw).or.(nb.ne.nbgw)) then
          write(6,*)'ERROR(readselfx): Inconsistent input parameters'
          write(6,*)'ib=',ib,'nb=',nb
          write(6,*)'ibgw=',ibgw,'nbgw=',nbgw
          stop
        end if
        do ie = ibgw, nbgw
          read(92,*) i, sxr, sxi
          selfex(ie,ikp)=cmplx(sxr,sxi,8)/hev
        enddo
      enddo    
      close(92)
      
      return
      end subroutine
!EOC      
