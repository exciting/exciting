!BOP
!
! !ROUTINE: testacont
!
! !INTERFACE:
       subroutine testacont

       use modmain
       use modgw
       
! !DESCRIPTION:
!
! This subroutine perform the analytic continuation of the selfenergy to
! real frequencies and calculates the quasi-particle energies. The
! selfenergy and exchange correlation potential are read from file, thus, 
! a previous run of the GW cycle is needed.        
!

! !EXTERNAL ROUTINES: 

      external readvxcnn
      external readselfx
      external readselfe
      external calceqp

!      
! !REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC

!     Allocate the arrays
      allocate(vxcnn(ibgw:nbgw,nkpt))
      allocate(selfex(ibgw:nbgw,nkpt))
      allocate(selfec(ibgw:nbgw,nkpt,nomeg))
!
!     Read the exchange correlation matrix elements
!
      call readvxcnn
!
!     Read the exchange and correlation selfenergy matrix elements
!            
      call readselfx
      call readselfc
!
!     Calculate the quasiparticle energies
!            
      call calceqp
      
      return
      
      end subroutine testacont
!EOC      

