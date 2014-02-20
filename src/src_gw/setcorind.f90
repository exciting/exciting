!BOP
!
! !ROUTINE: setcorind
!
! !INTERFACE:
      subroutine setcorind
!
! !DESCRIPTION:
!
! This subroutine sets a unique index for all the core states of all atoms
!
! !USES:

      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: is
      integer(4) :: ist
      integer(4) :: icg
      integer(4) :: lc 
      integer(4) :: mc
      
! !REVISION HISTORY:
!      
! Created 23.08.05 by RGA
!
!EOP
!BOC
      if (allocated(corind)) deallocate(corind)
      allocate(corind(ncg,5))
      icg=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do ist = 1, ncore(is)
            lc = spl(ist,is)
            do mc=-lc,lc
              icg=icg+1
              corind(icg,1)=is
              corind(icg,2)=ia
              corind(icg,3)=ist
              corind(icg,4)=lc
              corind(icg,5)=mc
            enddo
          enddo
        enddo
      enddo  
      ncg=icg
      return
      end subroutine setcorind
!EOC            
              
 
      
