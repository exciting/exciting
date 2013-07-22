!BOP
!
! !ROUTINE: setlocmixind
!
! !INTERFACE:
      subroutine setlocmixind
!
! !DESCRIPTION:
!
! This subroutine sets an array that stores the general index of the mixed
! function for a given mixed function of a given atom
!
! !USES:

      use modmain
      use modgw


! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias     ! (Counter) Runs over inequivalent atoms
      integer(4) :: im
      integer(4) :: imix     ! (Counter) runs over all core states
      integer(4) :: irm
      integer(4) :: is
      integer(4) :: l      ! Angular momentum of the core state
      integer(4) :: m
      
! !REVISION HISTORY:
!      
! Created 23.08.05 by RGA
!
!EOP
!BOC
      if(debug)then
        write(701,*)'------------------------------------------------------'
        write(701,*)'     Indexes of MT-sphere mixed basis functions'
        write(701,*)'------------------------------------------------------'
        write(701,*)'   chi_i=v_(aNL)Y_(LM)   (a = atom)'
        write(701,101)
      end if
      
      if (allocated(locmixind)) deallocate(locmixind)
      allocate(locmixind(natmtot,lmixmax))
      locmixind(:,:)=0
      im=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          imix=0
          do irm = 1, nmix(ias)
            l = bigl(ias,irm)
            do m=-l,l
              imix=imix+1
              im=im+1
              locmixind(ias,imix)=im
              if(debug)write(701,102)im,ias,irm,l,m
            enddo ! m
          enddo ! irm
        enddo ! ia
      enddo ! is
      
      if(debug)write(701,*)'------------------------------------------------------'
  101 format(5x,'i',5x,'a',5x,'N',5x,'L',5x,'M')
  102 format(5i6) 
      
      return
      end subroutine setlocmixind
!EOC            
              
 
      
