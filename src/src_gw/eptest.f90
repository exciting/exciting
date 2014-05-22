!BOP
!
! !ROUTINE: eptest
!
! !INTERFACE:
      subroutine eptest(ik,jk,iq)

! !DESCRIPTION:
!
!This subroutine calculates the matrices $M^i_{cm}(k,k')$ and
!$M^i_{nm}(k,k')$ for the given k and k' pair. 
! It is only used for testing purposes... for the calculation of the
!matrices for all pears see \verb"expand_prods"
!
! !USES:

       use modmain
       use modgw

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik     ! index of the k-point.
      integer(4), intent(in) :: jk     ! index of the k'-points.
      integer(4), intent(in) :: iq     ! index of the q-point

! !LOCAL VARIABLES:

! !EXTERNAL ROUTINES: 

      external calcminm
      external expand_evec

! !REVISION HISTORY:
!
! Created: May 2006 by RGA
! Revisited 29.04.2011 by DIN
!
!EOP
!BOC
!
      matsiz=locmatsiz+ngq(iq)
      write(*,101)iq,locmatsiz,ngq(iq),matsiz

      allocate(minmmat(matsiz,nstfv,nstfv))
      
      allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
      allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
      allocate(eveck(nmatmax,nstfv,nspnfv))
      allocate(eveckp(nmatmax,nstfv,nspnfv))
!
!     get the eigenvectors from file
      call getevecfvgw(jk,eveck)
      eveckp=conjg(eveck)
      call getevecfvgw(ik,eveck)
!
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')
!          
!     calculate the matrix elements $M^i_{nm}(\vec{k},\vec{q})$:
!
      call calcminm(ik,iq,0)
!
      deallocate(eveck)
      deallocate(eveckp)
      deallocate(eveckalm)
      deallocate(eveckpalm)
      deallocate(bradketa)
      deallocate(bradketlo)
      deallocate(bradketc)
!
  101 format(10x,'Data for q-point nr.:',i4,//,10x,'Mixed basis:',/,10x, &
     &           'Number of atomic basis functions:       ',i4,/,10x,   &
     &           'Number of interstitial basis functions: ',i4,/,10x,   &
     &           'Total number of basis functions:        ',i4,/)
!    

      return      
      end subroutine eptest
!EOC

