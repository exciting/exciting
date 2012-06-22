!BOP
!
! !ROUTINE: expand_prods
!
! !INTERFACE:
      subroutine expand_prods(ik,iq,flag)

! !DESCRIPTION:
!
! This subroutine calculates $M^i_{cm}(k,q)$, $M^i_{nc}(k,q)$, 
! and $M^i_{nm}(k,q)$ matrix elements.
!
! !USES:

      use modmain
      use modgw
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: ik, iq
      integer(4), intent(in) :: flag ! 0 - M^i_{nm}
                                     ! 1 - M^i_{nm} in the transformed 
                                     !     v-diag. basis set 

! !LOCAL VARIABLES:

      integer(4) :: jk
      real(8)    :: tstart, tend
 
! !EXTERNAL ROUTINES: 
      external calcminc
      external calcmicm
      external calcminm
      external expand_evec
      external getevecfv

! !REVISION HISTORY:
!
! Created: 17th. June 2004 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

      allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
      allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
      allocate(eveck(nmatmax,nstfv,nspnfv))
      allocate(eveckp(nmatmax,nstfv,nspnfv))

      jk=kqid(ik,iq) ! index of k-q vector

!     get the eigenvectors from file
      call getevecfv(vklnr(:,jk),vgklnr(:,:,:,jk),eveck)
      eveckp=conjg(eveck)
      call getevecfv(vklnr(:,ik),vgklnr(:,:,:,ik),eveck)

      call expand_evec(ik,'t')
      call expand_evec(jk,'c')

      if ((flag>0).and.(debug)) then
        write(fgw,*)
        write(fgw,*) 'WARNING(expand_prods): Minm matrix elements are transformed to &
     &                 the eigenvalues of the bare coulomb potential'
        write(fgw,*)
      end if
!
!     Calculate the matrix elements $M^i_{nm}(\vec{k},\vec{q})$:
!
      call calcminm(ik,iq,flag)
!          
!     Calculate the matrix elements M^i_{nm} where n is a core state
!
      if(wcore)then
        
        call calcmicm(ik,iq,flag)
        call calcminc(ik,iq,flag)

      endif

      deallocate(eveck)
      deallocate(eveckp)
      deallocate(eveckalm)
      deallocate(eveckpalm)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'EXPAND_PRODS')

      return
      end subroutine expand_prods
!EOC

