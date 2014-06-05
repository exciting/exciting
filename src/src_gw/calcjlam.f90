!BOP
!
! !ROUTINE: calcjlam
!
! !INTERFACE:
      subroutine calcjlam(iq,is,ias,blmax)
!
! !DESCRIPTION: 
!
!This subroutine calculates the matrix elements
! $< j_{\lambda}^{|\vec{G}+\vec{q}|}>_{aNL}$.
!
!
! !USES:

      use modmain
      use modgw
  
! !INPUT PARAMETERS:
!
      implicit none
      
      integer(4), intent(in) :: iq
      integer(4), intent(in) :: is
      integer(4), intent(in) :: ias     ! Atom for which the matrix 
!                                         elements are calculated.

      integer(4), intent(in) :: blmax

! !LOCAL VARIABLES:
!
      integer(4) :: nrp
      integer(4) :: ipw ! (Counter) runs over plane waves.
      integer(4) :: irm ! (Counter) runs over radial mixed basis
!                         functions
      integer(4) :: irp ! (Counter), runs over the radial mesh points.
      integer(4) :: l   ! Angular momentum quantum number of the 
!                         mixed atomic function irm.

      real(8) :: x      ! Argument of j_l, =|q+G|.r
      real(8) :: qglen  ! Length of q+G

      real(8) :: qg(3)       ! Coordinates of the vector q+G
      real(8) :: gvec(3)     ! Coordinates of the vector G

      real(8), allocatable :: fr(:),gr(:),cf(:,:)
      real(8), allocatable :: bessl(:,:),bestemp(:) 
                                    
! !EXTERNAL ROUTINES: 

      external fderiv
      external sbessel

! !INTRINSIC ROUTINES: 

      intrinsic dsqrt

!
! !REVISION HISTORY:
!
! Created: 17th. March 2004 by MF
! Last Modified: 30th. March 2004 by RGA
! Revisited: June 2011 by DIN
!
!EOP
!BOC
      if(allocated(jlam))deallocate(jlam)
      allocate(jlam(maxnmix,ngbarc(iq)))
      
      nrp=nrmt(is)
      allocate(fr(nrp),gr(nrp),cf(3,nrp))
      allocate(bessl(nrp,0:blmax),bestemp(0:blmax))

      do ipw=1,ngbarc(iq)
        gvec(1:3)=vgc(1:3,igqigb(ipw,iq))
        qg(1:3)=vqc(:,iq)+gvec(1:3)
        qglen=dsqrt(qg(1)*qg(1)+qg(2)*qg(2)+qg(3)*qg(3))
!
!       Calculate the spherical bessel function at each mesh point
!
        do irp = 1, nrp
          x=spr(irp,is)*qglen
          call sbessel(blmax,x,bestemp)
          bessl(irp,0:blmax) = bestemp(0:blmax)
        enddo ! irp
        do irm=1,nmix(ias)
          l=bigl(ias,irm)
          do irp = 1, nrp
            fr(irp) = bessl(irp,l)*umix(ias,irm,irp)*spr(irp,is)
          enddo ! irp
!         Integrate the wavefunctions:
          call fderiv(-1,nrp,spr(1:nrp,is),fr,gr,cf)
          jlam(irm,ipw) = gr(nrp)
        enddo ! irm
      enddo ! ipw
      
      deallocate(fr,gr,cf)
      deallocate(bessl,bestemp)
      
      end subroutine calcjlam
!EOC
