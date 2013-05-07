!BOP
!
! !ROUTINE: calcevecplot
!
! !INTERFACE:
      subroutine calcevecplot(ik,ib,ia1,is1,ia2,is2,evec)

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the product
! of two eigenvectors in the line joining the two given atoms 
! for ploting. 
!
! !USES:

      use modinput
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik  ! The k-point for which the (L)APW+lo
!                                     function is ploted
      
      integer(4), intent(in) :: ib  ! The band index of the function

      integer(4), intent(in) :: ia1
      integer(4), intent(in) :: ia2
      integer(4), intent(in) :: is1
      integer(4), intent(in) :: is2

! !OUTPUT PARAMETERS:
      
      complex(8), intent(out) :: evec(*)
      
! !LOCAL VARIABLES:

      integer(4) :: irc
      integer(4) :: jrc
      integer(4) :: krc

      real(8) :: rd(3)
      real(8) :: kgvecl(3)
      
      integer(4) :: i, j, ikp
      integer(4) :: igp
   
      real(8) :: rdlen
      real(8) :: ri(3)
      real(8) :: phs
      real(8) :: rr(99)
     
      real(8) :: phsat
      real(8) :: irlen
      real(8) :: rd2(3)
      real(8) :: kgvec(3)
     
      complex(8) :: pw(99)
      complex(8), allocatable :: apwalm(:,:,:,:,:)
      complex(8), allocatable :: evecmtlm(:,:)
      complex(8), allocatable :: evecmt(:)
      complex(8), allocatable :: evecfv(:,:,:)
      complex(8), allocatable :: yl(:)

! !EXTERNAL ROUTINES: 

      external match
      external ylm
      external wavefmt
      external zgemv


! !INTRINSIC ROUTINES: 
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified:  June 2006 by RGA
! Revisited: 29.04.2011 by DIN
!
!EOP
!BOC

!     Allocate the local arrays
      allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))

      allocate(yl(lmmaxapw))
      allocate(evecmtlm(lmmaxapw,nrmtmax))
      allocate(evecmt(nrcmtmax))
      allocate(evecfv(nmatmax,nstfv,nspnfv))

!     get eigenvalues
      if (.not.input%gw%skipgnd) then
        call getevecfvgw(ik,evecfv)
      else
        call getevecfv(vklnr(:,ik),vgklnr(:,:,:,ik),evecfv)
      end if

!     find the matching coefficients
      call match(ngknr(1,ik),gkcnr(:,1,ik),tpgkcnr(:,:,1,ik), &
     &  sfacgknr(:,:,1,ik),apwalm(:,:,:,:,1))
      
!     Set the cartesian coordinates of the vector
      if((is1.eq.is2).and.(ia1.eq.ia2))then
          rd(1:3)=avec(1:3,1)
      else  
          rd(1:3)=atposc(1:3,ia2,is2)-atposc(1:3,ia1,is1)
      endif
      rdlen=sqrt(rd(1)*rd(1)+rd(2)*rd(2)+rd(3)*rd(3))

!     calculate the values of the spherical harmonics for atom 1
      call ylm(rd,input%groundstate%lmaxapw,yl)

!     calculate the radial wavefunctions of atom 1
      call wavefmt(input%groundstate%lradstep,input%groundstate%lmaxapw,&
     &    is1,ia1,ngknr(1,ik),apwalm(:,:,:,:,1),evecfv(:,ib,1),lmmaxapw,evecmtlm)
!
      call zgemv('T',lmmaxapw,nrcmt(is1),zone,evecmtlm,lmmaxapw,yl,1, &
     &    zzero,evecmt,1)
!      
      do irc=1,nrcmt(is1)
         evec(irc)=evecmt(irc)
      enddo
      
!     Calculate the phase of the plane waves due to the change of origin
      irlen=rdlen-rmt(is1)-rmt(is2)
      do i=1,99
         rr(i)=rmt(is1)+dble(i)*irlen/1.0d+2
      enddo  
      pw(1:99)=zzero
      do igp = 1, ngknr(1,ik)
        kgvec(1:3)=vgkcnr(1:3,igp,1,ik)
        kgvecl(1:3)=vgklnr(1:3,igp,1,ik)
        phsat=2.0d0*pi*(kgvecl(1)*atposl(1,ia1,is1)+kgvecl(2)* &
     &        atposl(2,ia1,is1)+kgvecl(3)*atposl(3,ia1,is1))
        do i=1,99
          ri(1:3)=rd(1:3)*rr(i)/rdlen
          phs=phsat+(kgvec(1)*ri(1)+kgvec(2)*ri(2)+kgvec(3)*ri(3))
          pw(i)=pw(i)+evecfv(igp,ib,1)/sqrt(omega)* &
     &            cmplx(dcos(phs),dsin(phs),8)
        enddo  
      enddo  
      do i=1,99
        j=i+nrcmt(is1)
        evec(j)=pw(i)
      enddo
     
!     calculate the radial wavefunctions of atom 2 (if needed)
      if ((is1.eq.is2).and.(ia1.eq.ia2)) then
          do irc=1,nrcmt(is1)
            jrc=nrcmt(is1)-irc+1
            krc=irc+99+nrcmt(is1)
            evec(krc)=evecmt(jrc)
          enddo
      else
          rd2(1:3)=-1.0d0*rd(1:3)
          call ylm(rd2,input%groundstate%lmaxapw,yl)
          call wavefmt(input%groundstate%lradstep,input%groundstate%lmaxapw,&
         &  is2,ia2,ngknr(1,ik),apwalm(:,:,:,:,1),evecfv(:,ib,1),lmmaxapw,evecmtlm)
          call zgemv('T',lmmaxapw,nrcmt(is2),zone,evecmtlm,lmmaxapw,yl,1, &
         &  zzero,evecmt,1)
          do irc=1,nrcmt(is2)
            jrc=nrcmt(is2)-irc+1
            krc=irc+99+nrcmt(is1)
            evec(krc)=evecmt(jrc)
          enddo
      end if
      
      deallocate(apwalm)
      deallocate(yl)
      deallocate(evecmtlm)
      deallocate(evecmt)
      deallocate(evecfv)
      
      return
      end subroutine calcevecplot
!EOC      
