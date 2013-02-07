!BOP
!
! !ROUTINE: calcmpwmix
!
! !INTERFACE:
      subroutine calcmpwmix(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the matrix elements between mixed basis
! functions and plane waves.
!
! !USES:

      use modmain
      use modgw 
            
! !INPUT PARAMETERS:      
      implicit none
      
      integer(4), intent(in) :: iq
!
!
! !LOCAL VARIABLES:

      integer(4) :: ipin
      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: imix ! (Counter): runs over local mixed basis
!                                     functions
      integer(4) :: irm  ! (Counter): runs over radial mixed basis
!                                     functions      
      integer(4) :: ipw  ! (Counter): runs over plane waves
      integer(4) :: is
      integer(4) :: l1   ! angular momentum quantum number of the mixed
!                          basis function
      integer(4) :: m1   ! z component of the angular momentum of the mbf
      integer(4) :: lm1
      integer(4) :: ylsize 

      real(8) :: gvecl(3)
      real(8) :: gpr     ! Scalar product G.pos(atom)
      real(8) :: prefac
      
      real(8), dimension(3) :: qvec ! the q-vector which is zero
      real(8), dimension(3) :: qg1 ! the G-vector half rotated
      
      complex(8) :: expg
      complex(8) :: y((maxbigl+1)*(maxbigl+1))! spherical harmonics

      complex(8), allocatable :: sph(:,:)

 
! !EXTERNAL ROUTINES: 

      external calcjlam
      external ylm

! !INTRINSIC ROUTINES: 

      intrinsic sqrt
      intrinsic exp
      intrinsic cos
      intrinsic sin

!
! !REVISION HISTORY:
! 
! Created: 23. Nov. 2004 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC
      if(allocated(mpwmix))deallocate(mpwmix)
      allocate(mpwmix(matsiz,ngbarc(iq)))
      mpwmix(:,:)=zzero
      
      prefac=cmplx(4.0d0*pi*sqrt(vi),0.0d0,8)
      ylsize=(maxbigl+1)*(maxbigl+1)
      allocate(sph(ylsize,ngbarc(iq)))
      
      qvec(1:3)=vqc(1:3,iq) 
      
      if (Gamma) then
        mpwmix(1:matsiz,1)=wi0(1:matsiz)
        ipin=2
      else
        ipin=1
      endif    
!
!     Loop over inequivalent atoms:
!
      imix = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias=idxas(ia,is)
          call calcjlam(iq,is,ias,mbl(ias))
!
!         Calculate Y_lm(q+G) for all G and phases
!
          do ipw = 1, ngbarc(iq)
            qg1(1:3)=qvec(1:3)+vgc(1:3,igqigb(ipw,iq))
            call ylm(qg1,maxbigl,y)
            sph(1:ylsize,ipw)=y(1:ylsize)
          enddo
!
!         Loop over mixed functions:
!
          do irm = 1,nmix(ias)
            l1=bigl(ias,irm)
            do m1 = -l1,l1
              
              imix = imix+1
              lm1=l1*(l1+1)+m1+1

              do ipw = ipin, ngbarc(iq)
                gvecl(1:3)=dble(ivg(1:3,igqigb(ipw,iq)))
                gpr=gvecl(1)*atposl(1,ia,is)+ &
               &    gvecl(2)*atposl(2,ia,is)+ &
               &    gvecl(3)*atposl(3,ia,is)
                expg=cmplx(cos(2.0d0*pi*gpr),sin(2.0d0*pi*gpr),8)
                
                mpwmix(imix,ipw)=cmplx(prefac*jlam(irm,ipw),0.0d0,8)* &
               &    expg*(zi**l1)*conjg(sph(lm1,ipw))

              enddo ! ipw
            enddo ! m1  
          enddo ! irm
        
        enddo ! ia
      enddo ! is 
      
      do imix=locmatsiz+1,matsiz
        do ipw=1,ngbarc(iq)
          mpwmix(imix,ipw)=mpwipw(imix-locmatsiz,ipw)
        enddo
      enddo  
      deallocate(sph)

      return
      end subroutine calcmpwmix
!EOC
