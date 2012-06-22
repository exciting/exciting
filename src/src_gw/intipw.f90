!BOP
!
! !ROUTINE: intipw
!
! !INTERFACE:
      subroutine intipw
      
! !DESCRIPTION:

!This function calculates the integral of a plane wave with wave vector
!belonging to the reciprocal Bravais lattice in the
!interstitial region by the difference between the integral over the whole
!unit cell and the Muffin Tin spheres:
!
!\begin{equation}
!\frac{1}{\Omega}\int\limits_{I}{e^{i\vec{G}\cdot\vec{r}}d^3r}=%
!\frac{1}{\Omega}\int\limits_{\Omega}{e^{i\vec{G}\cdot\vec{r}}d^3r}-%
!\frac{1}{\Omega}\sum\limits_a{\int\limits_{MT_a}{e^{i\vec{G}\cdot\vec{r}}d^3r}}
!\end{equation}
!being:
!
!\begin{equation}
!\frac{1}{\Omega}\int\limits_{\Omega}{e^{i\vec{G}\cdot\vec{r}}d^3r}=\delta_{\vec{G},0}
!\end{equation}
!
!and
!
!\begin{equation}
!\frac{1}{\Omega}\int\limits_{MT_a}{e^{i\vec{G}\cdot\vec{r}}d^3r}=\frac{V_{MT_a}}{\Omega}\left[\delta_{\vec{G},0}+%
!3(1-\delta_{\vec{G},0})\frac{\sin(|\vec{G}|R_a)-(|\vec{G}|R_a)%
!\cos(|\vec{G}|R_a)}{(|\vec{G}|R_a)^3} e^{i\vec{G}\cdot\vec{r}_a}\right]
!\end{equation}
!
! !USES:
      
      use modmain
      use modgw
      
! !INPUT PARAMETERS:

      implicit none

! !LOCAL VARIABLES:

      integer(4) :: ig
      integer(4) :: i,ia,is

      real(8) :: glen   ! Length of g
      real(8) :: gr     ! glen times MT radius
      real(8) :: intmod ! Modulus of the MT integral
      real(8) :: phase  ! the phase of the MT integral
      real(8) :: mtintr ! The integral over the MT Sphere (real part)
      real(8) :: mtinti ! The integral over the MT Sphere (imaginary part)

      real(8) :: integr
      real(8) :: integi
      
      complex(8) :: integral

! !EXTERNAL ROUTINES:      
      


! !REVISION HISTORY:
!
! Created 30. April 2004 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC
! allocate global characteristic function arrays
      if (allocated(ipwint)) deallocate(ipwint)
      allocate(ipwint(ngrtot))

! begin loop over species
      do ig = 1, ngrtot
         integr=0.0d0
         integi=0.0d0
         glen=dsqrt(vgc(1,ig)*vgc(1,ig)+ &
        &           vgc(2,ig)*vgc(2,ig)+ &
        &           vgc(3,ig)*vgc(3,ig))
         if (glen.gt.1.0d-10) then
             do is=1,nspecies
                gr=glen*rmt(is)
                intmod=3.0d0*vmt(is)*(dsin(gr)-gr*dcos(gr))/(gr*gr*gr)
                mtintr=0.0d0
                mtinti=0.0d0
                do ia=1,natoms(is)
                   phase=0.0d0
                   do i=1,3
                      phase=phase+vgc(i,ig)*atposc(i,ia,is)
                   enddo
                   mtintr=mtintr+dcos(phase)
                   mtinti=mtinti+dsin(phase)
                end do ! ia
                integr=integr-intmod*mtintr
                integi=integi-intmod*mtinti
             end do ! is
             integral=cmplx(integr,integi,8)  
         else  
             integral=cmplx(1.0d0,0.0d0,8)
             do is=1,nspecies
                integral=integral-cmplx(vmt(is)*dble(natoms(is)),0.0d0,8)
             end do
         endif  
                        
         ipwint(ig)=integral
         
      end do ! ig
      
! Test     
!     sum = 0.0d0
!     do ig=1,ngq(1)
!        sum=sum+ipwint(ig)*conjg(ipwint(ig))
!     enddo      
!     write(*,*)'sumipwint =', sum
!     write(*,*) ivgig(0,0,0),ipwint(ivgig(0,0,0))
!     sum = 0.0d0
!     do ig=1,ngq(1)
!        sum=sum+cfunig(ig)*conjg(cfunig(ig))
!     enddo      
!     write(*,*)'sumcfunig =', sum
!     write(*,*) ivgig(0,0,0),cfunig(ivgig(0,0,0))
!     stop
      
      end subroutine intipw
      
!EOC      
              
            

