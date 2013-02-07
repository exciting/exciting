!
!BOP
! !ROUTINE: getdlmm
! !INTERFACE:
!
!
complex(8) function getdlmm(rot,l,m1,m2)
! !INPUT/OUTPUT PARAMETERS:
!   rot   : rotation matrix (in,real(3,3))
!   l     : angular momentum (in,integer)
!   m1    : 
!   m2    : 
! !DESCRIPTION:
!   Calculates the rotational matrix $D^l_{mm'}$ for the given rotation matrix $R$.
!   This is done by first the computing the Euler angles $(\alpha,\beta,\gamma)$
!   of $R$ (see routine {\tt euler}) and then generating the rotation
!   matrix for spherical harmonics, $D^l_{mm'}(\alpha,\beta,\gamma)$, with which
!   $$ Y_{lm}(\theta',\phi')=\sum_{m'}D^l_{mm'}(\alpha,\beta,\gamma)Y_{lm'}
!    (\theta,\phi), $$
!   where $(\theta',\phi')$ are the angles $(\theta,\phi)$ rotated by $R$. The
!   matrix $D$ is given explicitly by
!   \begin{align*}
!    D^l_{mm'}(\alpha,\beta,\gamma)=&\sum_i\frac{(-1)^i\sqrt{(l+m)!(l-m)!(l+m')!
!    (l-m')!}}{(l-m'-i)!(l+m-i)!i!(i+m'-m)!}\\
!    &\times\left(\cos\frac{\beta}{2}\right)^{2l+m-m'-2i}\left(\sin\frac{\beta}
!    {2}\right)^{2i+m'-m}e^{-i(m\alpha+m'\gamma)},
!   \end{align*}
!   where the sum runs over all $i$ which make the factorial arguments
!   non-negative. For improper rotations, i.e. those which are a combination of
!   a rotation and inversion, the rotation is first made proper with
!   $R\rightarrow-R$ and $D$ is modified with
!   $D^l_{mm'}\rightarrow(-1)^l D^l_{mm'}$.
!
! !REVISION HISTORY:
!   Created October 2011 by DIN based on rotzflm.f90
!   Important changes: In rotzflm.f90 D^l_{mm'}(R^-1) is calculated
!                      I need D^l_{mm'}(R)
!EOP
!BOC
      implicit none
! arguments
      real (8), intent (in) :: rot(3, 3)
      integer, intent (in) :: l
      integer, intent (in) :: m1
      integer, intent (in) :: m2
! local variables
      integer :: i, j, nm, p
      real(8) :: det, roti(3, 3), ang(3)
      real(8) :: cb, sb, sum, t1, t2, t3
      real(8) :: par, alpha
      complex(8), parameter :: zzero = (0.d0, 0.d0)
      complex(8), parameter :: zone = (1.d0, 0.d0)

! external functions
      real (8) :: factnm
      external factnm
!
!     first check that the input parameters are physical
!
      if (l .lt. 0) then
         write(*,*)
         write(*, '("Error(getdlmm): l < 0 : ", i8)') l
         write(*,*)
         stop
      end if

! find the determinant
      det = rot(1,2)*rot(2,3)*rot(3,1)-rot(1,3)*rot(2,2)*rot(3,1)+ &
     &      rot(1,3)*rot(2,1)*rot(3,2)-rot(1,1)*rot(2,3)*rot(3,2)+ &
     &      rot(1,1)*rot(2,2)*rot(3,3)-rot(1,2)*rot(2,1)*rot(3,3)

! make rotation proper
      if (det .gt. 0.d0) then
         p = 1
         roti(:,:) =  rot(:,:)
      else
         p = -1
         roti(:,:) = -rot(:,:)
      end if

! compute Euler angles of rotation matrix
      call euler(roti,ang)
      cb = cos(ang(2)/2.d0)
      sb = sin(ang(2)/2.d0)

! generate rotation coefficient
      sum = 0.d0
      do i = 0, min(l+m1,l-m2)
         if (((l+m1-i) .ge. 0) .and. ((l-m2-i) .ge. 0) .and. &
        &    ((i+m2-m1) .ge. 0)) then
            j = 2 * l + m1 - m2 - 2 * i
            if (j .eq. 0) then
               t1 = 1.d0
            else
               t1 = cb ** j
            end if
            j = 2 * i + m2 - m1
            if (j .eq. 0) then
               t2 = 1.d0
            else
               t2 = sb ** j
            end if
            t3 = t1*t2/(factnm(l+m1-i,1)*factnm(l-m2-i,1)*factnm(i,1)*factnm(i+m2-m1,1))
            if (mod(i,2) .ne. 0) t3 = -t3
            sum = sum + t3
         end if
      end do
      t1 = sqrt(factnm(l+m1,1)*factnm(l-m1,1)*factnm(l+m2,1)*factnm(l-m2,1))
      t2 = -dble(m1)*ang(1)-dble(m2)*ang(3)
      getdlmm = sum*t1*cmplx(cos(t2),sin(t2),8)
      if ((p.eq.-1).and.(mod(l,2).ne.0)) getdlmm = -getdlmm
      
      return
end function
!EOC


subroutine testdlmm()
   use modmain
   implicit none

   real(8) :: c(3,3), ci(3,3)
   real(8) :: r0(3), r1(3)
   integer :: lmax, l
   integer :: m1, m2
   integer :: lm1, lm2
   integer :: isym, i, lspl

   complex(8), allocatable :: ylm0(:), ylm1(:), ylm2(:)
   complex(8) :: dlmm
   
   complex(8) :: getdlmm
   external      getdlmm

   lmax=4

   allocate(ylm0((lmax+1)*(lmax+1)))
   allocate(ylm1((lmax+1)*(lmax+1)))
   allocate(ylm2((lmax+1)*(lmax+1)))
   ylm0=zzero
   ylm1=zzero
   ylm2=zzero

   ! generate Ylm for the given r
   r0=(/2.d0,1.d0,1.d0/)
   call ylm(r0,lmax,ylm0)

   ! apply the test to the matrix c1
   !c(1,:)=(/ 0.d0,  0.d0, 1.d0/)
   !c(2,:)=(/-1.d0,  0.d0, 0.d0/)
   !c(3,:)=(/ 0.d0, -1.d0, 0.d0/)
   isym=9
   lspl=lsplsymc(isym)
   c(:,:)=symlatc(:,:,lspl)
   
   write(*,*)'Matrix c'
   do i=1,3
     write(*,'(3f12.4)') c(i,:)
   enddo
   
   ! get Ylm(R^(-1)r)
   call r3minv(c,ci)
   call r3mv(ci,r0,r1)
   call ylm(r1,lmax,ylm1)

   ! get Ylm using the rotation matrix D^l_{mm'}
   lm1=0
   do l=0,lmax
      
     write(*,*)
     do m1=-l,l
       lm1=lm1+1

       do m2=-l,l
         lm2=l*l+l+m2+1        
         dlmm=getdlmm(c,l,m2,m1)
         ylm2(lm1)=ylm2(lm1)+dlmm*ylm0(lm2)
       end do
       
       write(*,'(2i3,6f14.8)') l, m1, ylm1(lm1), ylm2(lm1), ylm0(lm1)

     end do
     
   end do

   deallocate(ylm0,ylm1,ylm2)

end subroutine
