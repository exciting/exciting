
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rotzflm
! !INTERFACE:
subroutine rotzflm(rot,lmax,n,ld,zflm1,zflm2)
! !INPUT/OUTPUT PARAMETERS:
!   rot   : rotation matrix (in,real(3,3))
!   lmax  : maximum angular momentum (in,integer)
!   n     : number of functions to rotate (in,integer)
!   ld    : leading dimension (in,integer)
!   zflm1 : coefficients of complex spherical harmonic expansion for each
!           function (in,complex(ld,n))
!   zflm2 : coefficients of rotated functions (out,complex(ld,n))
! !DESCRIPTION:
!   Rotates a set of functions
!   $$ f_i({\bf r})=\sum_{lm}f_{lm}^iY_{lm}(\hat{\bf r}) $$
!   for all $i$, given the coefficients $f_{lm}^i$ and a rotation matrix $R$.
!   This is done by first the computing the Euler angles $(\alpha,\beta,\gamma)$
!   of $R^{-1}$ (see routine {\tt euler}) and then generating the rotation
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
!   $D^l_{mm'}\rightarrow(-1)^l D^l_{mm'}$. The routine may be used in-place, in
!   other words, {\tt zflm1} and {\tt zflm2} can refer to the same array.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rot(3,3)
integer, intent(in) :: lmax
integer, intent(in) :: n
integer, intent(in) :: ld
complex(8), intent(in) :: zflm1(ld,n)
complex(8), intent(out) :: zflm2(ld,n)
! local variables
integer lmmax,l,m1,m2,lm1,lm2
integer i,j,nm,p
real(8) roti(3,3),ang(3),cb,sb
real(8) sum,t1,t2,t3
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
! allocatable arrays
complex(8), allocatable :: d(:,:)
complex(8), allocatable :: zflm3(:)
! external functions
real(8) r3mdet,factnm
external r3mdet,factnm
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(rotzflm): lmax < 0 : ",I8)') lmax
  write(*,*)
  stop
end if
if (n.eq.0) return
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(rotzflm): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
allocate(d(lmmax,lmmax))
allocate(zflm3(lmmax))
! invert r because the function is to be rotated and not the coordinate system
call r3minv(rot,roti)
! determine if the transformation is proper or improper
if (r3mdet(rot).gt.0.d0) then
  p=1
  roti(:,:)=roti(:,:)
else
  p=-1
  roti(:,:)=-roti(:,:)
end if
! compute Euler angles of rotation matrix
call euler(roti,ang)
cb=cos(ang(2)/2.d0)
sb=sin(ang(2)/2.d0)
lm1=0
do l=0,lmax
! generate rotation operator for m-components of current l
  do m1=-l,l
    lm1=lm1+1
    lm2=l**2
    do m2=-l,l
      lm2=lm2+1
      sum=0.d0
      do i=0,min(l+m1,l-m2)
        if (((l+m1-i).ge.0).and.((l-m2-i).ge.0).and.((i+m2-m1).ge.0)) then
          j=2*l+m1-m2-2*i
          if (j.eq.0) then
            t1=1.d0
          else
            t1=cb**j
          end if
          j=2*i+m2-m1
          if (j.eq.0) then
            t2=1.d0
          else
            t2=sb**j
          end if
          t3=t1*t2/(factnm(l+m1-i,1)*factnm(l-m2-i,1)*factnm(i,1) &
           *factnm(i+m2-m1,1))
          if (mod(i,2).ne.0) t3=-t3
          sum=sum+t3
        end if
      end do
      t1=sqrt(factnm(l+m1,1)*factnm(l-m1,1)*factnm(l+m2,1)*factnm(l-m2,1))
      t2=-dble(m1)*ang(1)-dble(m2)*ang(3)
      d(lm1,lm2)=sum*t1*cmplx(cos(t2),sin(t2),8)
      if ((p.eq.-1).and.(mod(l,2).ne.0)) d(lm1,lm2)=-d(lm1,lm2)
    end do
  end do
! apply rotation operator
  nm=2*l+1
  lm2=l**2+1
  do i=1,n
! make a copy in case zflm1 and zflm2 refer to the same array
    zflm3(lm2:lm2+nm-1)=zflm1(lm2:lm2+nm-1,i)
    call zgemv('N',nm,nm,zone,d(lm2,lm2),lmmax,zflm3(lm2),1,zzero, &
     zflm2(lm2,i),1)
  end do
end do
deallocate(d,zflm3)
return
end subroutine
!EOC
