!BOP
!
! !ROUTINE: ylm
!
! !INTERFACE:
      subroutine ylm(v,lmax,y)
      
! !DESCRIPTION:
!
!  This subroutine calculates spherical harmonics up to l=lmax for a
! given vector in cartesian coordinates
!
!  The spherical harmonics (Condon and Shortley convention)
!  $Y_{0,0},Y_{1,-1},Y_{1,0},Y_{1,1},Y_{2,-2} ...Y_{LMAX,LMAX}$
!           for vector V (given in Cartesian coordinates)
!           are calculated. In the Condon Shortley convention the
!           spherical harmonics are defined as
!
!\begin{equation}                            
! Y_{l,m}=(-1)^m \sqrt{\frac{1}{2\pi}} P_l^m(cos(\theta))e^{im\phi}
!\end{equation}
!
!\noindent
!where $P_l^m(cos(\theta))$ is the normalized Associated Legendre
!           function. Thus,
!\begin{equation}                            
! Y_{l,-m}=(-1)^m Y^*_{l,m}
!\end{equation}
!
!The output is writen to the vector Y(:) such that
! Y(k)$=Y_{l,m}$ with $k=l^2+l+m+1$, thus,
! Y(1)$=Y_{0,0}$, Y(2)$=Y_{1,-1}$, Y(3)$=Y_{1,0}$, Y(4)$=Y_{1,1}$...
!Y(lmax*lmax+1)$=Y_{lmax,-lmax}... $Y((lmax+1)*(lmax+1))$=Y_{lmax,lmax}$
!
!   The basic algorithm used to calculate the spherical
!   harmonics for vector V is as follows:
!
!\begin{subequations}
!\begin{align}
! Y_{0,0}=&\sqrt{\frac{1}{4\pi}}\\
! Y_{1,0}=&\sqrt{\frac{3}{4\pi}}cos(\theta)\\
! Y_{1,1}=&-\sqrt{\frac{3}{8\pi}}sin(\theta)e^{i\phi}\\
! Y_{1,-1}=&-Y_{1,1}\\
! Y_{l,l}=&-\sqrt{\frac{2l+1}{2l}}sin(\theta)e^{i\phi}Y_{l-1,l-1}\\
! Y_{l,m}=&\sqrt{\frac{(2l-1)(2l+1)}{(l-m)(l+m)}}%
!cos(\theta)Y_{l-1,m}-\sqrt{\frac{(l-1+m)(l-1-m)(2l+1)}{(2l-3)(l-m)(l+m)}}%
!Y_{l-2,m}
!\end{align}
!\end{subequations}
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!

! !INPUT PARAMETERS:      

      implicit none

      integer(4) :: lmax ! maximum l for which the spherical 
!                          harmonics are calculated
!      
      real(8) :: v(3)    ! vector (cartesian coordinates) for which
!                          the spherical harmonics are calculated
!
! !OUTPUT PARAMETERS: 
      complex(8) :: y(*) ! the values of the spherical harmonics 
!                          dimension (lmax+1)^2
!      
!
!EOP
!
!BOC
      integer            i2l, i4l2, index, index2, l, m, msign
      double precision   a, b, c, ab, abc, abmax, abcmax
      double precision   d4ll1c, d2l13, pi
      double precision   costh, sinth, cosph, sinph
      double precision   temp1, temp2, temp3
      double precision   yllr, yll1r, yl1l1r, ylmr
      double precision   ylli, yll1i, yl1l1i, ylmi
!
      pi = (4.0d+0)*atan(1.0d+0)
!
!        y(0,0)
!
      yllr = 1.0d+0/sqrt(4.0d+0*pi)
      ylli = 0.0d+0
      y(1) = cmplx(yllr,ylli,8)
!
!        continue only if spherical harmonics for (l .gt. 0) are desired
!
      if (lmax .le. 0) goto 999
!
!        calculate sin(phi), cos(phi), sin(theta), cos(theta)
!        theta, phi ... polar angles of vector v
!
      abmax  = max(abs(v(1)),abs(v(2)))
      if (abmax .gt. 0.0d+0) then
         a = v(1)/abmax
         b = v(2)/abmax
         ab = sqrt(a*a+b*b)
         cosph = a/ab
         sinph = b/ab
      else
         cosph = 1.0d+0
         sinph = 0.0d+0
      endif
      abcmax = max(abmax,abs(v(3)))
      if (abcmax .gt. 0.0d+0) then
         a = v(1)/abcmax
         b = v(2)/abcmax
         c = v(3)/abcmax
         ab = a*a + b*b
         abc = sqrt(ab + c*c)
         costh = c/abc
         sinth = sqrt(ab)/abc
      else
         costh = 1.0d+0
         sinth = 0.0d+0
      endif
!
!        y(1,0)
!
      y(3) = cmplx(sqrt(3.0d+0)*yllr*costh,0.0d+0,8)
!
!        y(1,1) ( = -conjg(y(1,-1)))
!
      temp1 = -sqrt(1.5d+0)*yllr*sinth
      y(4) = cmplx(temp1*cosph,temp1*sinph,8)
      y(2) = -conjg(y(4))
!
      do 20 l = 2, lmax
         index  = l*l+1
         index2 = index + 2*l
         msign  = 1 - 2*mod(l,2)
!
!        yll = y(l,l) = f(y(l-1,l-1)) ... formula 1
!
         yl1l1r = dble(y(index-1))
         yl1l1i = aimag(y(index-1))
         temp1 = -sqrt(dble(2*l+1)/dble(2*l))*sinth
         yllr = temp1*(cosph*yl1l1r - sinph*yl1l1i)
         ylli = temp1*(cosph*yl1l1i + sinph*yl1l1r)
         y(index2) = cmplx(yllr,ylli,8)
         y(index)  = msign*conjg(y(index2))
         index2 = index2 - 1
         index  = index  + 1
!
!        yll1 = y(l,l-1) = f(y(l-1,l-1)) ... formula 2
!               (the coefficient for y(l-2,l-1) in formula 2 is zero)
!
         temp2 = sqrt(dble(2*l+1))*costh
         yll1r = temp2*yl1l1r
         yll1i = temp2*yl1l1i
         y(index2) = cmplx(yll1r,yll1i,8)
         y(index)  = -msign*conjg(y(index2))
         index2 = index2 - 1
         index  = index  + 1
!
         i4l2 = index2 - 4*l + 2
         i2l  = index2 - 2*l
         d4ll1c = costh*sqrt(dble(4*l*l-1))
         d2l13  = -sqrt(dble(2*l+1)/dble(2*l-3))
!
         do 10 m = l - 2, 0, -1
!
!        ylm = y(l,m) = f(y(l-2,m),y(l-1,m)) ... formula 2
!
            temp1 = 1.0d+0/sqrt(dble((l+m)*(l-m)))
            temp2 = d4ll1c*temp1
            temp3 = d2l13*sqrt(dble((l+m-1)*(l-m-1)))*temp1
            ylmr = temp2*dble(y(i2l))  + temp3*dble(y(i4l2))
            ylmi = temp2*aimag(y(i2l)) + temp3*aimag(y(i4l2))
            y(index2) = cmplx(ylmr,ylmi,8)
            y(index)  = msign*conjg(y(index2))
!
            msign  = -msign
            index2 = index2 - 1
            index  = index  + 1
            i4l2   = i4l2   - 1
            i2l    = i2l    - 1
   10    continue
   20 continue
!
  999 return
!
!        end of 'ylm'
!
      end
!EOC
