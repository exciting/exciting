      subroutine zxpay(n,dx,incx,da,dy,incy)
c
c     modified by:  D.R. Fokkema
c     01/06/94
c
c     a vector plus constant times a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      double complex dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = da*dy(iy) + dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = da*dy(i) + dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = da*dy(i) + dx(i)
        dy(i + 1) = da*dy(i + 1) + dx(i + 1)
        dy(i + 2) = da*dy(i + 2) + dx(i + 2)
        dy(i + 3) = da*dy(i + 3) + dx(i + 3)
   50 continue
      return
      end
