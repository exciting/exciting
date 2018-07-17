module mod_optkgrid
  use modinput
  use mod_constants, only: pi
  implicit none

  real(8) :: brot(3,3), eye(3,3), vert(3,3), y0, z0
  integer :: on(3)
! methods
  contains

    subroutine getoptkgrid( r, bvec, kgrid, opt, ropt)
      ! This routine finds 3 integers to divide the reciprocal lattice vectors
      ! such that the points of the resulting Monkhorst-Pack grid are as equally
      ! spaced as possible. This is achieved when a k-point and its three nearest
      ! neighbors form a tetrahedron that is as regular as possible. The input 
      ! parameter r determines the aimed radius of the insphere of this tetrahedron.
      real(8), intent( in) :: r, bvec(3,3)
      integer, intent( out) :: kgrid(3)
      real(8), intent( out) :: opt, ropt

      integer :: i, j, k
      real(8) :: rot(3,3), m1(3,3), v1(3), v2(3), km(3,3), c, s, p, t, phi(3), reg(3,2), optmax

      eye = 0.d0
      eye(1,1) = 1.d0
      eye(2,2) = 1.d0
      eye(3,3) = 1.d0

      reg(1,:) = 0.5d0
      reg(2,1) = dsqrt( 3.d0)/2.d0
      reg(3,1) = 0.d0
      reg(2,2) = dsqrt( 3.d0)/6.d0
      reg(3,2) = dsqrt( 6.d0)/3.d0

      ! normalize
      brot(:,1) = bvec(:,1)/norm2( bvec(:,1))
      brot(:,2) = bvec(:,2)/norm2( bvec(:,2))
      brot(:,3) = bvec(:,3)/norm2( bvec(:,3))

      ! rotate first vector on x-axis
      call r3cross( brot(:,1), eye(:,1), v1)
      if( norm2( v1) .gt. 1.d-14) then
        v1 = v1/norm2( v1)
        c = dot_product( brot(:,1), eye(:,1))/norm2( brot(:,1))
        s = dsqrt( 1.d0 - c*c)
        km = 0.d0
        km(1,2) = -v1(3)
        km(2,1) =  v1(3)
        km(1,3) =  v1(2)
        km(3,1) = -v1(2)
        km(2,3) = -v1(1)
        km(3,2) =  v1(1)
        call r3mm( km, km, m1)
        rot = eye + s*km + (1.d0-c)*m1
        call r3mm( rot, brot, m1)
        brot = m1
      end if

      ! rotate second vector in x-y-plane
      v1 = eye(:,1)
      p = atan2( brot(3,2), brot(2,2))
      c = dcos( -p)
      s = dsin( -p)
      km = 0.d0
      km(1,2) = -v1(3)
      km(2,1) =  v1(3)
      km(1,3) =  v1(2)
      km(3,1) = -v1(2)
      km(2,3) = -v1(1)
      km(3,2) =  v1(1)
      call r3mm( km, km, m1)
      rot = eye + s*km + (1.d0-c)*m1
      call r3mm( rot, brot, m1)
      brot = m1

      ! map to correct octants
      if( abs( atan2( brot(2,2), brot(1,2))) .gt. 0.5d0*pi) brot(:,2) = -brot(:,2)
      if( brot(2,2) .lt. 0.d0) brot(2:3,:) = -brot(2:3,:)
      if( abs( atan2( brot(2,3), brot(1,3))) .gt. 0.5d0*pi) brot(:,3) = -brot(:,3)

      ! find starting vertices and onsets
      vert(:,1) = brot(:,1)
      reg(3,2) = sign( 1.d0, brot(3,3))*reg(3,2)
      y0 = reg(2,1)/brot(2,2)
      on(1) = nint( reg(1,1) - y0*brot(1,2))
      y0 = dot_product( reg(:,1), brot(:,2)) - on(1)*dot_product( brot(:,1), brot(:,2))
      vert(:,2) = on(1)*brot(:,1) + y0*brot(:,2)

      z0 = reg(3,2)/brot(3,3)
      v1 = reg(:,2) - z0*brot(:,3)
      on(2) = nint( v1(1) - v1(2)/vert(2,2)*vert(1,2))
      on(3) = nint( v1(2)/vert(2,2))
      z0 = dot_product( reg(:,2), brot(:,3)) - on(2)*dot_product( brot(:,1), brot(:,3)) - on(3)*dot_product( vert(:,2), brot(:,3))
      vert(:,3) = on(2)*brot(:,1) + on(3)*vert(:,2) + z0*brot(:,3)

      ! find optimal relations
      v1 = (/1.d0, y0, z0/)
      call optkgrid_findmin( v1, 10, (/-1.d0, 1.d0/), v2)
      call optkgrid_insphere( v2, p)
      optmax = 48.d0*dsqrt(3.d0)*p**3/(abs( brot(2,2)*brot(3,3)*v2(1)*v2(2)*v2(3)))
      
      ! find optimal integer grid
      v1(1) = norm2( bvec(:,1))*p/(v1(1)*r)
      v1(2) = norm2( bvec(:,2))*p/(v1(2)*r)
      v1(3) = norm2( bvec(:,3))*p/(v1(3)*r)
      m1(1,1) = dble( floor( v1(1)))
      m1(2,1) = dble( ceiling( v1(1)))
      m1(1,2) = dble( floor( v1(2)))
      m1(2,2) = dble( ceiling( v1(2)))
      m1(1,3) = dble( floor( v1(3)))
      m1(2,3) = dble( ceiling( v1(3)))
      s = 1.d100
      do i = 1, 2
        do j = 1, 2
          do k = 1, 2
            if( abs( m1(i,1)*m1(j,2)*m1(k,3)) .gt. 0.9d0) then
              v1(1) = norm2( bvec(:,1))/m1(i,1)
              v1(2) = norm2( bvec(:,2))/m1(j,2)
              v1(3) = norm2( bvec(:,3))/m1(k,3)
              call optkgrid_insphere( v1, p)
              c = 48.d0*dsqrt(3.d0)*p**3/(abs( brot(2,2)*brot(3,3)*v1(1)*v1(2)*v1(3)))
              t = abs( c - optmax)/optmax + abs( p - r)/r
              if( t .lt. s) then
                s = t
                opt = c
                ropt = p
                kgrid = abs( nint( (/m1(i,1), m1(j,2), m1(k,3)/)))
              end if
            end if
          end do
        end do
      end do
        
      return
    end subroutine getoptkgrid
    
    subroutine optkgrid_insphere( s, r)
      real(8), intent( in) :: s(3)
      real(8), intent( out) :: r

      real(8) :: n(3,3), v1(3), v2(3), minv(3,3)

      if( minval( abs( s)) .le. 1.d-23) then
        r = 0.d0
        return
      end if
      ! apply scaling
      vert(:,1) = s(1)*brot(:,1)
      vert(:,2) = on(1)*s(1)*brot(:,1) + s(2)*brot(:,2)
      vert(:,3) = on(2)*s(1)*brot(:,1) + on(3)*vert(:,2) + s(3)*brot(:,3)

      ! find outward normals
      call r3cross( vert(:,2)-vert(:,1), vert(:,3)-vert(:,1), n(:,1))
      n(:,1) = -n(:,1)*dot_product( n(:,1), -vert(:,1))
      n(:,1) = n(:,1)/norm2( n(:,1))
      call r3cross( vert(:,3)-vert(:,2), -vert(:,2), n(:,2))
      n(:,2) = -n(:,2)*dot_product( n(:,2), vert(:,1)-vert(:,2))
      n(:,2) = n(:,2)/norm2( n(:,2))
      call r3cross( -vert(:,3), vert(:,1)-vert(:,3), n(:,3))
      n(:,3) = -n(:,3)*dot_product( n(:,3), vert(:,2)-vert(:,3))
      n(:,3) = n(:,3)/norm2( n(:,3))

      ! solve linear system for center and radius
      v1(1) = dot_product( vert(:,1), n(:,1))
      v1(2) = dot_product( vert(:,2), n(:,2))
      v1(3) = dot_product( vert(:,3), n(:,3))
      n = transpose( n)
      n(:,3) = n(:,3) + sign( 1.d0, vert(3,3))
      call r3minv( n, minv)
      call r3mv( minv, v1, v2)

      r = abs( v2(3))
      return
    end subroutine optkgrid_insphere

    subroutine optkgrid_findmin( s0, n, intv0, s)
      integer, intent( in) :: n
      real(8), intent( in) :: s0(3), intv0(2)
      real(8), intent( out) :: s(3)

      real(8) :: w, st(3), intv(2,3)

      !st = s0
      !call optkgrid_findmin1d( 2, st, n, intv0, st(2))
      !call optkgrid_findmin1d( 3, st, n, intv0, st(3))
      
      st = s0
      s = 0.d0
      intv(:,1) = intv0
      intv(:,2) = intv0
      intv(:,3) = intv0

      do while( norm2( st - s) .gt. input%structure%epslat)
        s = st
        call optkgrid_findmin1d( 2, st, n, intv(:,2), w)
        st(2) = w
        call optkgrid_findmin1d( 3, st, n, intv(:,3), w)
        st(3) = w
        w = intv(2,2) - intv(1,2)
        intv(1,2) = -0.25d0*w
        intv(2,2) =  0.25d0*w
        intv(1,3) = -0.25d0*w
        intv(2,3) =  0.25d0*w
      end do

      s = st
      return
    end subroutine optkgrid_findmin

    subroutine optkgrid_findmin1d( d, s0, n, intv0, x)
      integer, intent( in) :: d, n
      real(8), intent( in) :: s0(3), intv0(2)
      real(8), intent( out) :: x

      integer :: i
      real(8) :: s(3), intv(2), dx, c, ddx, r, t, x0
      
      s = s0
      intv = intv0
      x0 = s0(d)
      x = x0

      ddx = 1.d0

      do while( ddx .gt. input%structure%epslat)
        c = 0.d0
        dx = (intv(2)-intv(1))/n

        do i = 0, n
          s(d) = x0 + intv(1) + i*dx
          call optkgrid_insphere( s, r)
          t = 48.d0*dsqrt(3.d0)*r**3/(abs( brot(2,2)*brot(3,3)*s(1)*s(2)*s(3)))
          if( t .gt. c) then
            ddx = abs( x - s(d))
            c = t
            x = s(d)
          end if
        end do
        x0 = x
        intv(1) = -dx
        intv(2) =  dx
      end do

      return
    end subroutine optkgrid_findmin1d

end module mod_optkgrid
