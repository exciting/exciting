!*******************************************************************************
!
!! BLEND_103 extends scalar point data into a cube.
!
!  Diagram:
!
!    011--------------111 
!      |               |
!      |               | 
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               | 
!      |               |
!    000--------------100 
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the 
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
subroutine blend103(xx,dens,grdens)
!----------------------------------
use param
!----------------------------------
implicit none

! input
real*8  :: xx(3)      ! input point
! output
real*8  :: dens       ! interpolated value of density
real*8  :: grdens     ! interpolated value of grad(n)

integer :: i,j,k
real*8  :: r,s,t
real*8  :: V000, V100, V010, V001, V110, V101, V011, V111

! location of indexes
  i = 1 + int(xx(1)*nx)
  j = 1 + int(xx(2)*ny)
  k = 1 + int(xx(3)*nz)
  !print*, 'i,j,k', i, j, k

! In the case xx==1.0 (i.e., i=nx+1, or j=ny+1, or k=nz+1)
  if (i == nx+1) i=nx
  if (j == ny+1) j=ny
  if (k == nz+1) k=nz
  
  r = xx(1)*dble(nx)-dble(i-1)
  s = xx(2)*dble(ny)-dble(j-1)
  t = xx(3)*dble(nz)-dble(k-1)
  !print*, 'r,s,t', r, s, t

! Density interpolation
  V000 = density(i  ,j  ,k  )
  V100 = density(i+1,j  ,k  )
  V010 = density(i  ,j+1,k  )
  V110 = density(i+1,j+1,k  )
  V001 = density(i  ,j  ,k+1)
  V101 = density(i+1,j  ,k+1)
  V011 = density(i  ,j+1,k+1)
  V111 = density(i+1,j+1,k+1)

  dens  = 								      &
    1	      * ( + V000 )						      &
  + r	      * ( - V000 + V100 )					      &
  +	s     * ( - V000	+ V010 )				      &
  +	    t * ( - V000	       + V001 ) 			      &
  + r * s     * ( + V000 - V100 - V010  		     + V110 )	      &
  + r	  * t * ( + V000 - V100        - V001	     + V101 )		      &
  +	s * t * ( + V000	- V010 - V001 + V011 )			      &
  + r * s * t * ( - V000 + V100 + V010 + V001 - V011 - V101 - V110 + V111 )

! Squared gradient density interpolation
  V000 = graddensity(i  ,j  ,k  )
  V100 = graddensity(i+1,j  ,k  )
  V010 = graddensity(i  ,j+1,k  )
  V110 = graddensity(i+1,j+1,k  )
  V001 = graddensity(i  ,j  ,k+1)
  V101 = graddensity(i+1,j  ,k+1)
  V011 = graddensity(i  ,j+1,k+1)
  V111 = graddensity(i+1,j+1,k+1)

  grdens  = 								      &
    1	      * ( + V000 )						      &
  + r	      * ( - V000 + V100 )					      & 
  +	s     * ( - V000	+ V010 )				      &
  +	    t * ( - V000	       + V001 ) 			      &
  + r * s     * ( + V000 - V100 - V010  		     + V110 )	      &
  + r	  * t * ( + V000 - V100        - V001	     + V101 )		      &
  +	s * t * ( + V000	- V010 - V001 + V011 )			      &
  + r * s * t * ( - V000 + V100 + V010 + V001 - V011 - V101 - V110 + V111 )

return
end subroutine blend103
