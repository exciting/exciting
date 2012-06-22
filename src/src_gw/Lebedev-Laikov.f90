!BOP
!
! !ROUTINE: gen_oh
!
! !INTERFACE:
       subroutine gen_oh(code, num, x, y, z, w, a, b, v)

! !INPUT PARAMETERS:
       
       implicit none

       integer(4) :: code ! Depending on code, there are 6...48 different 
!                           but equivalent points.
!                           code=1: (0,0,1) etc                        ( 6 p)
!                           code=2: (0,a,a) etc, a=1/sqrt(2)           (12 p)
!                           code=3: (a,a,a) etc, a=1/sqrt(3)           ( 8 p)
!                           code=4: (a,a,b) etc, b=sqrt(1-2 a^2)       (24 p)
!                           code=5: (a,b,0) etc, b=sqrt(1-a^2), a input(24 p)
!                           code=6: (a,b,c) etc, c=sqrt(1-a^2-b^2), 
!                                   a/b input  ( 48 p)
       
       integer(4) :: num ! number of equivalent points
       
       real(8) :: x(*)   ! x-coordinates of the gridpoints

       real(8) :: y(*)   ! y-coordinates of the gridpoints

       real(8) :: z(*)   ! z-coordinates of the gridpoints

       real(8) :: w(*)   ! weight of the grid point

       real(8) :: a      ! coordinates of the point in the sphere
       
       real(8) :: b      ! coordinates of the point in the sphere
       
       real(8) :: v      ! weight of the grid point

! !DESCRIPTION:
!
!\underline{Lebedev grids of orders n=6m+5 where m=0,1,...,21 in 16 digit
!precision}
!
!
!The file Lebedev-Laikov.F implements a set of subroutines providing 
!Lebedev-Laikov grids of order n=2m+1, where m=1,2,...,15, and additionally
!grids of order n=6m+5, where m=5,6,...,21. The parameters ensure 
!that angular integration of polynomials $x^k y^l z^m$, where $k+l+m \le
!131$ 
!can be performed with a relative accuracy of $2^{-14}$ [1]. Note that the weights
!are normalised to add up to 1.0.
!
!For each order n a separate subroutine is provided named 
!LD. The parameters X, Y, Z are arrays for the 
!cartesian components of each point, and the parameter W is an array for the
!!weights. The subroutines increase the integer(4) :: parameter N by number of grid
!points generated. All these routines use the subroutine gen\_oh which takes care 
!of the octahedral symmetry of the grids.
!
!Christoph van Wuellen (Ruhr-Universitaet, Bochum, Germany) generated the 
!routines in Lebedev-Laikov.F by translating the original C-routines kindly 
!provided by Dmitri Laikov (Moscow State University, Moscow, Russia). We 
!are in debt to Dmitri Laikov for giving us permission to make these routines
!publically available.
!
!   This subroutine is part of a set of subroutines that generate
!   Lebedev grids [1-6] for integration on a sphere. The original 
!   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!   translated into fortran by Dr. Christoph van Wuellen.
!   This subroutine was translated from C to fortran77 by hand.
!
!   Users of this code are asked to include reference [1] in their
!   publications, and in the user- and programmers-manuals 
!   describing their codes.
!
!   This code was distributed through CCL (http://www.ccl.net/).
!
!\begin{enumerate}
!   \item V.I. Lebedev, and D.N. Laikov
!       "A quadrature formula for the sphere of the 131st
!        algebraic order of accuracy"
!       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
!   \item V.I. Lebedev
!       "A quadrature formula for the sphere of 59th algebraic
!        order of accuracy"
!       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!
!   \item V.I. Lebedev, and A.L. Skorokhodov
!       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!
!   \item V.I. Lebedev
!       "Spherical quadrature formulas exact to orders 25-29"
!       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!
!   \item V.I. Lebedev
!       "Quadratures on a sphere"
!       Computational Mathematics and Mathematical Physics, Vol. 16,
!       1976, pp. 10-24. 
!
!   \item V.I. Lebedev
!       "Values of the nodes and weights of ninth to seventeenth 
!        order Gauss-Markov quadrature formulae invariant under the
!        octahedron group with inversion"
!       Computational Mathematics and Mathematical Physics, Vol. 15,
!       1975, pp. 44-51.
!\end{enumerate}
!
!
!    Given a point on a sphere (specified by a and b), generate all
!    the equivalent points under Oh symmetry, making grid points with
!    weight v.
!    The variable num is increased by the number of different points
!    generated.
!
!{\sf Author: } Huub van Dam
!
!Daresbury Laboratory, Daresbury, United Kingdom
!
!
! !REVISION HISTORY:
!  Created: April, 2000       
!  Last modified: Feb. 2004 (RGA)
!

! !LOCAL VARIABLES:

       real(8) :: c ! for code = 6 => c=sqrt(1-a^2-b^2)
       
!EOP       
!
!BOC
       select case (code)

         case(1)

           a=1.0d0
           x(1) =  a
           y(1) =  0.0d0
           z(1) =  0.0d0
           w(1) =  v
           x(2) = -a
           y(2) =  0.0d0
           z(2) =  0.0d0
           w(2) =  v
           x(3) =  0.0d0
           y(3) =  a
           z(3) =  0.0d0
           w(3) =  v
           x(4) =  0.0d0
           y(4) = -a
           z(4) =  0.0d0
           w(4) =  v
           x(5) =  0.0d0
           y(5) =  0.0d0
           z(5) =  a
           w(5) =  v
           x(6) =  0.0d0
           y(6) =  0.0d0
           z(6) = -a
           w(6) =  v
           num=num+6

         case (2)

           a=sqrt(0.5d0)
           x( 1) =  0d0
           y( 1) =  a
           z( 1) =  a
           w( 1) =  v
           x( 2) =  0d0
           y( 2) = -a
           z( 2) =  a
           w( 2) =  v
           x( 3) =  0d0
           y( 3) =  a
           z( 3) = -a
           w( 3) =  v
           x( 4) =  0d0
           y( 4) = -a
           z( 4) = -a
           w( 4) =  v
           x( 5) =  a
           y( 5) =  0d0
           z( 5) =  a
           w( 5) =  v
           x( 6) = -a
           y( 6) =  0d0
           z( 6) =  a
           w( 6) =  v
           x( 7) =  a
           y( 7) =  0d0
           z( 7) = -a
           w( 7) =  v
           x( 8) = -a
           y( 8) =  0d0
           z( 8) = -a
           w( 8) =  v
           x( 9) =  a
           y( 9) =  a
           z( 9) =  0d0
           w( 9) =  v
           x(10) = -a
           y(10) =  a
           z(10) =  0d0
           w(10) =  v
           x(11) =  a
           y(11) = -a
           z(11) =  0d0
           w(11) =  v
           x(12) = -a
           y(12) = -a
           z(12) =  0d0
           w(12) =  v
           num=num+12

         case(3)

           a = sqrt(1d0/3d0)
           x(1) =  a
           y(1) =  a
           z(1) =  a
           w(1) =  v
           x(2) = -a
           y(2) =  a
           z(2) =  a
           w(2) =  v
           x(3) =  a
           y(3) = -a
           z(3) =  a
           w(3) =  v
           x(4) = -a
           y(4) = -a
           z(4) =  a
           w(4) =  v
           x(5) =  a
           y(5) =  a
           z(5) = -a
           w(5) =  v
           x(6) = -a
           y(6) =  a
           z(6) = -a
           w(6) =  v
           x(7) =  a
           y(7) = -a
           z(7) = -a
           w(7) =  v
           x(8) = -a
           y(8) = -a
           z(8) = -a
           w(8) =  v
           num=num+8

         case(4)
    
           b = sqrt(1d0 - 2d0*a*a)
           x( 1) =  a
           y( 1) =  a
           z( 1) =  b
           w( 1) =  v
           x( 2) = -a
           y( 2) =  a
           z( 2) =  b
           w( 2) =  v
           x( 3) =  a
           y( 3) = -a
           z( 3) =  b
           w( 3) =  v
           x( 4) = -a
           y( 4) = -a
           z( 4) =  b
           w( 4) =  v
           x( 5) =  a
           y( 5) =  a
           z( 5) = -b
           w( 5) =  v
           x( 6) = -a
           y( 6) =  a
           z( 6) = -b
           w( 6) =  v
           x( 7) =  a
           y( 7) = -a
           z( 7) = -b
           w( 7) =  v
           x( 8) = -a
           y( 8) = -a
           z( 8) = -b
           w( 8) =  v
           x( 9) =  a
           y( 9) =  b
           z( 9) =  a
           w( 9) =  v
           x(10) = -a
           y(10) =  b
           z(10) =  a
           w(10) =  v
           x(11) =  a
           y(11) = -b
           z(11) =  a
           w(11) =  v
           x(12) = -a
           y(12) = -b
           z(12) =  a
           w(12) =  v
           x(13) =  a
           y(13) =  b
           z(13) = -a
           w(13) =  v
           x(14) = -a
           y(14) =  b
           z(14) = -a
           w(14) =  v
           x(15) =  a
           y(15) = -b
           z(15) = -a
           w(15) =  v
           x(16) = -a
           y(16) = -b
           z(16) = -a
           w(16) =  v
           x(17) =  b
           y(17) =  a
           z(17) =  a
           w(17) =  v
           x(18) = -b
           y(18) =  a
           z(18) =  a
           w(18) =  v
           x(19) =  b
           y(19) = -a
           z(19) =  a
           w(19) =  v
           x(20) = -b
           y(20) = -a
           z(20) =  a
           w(20) =  v
           x(21) =  b
           y(21) =  a
           z(21) = -a
           w(21) =  v
           x(22) = -b
           y(22) =  a
           z(22) = -a
           w(22) =  v
           x(23) =  b
           y(23) = -a
           z(23) = -a
           w(23) =  v
           x(24) = -b
           y(24) = -a
           z(24) = -a
           w(24) =  v
           num=num+24

         case(5)

           b=sqrt(1d0-a*a)
           x( 1) =  a
           y( 1) =  b
           z( 1) =  0d0
           w( 1) =  v
           x( 2) = -a
           y( 2) =  b
           z( 2) =  0d0
           w( 2) =  v
           x( 3) =  a
           y( 3) = -b
           z( 3) =  0d0
           w( 3) =  v
           x( 4) = -a
           y( 4) = -b
           z( 4) =  0d0
           w( 4) =  v
           x( 5) =  b
           y( 5) =  a
           z( 5) =  0d0
           w( 5) =  v
           x( 6) = -b
           y( 6) =  a
           z( 6) =  0d0
           w( 6) =  v
           x( 7) =  b
           y( 7) = -a
           z( 7) =  0d0
           w( 7) =  v
           x( 8) = -b
           y( 8) = -a
           z( 8) =  0d0
           w( 8) =  v
           x( 9) =  a
           y( 9) =  0d0
           z( 9) =  b
           w( 9) =  v
           x(10) = -a
           y(10) =  0d0
           z(10) =  b
           w(10) =  v
           x(11) =  a
           y(11) =  0d0
           z(11) = -b
           w(11) =  v
           x(12) = -a
           y(12) =  0d0
           z(12) = -b
           w(12) =  v
           x(13) =  b
           y(13) =  0d0
           z(13) =  a
           w(13) =  v
           x(14) = -b
           y(14) =  0d0
           z(14) =  a
           w(14) =  v
           x(15) =  b
           y(15) =  0d0
           z(15) = -a
           w(15) =  v
           x(16) = -b
           y(16) =  0d0
           z(16) = -a
           w(16) =  v
           x(17) =  0d0
           y(17) =  a
           z(17) =  b
           w(17) =  v
           x(18) =  0d0
           y(18) = -a
           z(18) =  b
           w(18) =  v
           x(19) =  0d0
           y(19) =  a
           z(19) = -b
           w(19) =  v
           x(20) =  0d0
           y(20) = -a
           z(20) = -b
           w(20) =  v
           x(21) =  0d0
           y(21) =  b
           z(21) =  a
           w(21) =  v
           x(22) =  0d0
           y(22) = -b
           z(22) =  a
           w(22) =  v
           x(23) =  0d0
           y(23) =  b
           z(23) = -a
           w(23) =  v
           x(24) =  0d0
           y(24) = -b
           z(24) = -a
           w(24) =  v
           num=num+24
    
         case(6)

           c=sqrt(1d0 - a*a - b*b)
           x( 1) =  a
           y( 1) =  b
           z( 1) =  c
           w( 1) =  v
           x( 2) = -a
           y( 2) =  b
           z( 2) =  c
           w( 2) =  v
           x( 3) =  a
           y( 3) = -b
           z( 3) =  c
           w( 3) =  v
           x( 4) = -a
           y( 4) = -b
           z( 4) =  c
           w( 4) =  v
           x( 5) =  a
           y( 5) =  b
           z( 5) = -c
           w( 5) =  v
           x( 6) = -a
           y( 6) =  b
           z( 6) = -c
           w( 6) =  v
           x( 7) =  a
           y( 7) = -b
           z( 7) = -c
           w( 7) =  v
           x( 8) = -a
           y( 8) = -b
           z( 8) = -c
           w( 8) =  v
           x( 9) =  a
           y( 9) =  c
           z( 9) =  b
           w( 9) =  v
           x(10) = -a
           y(10) =  c
           z(10) =  b
           w(10) =  v
           x(11) =  a
           y(11) = -c
           z(11) =  b
           w(11) =  v
           x(12) = -a
           y(12) = -c
           z(12) =  b
           w(12) =  v
           x(13) =  a
           y(13) =  c
           z(13) = -b
           w(13) =  v
           x(14) = -a
           y(14) =  c
           z(14) = -b
           w(14) =  v
           x(15) =  a
           y(15) = -c
           z(15) = -b
           w(15) =  v
           x(16) = -a
           y(16) = -c
           z(16) = -b
           w(16) =  v
           x(17) =  b
           y(17) =  a
           z(17) =  c
           w(17) =  v
           x(18) = -b
           y(18) =  a
           z(18) =  c
           w(18) =  v
           x(19) =  b
           y(19) = -a
           z(19) =  c
           w(19) =  v
           x(20) = -b
           y(20) = -a
           z(20) =  c
           w(20) =  v
           x(21) =  b
           y(21) =  a
           z(21) = -c
           w(21) =  v
           x(22) = -b
           y(22) =  a
           z(22) = -c
           w(22) =  v
           x(23) =  b
           y(23) = -a
           z(23) = -c
           w(23) =  v
           x(24) = -b
           y(24) = -a
           z(24) = -c
           w(24) =  v
           x(25) =  b
           y(25) =  c
           z(25) =  a
           w(25) =  v
           x(26) = -b
           y(26) =  c
           z(26) =  a
           w(26) =  v
           x(27) =  b
           y(27) = -c
           z(27) =  a
           w(27) =  v
           x(28) = -b
           y(28) = -c
           z(28) =  a
           w(28) =  v
           x(29) =  b
           y(29) =  c
           z(29) = -a
           w(29) =  v
           x(30) = -b
           y(30) =  c
           z(30) = -a
           w(30) =  v
           x(31) =  b
           y(31) = -c
           z(31) = -a
           w(31) =  v
           x(32) = -b
           y(32) = -c
           z(32) = -a
           w(32) =  v
           x(33) =  c
           y(33) =  a
           z(33) =  b
           w(33) =  v
           x(34) = -c
           y(34) =  a
           z(34) =  b
           w(34) =  v
           x(35) =  c
           y(35) = -a
           z(35) =  b
           w(35) =  v
           x(36) = -c
           y(36) = -a
           z(36) =  b
           w(36) =  v
           x(37) =  c
           y(37) =  a
           z(37) = -b
           w(37) =  v
           x(38) = -c
           y(38) =  a
           z(38) = -b
           w(38) =  v
           x(39) =  c
           y(39) = -a
           z(39) = -b
           w(39) =  v
           x(40) = -c
           y(40) = -a
           z(40) = -b
           w(40) =  v
           x(41) =  c
           y(41) =  b
           z(41) =  a
           w(41) =  v
           x(42) = -c
           y(42) =  b
           z(42) =  a
           w(42) =  v
           x(43) =  c
           y(43) = -b
           z(43) =  a
           w(43) =  v
           x(44) = -c
           y(44) = -b
           z(44) =  a
           w(44) =  v
           x(45) =  c
           y(45) =  b
           z(45) = -a
           w(45) =  v
           x(46) = -c
           y(46) =  b
           z(46) = -a
           w(46) =  v
           x(47) =  c
           y(47) = -b
           z(47) = -a
           w(47) =  v
           x(48) = -c
           y(48) = -b
           z(48) = -a
           w(48) =  v
           num=num+48

         case default

           write (*,*) 'Gen_Oh: Invalid Code'
           stop 

       end select

       return

       end
!EOC
!
!BOP
!
! !ROUTINE: ldmmmm(x,y,z)
!
! !INTERFACE:
!       subroutine ldmmmm(x,y,z,w,n)
!
! !INPUT PARAMETERS:
!
!       real(8) :: x(mmmm) ! x-coordinates of the gridpoints
!
!       real(8) :: y(mmmm) ! y-coordinates of the gridpoints
!
!       real(8) :: z(mmmm) ! z-coordinates of the gridpoints
!
!       real(8) :: w(mmmm) ! weight of the grid point
!
!       integer(4) :: n    ! number of grid points
!
! !DESCRIPTION:
!
!  \underline{Lebedev mmmm-point angular grid}
!
!   This is a set of subroutines that generate
!   Lebedev grids [1-6] for integration on a sphere. The original 
!   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!   translated into fortran by Dr. Christoph van Wuellen.
!   This subroutine was translated using a C to fortran77 conversion
!   tool written by Dr. Christoph van Wuellen.
!
!   mmmm stands for the number of grid points, its possible values are
!   6, 14, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590,
!   770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334,
!   4802, 5294, 5810.
!
!   Users of this code are asked to include reference [1] in their
!   publications, and in the user- and programmers-manuals 
!   describing their codes.
!
!   This code was distributed through CCL (http://www.ccl.net/).
!
!\begin{enumerate}
!   \item V.I. Lebedev, and D.N. Laikov
!       "A quadrature formula for the sphere of the 131st
!        algebraic order of accuracy"
!       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
!   \item V.I. Lebedev
!       "A quadrature formula for the sphere of 59th algebraic
!        order of accuracy"
!       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!
!   \item V.I. Lebedev, and A.L. Skorokhodov
!       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!
!   \item V.I. Lebedev
!       "Spherical quadrature formulas exact to orders 25-29"
!       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!
!   \item V.I. Lebedev
!       "Quadratures on a sphere"
!       Computational Mathematics and Mathematical Physics, Vol. 16,
!       1976, pp. 10-24. 
!
!   \item V.I. Lebedev
!       "Values of the nodes and weights of ninth to seventeenth 
!        order Gauss-Markov quadrature formulae invariant under the
!        octahedron group with inversion"
!       Computational Mathematics and Mathematical Physics, Vol. 15,
!       1975, pp. 44-51.
!\end{enumerate}
!
! !LOCAL VARIABLES:
!
!       real(8) :: a,b,v
!
!EOP
!
!BOC       
! 
!   Only the 6- and 14-point routines are shown as an example
!
!    Lebedev 6-point angular grid
!
       subroutine ld0006(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(   6)
! y-coordinates of the gridpoints
       real(8) :: y(   6)
! z-coordinates of the gridpoints
       real(8) :: z(   6)
! weight of the grid point
       real(8) :: w(   6)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1666666666666667d+0
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 14-point angular grid
!
       subroutine ld0014(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(  14)
! y-coordinates of the gridpoints
       real(8) :: y(  14)
! z-coordinates of the gridpoints
       real(8) :: z(  14)
! weight of the grid point
       real(8) :: w(  14)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.6666666666666667d-1
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.7500000000000000d-1
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!EOC
!
!    Lebedev 26-point angular grid
!
       subroutine ld0026(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(  26)
! y-coordinates of the gridpoints
       real(8) :: y(  26)
! z-coordinates of the gridpoints
       real(8) :: z(  26)
! weight of the grid point
       real(8) :: w(  26)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.4761904761904762d-1
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3809523809523810d-1
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3214285714285714d-1
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 38-point angular grid
!
       subroutine ld0038(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(  38)
! y-coordinates of the gridpoints
       real(8) :: y(  38)
! z-coordinates of the gridpoints
       real(8) :: z(  38)
! weight of the grid point
       real(8) :: w(  38)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v

       n=1
       v=0.9523809523809524d-2
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3214285714285714d-1
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4597008433809831d+0
       v=0.2857142857142857d-1
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 50-point angular grid
!
       subroutine ld0050(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(  50)
! y-coordinates of the gridpoints
       real(8) :: y(  50)
! z-coordinates of the gridpoints
       real(8) :: z(  50)
! weight of the grid point
       real(8) :: w(  50)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1269841269841270d-1
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2257495590828924d-1
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2109375000000000d-1
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3015113445777636d+0
       v=0.2017333553791887d-1
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 74-point angular grid
!
       subroutine ld0074(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(  74)
! y-coordinates of the gridpoints
       real(8) :: y(  74)
! z-coordinates of the gridpoints
       real(8) :: z(  74)
! weight of the grid point
       real(8) :: w(  74)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.5130671797338464d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1660406956574204d-1
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=-0.2958603896103896d-1
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4803844614152614d+0
       v=0.2657620708215946d-1
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3207726489807764d+0
       v=0.1652217099371571d-1
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 86-point angular grid
!
       subroutine ld0086(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(  86)
! y-coordinates of the gridpoints
       real(8) :: y(  86)
! z-coordinates of the gridpoints
       real(8) :: z(  86)
! weight of the grid point
       real(8) :: w(  86)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1154401154401154d-1
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1194390908585628d-1
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3696028464541502d+0
       v=0.1111055571060340d-1
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6943540066026664d+0
       v=0.1187650129453714d-1
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3742430390903412d+0
       v=0.1181230374690448d-1
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 110-point angular grid
!
       subroutine ld0110(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 110)
! y-coordinates of the gridpoints
       real(8) :: y( 110)
! z-coordinates of the gridpoints
       real(8) :: z( 110)
! weight of the grid point
       real(8) :: w( 110)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.3828270494937162d-2
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.9793737512487512d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1851156353447362d+0
       v=0.8211737283191111d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6904210483822922d+0
       v=0.9942814891178103d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3956894730559419d+0
       v=0.9595471336070963d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4783690288121502d+0
       v=0.9694996361663028d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 146-point angular grid
!
       subroutine ld0146(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 146)
! y-coordinates of the gridpoints
       real(8) :: y( 146)
! z-coordinates of the gridpoints
       real(8) :: z( 146)
! weight of the grid point
       real(8) :: w( 146)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.5996313688621381d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.7372999718620756d-2
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.7210515360144488d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6764410400114264d+0
       v=0.7116355493117555d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4174961227965453d+0
       v=0.6753829486314477d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1574676672039082d+0
       v=0.7574394159054034d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1403553811713183d+0
       b=0.4493328323269557d+0
       v=0.6991087353303262d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 170-point angular grid
!
       subroutine ld0170(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 170)
! y-coordinates of the gridpoints
       real(8) :: y( 170)
! z-coordinates of the gridpoints
       real(8) :: z( 170)
! weight of the grid point
       real(8) :: w( 170)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.5544842902037365d-2
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.6071332770670752d-2
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.6383674773515093d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2551252621114134d+0
       v=0.5183387587747790d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6743601460362766d+0
       v=0.6317929009813725d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4318910696719410d+0
       v=0.6201670006589077d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2613931360335988d+0
       v=0.5477143385137348d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4990453161796037d+0
       b=0.1446630744325115d+0
       v=0.5968383987681156d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 194-point angular grid
!
       subroutine ld0194(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 194)
! y-coordinates of the gridpoints
       real(8) :: y( 194)
! z-coordinates of the gridpoints
       real(8) :: z( 194)
! weight of the grid point
       real(8) :: w( 194)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1782340447244611d-2
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.5716905949977102d-2
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.5573383178848738d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6712973442695226d+0
       v=0.5608704082587997d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2892465627575439d+0
       v=0.5158237711805383d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4446933178717437d+0
       v=0.5518771467273614d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1299335447650067d+0
       v=0.4106777028169394d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3457702197611283d+0
       v=0.5051846064614808d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1590417105383530d+0
       b=0.8360360154824589d+0
       v=0.5530248916233094d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 230-point angular grid
!
       subroutine ld0230(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 230)
! y-coordinates of the gridpoints
       real(8) :: y( 230)
! z-coordinates of the gridpoints
       real(8) :: z( 230)
! weight of the grid point
       real(8) :: w( 230)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=-0.5522639919727325d-1
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.4450274607445226d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4492044687397611d+0
       v=0.4496841067921404d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2520419490210201d+0
       v=0.5049153450478750d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6981906658447242d+0
       v=0.3976408018051883d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6587405243460960d+0
       v=0.4401400650381014d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4038544050097660d-1
       v=0.1724544350544401d-1
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5823842309715585d+0
       v=0.4231083095357343d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3545877390518688d+0
       v=0.5198069864064399d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2272181808998187d+0
       b=0.4864661535886647d+0
       v=0.4695720972568883d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 266-point angular grid
!
       subroutine ld0266(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 266)
! y-coordinates of the gridpoints
       real(8) :: y( 266)
! z-coordinates of the gridpoints
       real(8) :: z( 266)
! weight of the grid point
       real(8) :: w( 266)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=-0.1313769127326952d-2
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=-0.2522728704859336d-2
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.4186853881700583d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7039373391585475d+0
       v=0.5315167977810885d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1012526248572414d+0
       v=0.4047142377086219d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4647448726420539d+0
       v=0.4112482394406990d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3277420654971629d+0
       v=0.3595584899758782d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6620338663699974d+0
       v=0.4256131351428158d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8506508083520399d+0
       v=0.4229582700647240d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3233484542692899d+0
       b=0.1153112011009701d+0
       v=0.4080914225780505d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2314790158712601d+0
       b=0.5244939240922365d+0
       v=0.4071467593830964d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 302-point angular grid
!
       subroutine ld0302(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 302)
! y-coordinates of the gridpoints
       real(8) :: y( 302)
! z-coordinates of the gridpoints
       real(8) :: z( 302)
! weight of the grid point
       real(8) :: w( 302)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.8545911725128148d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3599119285025571d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3515640345570105d+0
       v=0.3449788424305883d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6566329410219612d+0
       v=0.3604822601419882d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4729054132581005d+0
       v=0.3576729661743367d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9618308522614784d-1
       v=0.2352101413689164d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2219645236294178d+0
       v=0.3108953122413675d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7011766416089545d+0
       v=0.3650045807677255d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2644152887060663d+0
       v=0.2982344963171804d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5718955891878961d+0
       v=0.3600820932216460d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2510034751770465d+0
       b=0.8000727494073952d+0
       v=0.3571540554273387d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1233548532583327d+0
       b=0.4127724083168531d+0
       v=0.3392312205006170d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 350-point angular grid
!
       subroutine ld0350(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 350)
! y-coordinates of the gridpoints
       real(8) :: y( 350)
! z-coordinates of the gridpoints
       real(8) :: z( 350)
! weight of the grid point
       real(8) :: w( 350)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.3006796749453936d-2
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3050627745650771d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7068965463912316d+0
       v=0.1621104600288991d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4794682625712025d+0
       v=0.3005701484901752d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1927533154878019d+0
       v=0.2990992529653774d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6930357961327123d+0
       v=0.2982170644107595d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3608302115520091d+0
       v=0.2721564237310992d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6498486161496169d+0
       v=0.3033513795811141d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1932945013230339d+0
       v=0.3007949555218533d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3800494919899303d+0
       v=0.2881964603055307d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2899558825499574d+0
       b=0.7934537856582316d+0
       v=0.2958357626535696d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9684121455103957d-1
       b=0.8280801506686862d+0
       v=0.3036020026407088d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1833434647041659d+0
       b=0.9074658265305127d+0
       v=0.2832187403926303d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 434-point angular grid
!
       subroutine ld0434(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 434)
! y-coordinates of the gridpoints
       real(8) :: y( 434)
! z-coordinates of the gridpoints
       real(8) :: z( 434)
! weight of the grid point
       real(8) :: w( 434)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.5265897968224436d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2548219972002607d-2
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2512317418927307d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6909346307509111d+0
       v=0.2530403801186355d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1774836054609158d+0
       v=0.2014279020918528d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4914342637784746d+0
       v=0.2501725168402936d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6456664707424256d+0
       v=0.2513267174597564d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2861289010307638d+0
       v=0.2302694782227416d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7568084367178018d-1
       v=0.1462495621594614d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3927259763368002d+0
       v=0.2445373437312980d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8818132877794288d+0
       v=0.2417442375638981d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9776428111182649d+0
       v=0.1910951282179532d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2054823696403044d+0
       b=0.8689460322872412d+0
       v=0.2416930044324775d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5905157048925271d+0
       b=0.7999278543857286d+0
       v=0.2512236854563495d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5550152361076807d+0
       b=0.7717462626915901d+0
       v=0.2496644054553086d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9371809858553722d+0
       b=0.3344363145343455d+0
       v=0.2236607760437849d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 590-point angular grid
!
       subroutine ld0590(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 590)
! y-coordinates of the gridpoints
       real(8) :: y( 590)
! z-coordinates of the gridpoints
       real(8) :: z( 590)
! weight of the grid point
       real(8) :: w( 590)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.3095121295306187d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1852379698597489d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7040954938227469d+0
       v=0.1871790639277744d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6807744066455243d+0
       v=0.1858812585438317d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6372546939258752d+0
       v=0.1852028828296213d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5044419707800358d+0
       v=0.1846715956151242d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4215761784010967d+0
       v=0.1818471778162769d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3317920736472123d+0
       v=0.1749564657281154d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2384736701421887d+0
       v=0.1617210647254411d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1459036449157763d+0
       v=0.1384737234851692d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6095034115507196d-1
       v=0.9764331165051050d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6116843442009876d+0
       v=0.1857161196774078d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3964755348199858d+0
       v=0.1705153996395864d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1724782009907724d+0
       v=0.1300321685886048d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5610263808622060d+0
       b=0.3518280927733519d+0
       v=0.1842866472905286d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4742392842551980d+0
       b=0.2634716655937950d+0
       v=0.1802658934377451d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5984126497885380d+0
       b=0.1816640840360209d+0
       v=0.1849830560443660d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3791035407695563d+0
       b=0.1720795225656878d+0
       v=0.1713904507106709d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2778673190586244d+0
       b=0.8213021581932511d-1
       v=0.1555213603396808d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5033564271075117d+0
       b=0.8999205842074875d-1
       v=0.1802239128008525d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 770-point angular grid
!
       subroutine ld0770(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 770)
! y-coordinates of the gridpoints
       real(8) :: y( 770)
! z-coordinates of the gridpoints
       real(8) :: z( 770)
! weight of the grid point
       real(8) :: w( 770)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.2192942088181184d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1436433617319080d-2
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1421940344335877d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5087204410502360d-1
       v=0.6798123511050502d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1228198790178831d+0
       v=0.9913184235294912d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2026890814408786d+0
       v=0.1180207833238949d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2847745156464294d+0
       v=0.1296599602080921d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3656719078978026d+0
       v=0.1365871427428316d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4428264886713469d+0
       v=0.1402988604775325d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5140619627249735d+0
       v=0.1418645563595609d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6306401219166803d+0
       v=0.1421376741851662d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6716883332022612d+0
       v=0.1423996475490962d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6979792685336881d+0
       v=0.1431554042178567d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1446865674195309d+0
       v=0.9254401499865368d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3390263475411216d+0
       v=0.1250239995053509d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5335804651263506d+0
       v=0.1394365843329230d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6944024393349413d-1
       b=0.2355187894242326d+0
       v=0.1127089094671749d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2269004109529460d+0
       b=0.4102182474045730d+0
       v=0.1345753760910670d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8025574607775339d-1
       b=0.6214302417481605d+0
       v=0.1424957283316783d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1467999527896572d+0
       b=0.3245284345717394d+0
       v=0.1261523341237750d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1571507769824727d+0
       b=0.5224482189696630d+0
       v=0.1392547106052696d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2365702993157246d+0
       b=0.6017546634089558d+0
       v=0.1418761677877656d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7714815866765732d-1
       b=0.4346575516141163d+0
       v=0.1338366684479554d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3062936666210730d+0
       b=0.4908826589037616d+0
       v=0.1393700862676131d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3822477379524787d+0
       b=0.5648768149099500d+0
       v=0.1415914757466932d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 974-point angular grid
!
       subroutine ld0974(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x( 974)
! y-coordinates of the gridpoints
       real(8) :: y( 974)
! z-coordinates of the gridpoints
       real(8) :: z( 974)
! weight of the grid point
       real(8) :: w( 974)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1438294190527431d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1125772288287004d-2
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4292963545341347d-1
       v=0.4948029341949241d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1051426854086404d+0
       v=0.7357990109125470d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1750024867623087d+0
       v=0.8889132771304384d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2477653379650257d+0
       v=0.9888347838921435d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3206567123955957d+0
       v=0.1053299681709471d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3916520749849983d+0
       v=0.1092778807014578d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4590825874187624d+0
       v=0.1114389394063227d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5214563888415861d+0
       v=0.1123724788051555d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6253170244654199d+0
       v=0.1125239325243814d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6637926744523170d+0
       v=0.1126153271815905d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6910410398498301d+0
       v=0.1130286931123841d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7052907007457760d+0
       v=0.1134986534363955d-2
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1236686762657990d+0
       v=0.6823367927109931d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2940777114468387d+0
       v=0.9454158160447096d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4697753849207649d+0
       v=0.1074429975385679d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6334563241139567d+0
       v=0.1129300086569132d-2
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5974048614181342d-1
       b=0.2029128752777523d+0
       v=0.8436884500901954d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1375760408473636d+0
       b=0.4602621942484054d+0
       v=0.1075255720448885d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3391016526336286d+0
       b=0.5030673999662036d+0
       v=0.1108577236864462d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1271675191439820d+0
       b=0.2817606422442134d+0
       v=0.9566475323783357d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2693120740413512d+0
       b=0.4331561291720157d+0
       v=0.1080663250717391d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1419786452601918d+0
       b=0.6256167358580814d+0
       v=0.1126797131196295d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6709284600738255d-1
       b=0.3798395216859157d+0
       v=0.1022568715358061d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7057738183256172d-1
       b=0.5517505421423520d+0
       v=0.1108960267713108d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2783888477882155d+0
       b=0.6029619156159187d+0
       v=0.1122790653435766d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1979578938917407d+0
       b=0.3589606329589096d+0
       v=0.1032401847117460d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2087307061103274d+0
       b=0.5348666438135476d+0
       v=0.1107249382283854d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4055122137872836d+0
       b=0.5674997546074373d+0
       v=0.1121780048519972d-2
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 1202-point angular grid
!
       subroutine ld1202(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(1202)
! y-coordinates of the gridpoints
       real(8) :: y(1202)
! z-coordinates of the gridpoints
       real(8) :: z(1202)
! weight of the grid point
       real(8) :: w(1202)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1105189233267572d-3
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.9205232738090741d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.9133159786443561d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3712636449657089d-1
       v=0.3690421898017899d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9140060412262223d-1
       v=0.5603990928680660d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1531077852469906d+0
       v=0.6865297629282609d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2180928891660612d+0
       v=0.7720338551145630d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2839874532200175d+0
       v=0.8301545958894795d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3491177600963764d+0
       v=0.8686692550179628d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4121431461444309d+0
       v=0.8927076285846890d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4718993627149127d+0
       v=0.9060820238568219d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5273145452842337d+0
       v=0.9119777254940867d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6209475332444019d+0
       v=0.9128720138604181d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6569722711857291d+0
       v=0.9130714935691735d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6841788309070143d+0
       v=0.9152873784554116d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7012604330123631d+0
       v=0.9187436274321654d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1072382215478166d+0
       v=0.5176977312965694d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2582068959496968d+0
       v=0.7331143682101417d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4172752955306717d+0
       v=0.8463232836379928d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5700366911792503d+0
       v=0.9031122694253992d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9827986018263947d+0
       b=0.1771774022615325d+0
       v=0.6485778453163257d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9624249230326228d+0
       b=0.2475716463426288d+0
       v=0.7435030910982369d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9402007994128811d+0
       b=0.3354616289066489d+0
       v=0.7998527891839054d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9320822040143202d+0
       b=0.3173615246611977d+0
       v=0.8101731497468018d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9043674199393299d+0
       b=0.4090268427085357d+0
       v=0.8483389574594331d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8912407560074747d+0
       b=0.3854291150669224d+0
       v=0.8556299257311812d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8676435628462708d+0
       b=0.4932221184851285d+0
       v=0.8803208679738260d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8581979986041619d+0
       b=0.4785320675922435d+0
       v=0.8811048182425720d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8396753624049856d+0
       b=0.4507422593157064d+0
       v=0.8850282341265444d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8165288564022188d+0
       b=0.5632123020762100d+0
       v=0.9021342299040653d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8015469370783529d+0
       b=0.5434303569693900d+0
       v=0.9010091677105086d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7773563069070351d+0
       b=0.5123518486419871d+0
       v=0.9022692938426915d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7661621213900394d+0
       b=0.6394279634749102d+0
       v=0.9158016174693465d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7553584143533510d+0
       b=0.6269805509024392d+0
       v=0.9131578003189435d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7344305757559503d+0
       b=0.6031161693096310d+0
       v=0.9107813579482705d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7043837184021765d+0
       b=0.5693702498468441d+0
       v=0.9105760258970126d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 1454-point angular grid
!
       subroutine ld1454(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(1454)
! y-coordinates of the gridpoints
       real(8) :: y(1454)
! z-coordinates of the gridpoints
       real(8) :: z(1454)
! weight of the grid point
       real(8) :: w(1454)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.7777160743261247d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.7557646413004701d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3229290663413854d-1
       v=0.2841633806090617d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8036733271462222d-1
       v=0.4374419127053555d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1354289960531653d+0
       v=0.5417174740872172d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1938963861114426d+0
       v=0.6148000891358593d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2537343715011275d+0
       v=0.6664394485800705d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3135251434752570d+0
       v=0.7025039356923220d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3721558339375338d+0
       v=0.7268511789249627d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4286809575195696d+0
       v=0.7422637534208629d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4822510128282994d+0
       v=0.7509545035841214d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5320679333566263d+0
       v=0.7548535057718401d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6172998195394274d+0
       v=0.7554088969774001d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6510679849127481d+0
       v=0.7553147174442808d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6777315251687360d+0
       v=0.7564767653292297d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6963109410648741d+0
       v=0.7587991808518730d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7058935009831749d+0
       v=0.7608261832033027d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9955546194091857d+0
       v=0.4021680447874916d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9734115901794209d+0
       v=0.5804871793945964d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9275693732388626d+0
       v=0.6792151955945159d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8568022422795103d+0
       v=0.7336741211286294d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7623495553719372d+0
       v=0.7581866300989608d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5707522908892223d+0
       b=0.4387028039889501d+0
       v=0.7538257859800743d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5196463388403083d+0
       b=0.3858908414762617d+0
       v=0.7483517247053123d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4646337531215351d+0
       b=0.3301937372343854d+0
       v=0.7371763661112059d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4063901697557691d+0
       b=0.2725423573563777d+0
       v=0.7183448895756934d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3456329466643087d+0
       b=0.2139510237495250d+0
       v=0.6895815529822191d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2831395121050332d+0
       b=0.1555922309786647d+0
       v=0.6480105801792886d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2197682022925330d+0
       b=0.9892878979686097d-1
       v=0.5897558896594636d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1564696098650355d+0
       b=0.4598642910675510d-1
       v=0.5095708849247346d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6027356673721295d+0
       b=0.3376625140173426d+0
       v=0.7536906428909755d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5496032320255096d+0
       b=0.2822301309727988d+0
       v=0.7472505965575118d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4921707755234567d+0
       b=0.2248632342592540d+0
       v=0.7343017132279698d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4309422998598483d+0
       b=0.1666224723456479d+0
       v=0.7130871582177445d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3664108182313672d+0
       b=0.1086964901822169d+0
       v=0.6817022032112776d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2990189057758436d+0
       b=0.5251989784120085d-1
       v=0.6380941145604121d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6268724013144998d+0
       b=0.2297523657550023d+0
       v=0.7550381377920310d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5707324144834607d+0
       b=0.1723080607093800d+0
       v=0.7478646640144802d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5096360901960365d+0
       b=0.1140238465390513d+0
       v=0.7335918720601220d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4438729938312456d+0
       b=0.5611522095882537d-1
       v=0.7110120527658118d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6419978471082389d+0
       b=0.1164174423140873d+0
       v=0.7571363978689501d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5817218061802611d+0
       b=0.5797589531445219d-1
       v=0.7489908329079234d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 1730-point angular grid
!
       subroutine ld1730(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(1730)
! y-coordinates of the gridpoints
       real(8) :: y(1730)
! z-coordinates of the gridpoints
       real(8) :: z(1730)
! weight of the grid point
       real(8) :: w(1730)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.6309049437420976d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.6398287705571748d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.6357185073530720d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2860923126194662d-1
       v=0.2221207162188168d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7142556767711522d-1
       v=0.3475784022286848d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1209199540995559d+0
       v=0.4350742443589804d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1738673106594379d+0
       v=0.4978569136522127d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2284645438467734d+0
       v=0.5435036221998053d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2834807671701512d+0
       v=0.5765913388219542d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3379680145467339d+0
       v=0.6001200359226003d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3911355454819537d+0
       v=0.6162178172717512d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4422860353001403d+0
       v=0.6265218152438485d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4907781568726057d+0
       v=0.6323987160974212d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5360006153211468d+0
       v=0.6350767851540569d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6142105973596603d+0
       v=0.6354362775297107d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6459300387977504d+0
       v=0.6352302462706235d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6718056125089225d+0
       v=0.6358117881417972d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6910888533186254d+0
       v=0.6373101590310117d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7030467416823252d+0
       v=0.6390428961368665d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8354951166354646d-1
       v=0.3186913449946576d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2050143009099486d+0
       v=0.4678028558591711d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3370208290706637d+0
       v=0.5538829697598626d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4689051484233963d+0
       v=0.6044475907190476d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5939400424557334d+0
       v=0.6313575103509012d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1394983311832261d+0
       b=0.4097581162050343d-1
       v=0.4078626431855630d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1967999180485014d+0
       b=0.8851987391293348d-1
       v=0.4759933057812725d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2546183732548967d+0
       b=0.1397680182969819d+0
       v=0.5268151186413440d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3121281074713875d+0
       b=0.1929452542226526d+0
       v=0.5643048560507316d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3685981078502492d+0
       b=0.2467898337061562d+0
       v=0.5914501076613073d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4233760321547856d+0
       b=0.3003104124785409d+0
       v=0.6104561257874195d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4758671236059246d+0
       b=0.3526684328175033d+0
       v=0.6230252860707806d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5255178579796463d+0
       b=0.4031134861145713d+0
       v=0.6305618761760796d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5718025633734589d+0
       b=0.4509426448342351d+0
       v=0.6343092767597889d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2686927772723415d+0
       b=0.4711322502423248d-1
       v=0.5176268945737826d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3306006819904809d+0
       b=0.9784487303942695d-1
       v=0.5564840313313692d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3904906850594983d+0
       b=0.1505395810025273d+0
       v=0.5856426671038980d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4479957951904390d+0
       b=0.2039728156296050d+0
       v=0.6066386925777091d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5027076848919780d+0
       b=0.2571529941121107d+0
       v=0.6208824962234458d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5542087392260217d+0
       b=0.3092191375815670d+0
       v=0.6296314297822907d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6020850887375187d+0
       b=0.3593807506130276d+0
       v=0.6340423756791859d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4019851409179594d+0
       b=0.5063389934378671d-1
       v=0.5829627677107342d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4635614567449800d+0
       b=0.1032422269160612d+0
       v=0.6048693376081110d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5215860931591575d+0
       b=0.1566322094006254d+0
       v=0.6202362317732461d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5758202499099271d+0
       b=0.2098082827491099d+0
       v=0.6299005328403779d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6259893683876795d+0
       b=0.2618824114553391d+0
       v=0.6347722390609353d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5313795124811891d+0
       b=0.5263245019338556d-1
       v=0.6203778981238834d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5893317955931995d+0
       b=0.1061059730982005d+0
       v=0.6308414671239979d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6426246321215801d+0
       b=0.1594171564034221d+0
       v=0.6362706466959498d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6511904367376113d+0
       b=0.5354789536565540d-1
       v=0.6375414170333233d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 2030-point angular grid
!
       subroutine ld2030(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(2030)
! y-coordinates of the gridpoints
       real(8) :: y(2030)
! z-coordinates of the gridpoints
       real(8) :: z(2030)
! weight of the grid point
       real(8) :: w(2030)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.4656031899197431d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.5421549195295507d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2540835336814348d-1
       v=0.1778522133346553d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6399322800504915d-1
       v=0.2811325405682796d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1088269469804125d+0
       v=0.3548896312631459d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1570670798818287d+0
       v=0.4090310897173364d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2071163932282514d+0
       v=0.4493286134169965d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2578914044450844d+0
       v=0.4793728447962723d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3085687558169623d+0
       v=0.5015415319164265d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3584719706267024d+0
       v=0.5175127372677937d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4070135594428709d+0
       v=0.5285522262081019d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4536618626222638d+0
       v=0.5356832703713962d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4979195686463577d+0
       v=0.5397914736175170d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5393075111126999d+0
       v=0.5416899441599930d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6115617676843916d+0
       v=0.5419308476889938d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6414308435160159d+0
       v=0.5416936902030596d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6664099412721607d+0
       v=0.5419544338703164d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6859161771214913d+0
       v=0.5428983656630975d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6993625593503890d+0
       v=0.5442286500098193d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7062393387719380d+0
       v=0.5452250345057301d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7479028168349763d-1
       v=0.2568002497728530d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1848951153969366d+0
       v=0.3827211700292145d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3059529066581305d+0
       v=0.4579491561917824d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4285556101021362d+0
       v=0.5042003969083574d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5468758653496526d+0
       v=0.5312708889976025d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6565821978343439d+0
       v=0.5438401790747117d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1253901572367117d+0
       b=0.3681917226439641d-1
       v=0.3316041873197344d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1775721510383941d+0
       b=0.7982487607213301d-1
       v=0.3899113567153771d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2305693358216114d+0
       b=0.1264640966592335d+0
       v=0.4343343327201309d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2836502845992063d+0
       b=0.1751585683418957d+0
       v=0.4679415262318919d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3361794746232590d+0
       b=0.2247995907632670d+0
       v=0.4930847981631031d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3875979172264824d+0
       b=0.2745299257422246d+0
       v=0.5115031867540091d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4374019316999074d+0
       b=0.3236373482441118d+0
       v=0.5245217148457367d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4851275843340022d+0
       b=0.3714967859436741d+0
       v=0.5332041499895321d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5303391803806868d+0
       b=0.4175353646321745d+0
       v=0.5384583126021542d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5726197380596287d+0
       b=0.4612084406355461d+0
       v=0.5411067210798852d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2431520732564863d+0
       b=0.4258040133043952d-1
       v=0.4259797391468714d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3002096800895869d+0
       b=0.8869424306722721d-1
       v=0.4604931368460021d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3558554457457432d+0
       b=0.1368811706510655d+0
       v=0.4871814878255202d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4097782537048887d+0
       b=0.1860739985015033d+0
       v=0.5072242910074885d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4616337666067458d+0
       b=0.2354235077395853d+0
       v=0.5217069845235350d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5110707008417874d+0
       b=0.2842074921347011d+0
       v=0.5315785966280310d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5577415286163795d+0
       b=0.3317784414984102d+0
       v=0.5376833708758905d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6013060431366950d+0
       b=0.3775299002040700d+0
       v=0.5408032092069521d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3661596767261781d+0
       b=0.4599367887164592d-1
       v=0.4842744917904866d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4237633153506581d+0
       b=0.9404893773654421d-1
       v=0.5048926076188130d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4786328454658452d+0
       b=0.1431377109091971d+0
       v=0.5202607980478373d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5305702076789774d+0
       b=0.1924186388843570d+0
       v=0.5309932388325743d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5793436224231788d+0
       b=0.2411590944775190d+0
       v=0.5377419770895208d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6247069017094747d+0
       b=0.2886871491583605d+0
       v=0.5411696331677717d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4874315552535204d+0
       b=0.4804978774953206d-1
       v=0.5197996293282420d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5427337322059053d+0
       b=0.9716857199366665d-1
       v=0.5311120836622945d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5943493747246700d+0
       b=0.1465205839795055d+0
       v=0.5384309319956951d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6421314033564943d+0
       b=0.1953579449803574d+0
       v=0.5421859504051886d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6020628374713980d+0
       b=0.4916375015738108d-1
       v=0.5390948355046314d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6529222529856881d+0
       b=0.9861621540127005d-1
       v=0.5433312705027845d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 2354-point angular grid
!
       subroutine ld2354(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(2354)
! y-coordinates of the gridpoints
       real(8) :: y(2354)
! z-coordinates of the gridpoints
       real(8) :: z(2354)
! weight of the grid point
       real(8) :: w(2354)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.3922616270665292d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.4703831750854424d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.4678202801282136d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2290024646530589d-1
       v=0.1437832228979900d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5779086652271284d-1
       v=0.2303572493577644d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9863103576375984d-1
       v=0.2933110752447454d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1428155792982185d+0
       v=0.3402905998359838d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1888978116601463d+0
       v=0.3759138466870372d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2359091682970210d+0
       v=0.4030638447899798d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2831228833706171d+0
       v=0.4236591432242211d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3299495857966693d+0
       v=0.4390522656946746d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3758840802660796d+0
       v=0.4502523466626247d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4204751831009480d+0
       v=0.4580577727783541d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4633068518751051d+0
       v=0.4631391616615899d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5039849474507313d+0
       v=0.4660928953698676d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5421265793440747d+0
       v=0.4674751807936953d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6092660230557310d+0
       v=0.4676414903932920d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6374654204984869d+0
       v=0.4674086492347870d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6615136472609892d+0
       v=0.4674928539483207d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6809487285958127d+0
       v=0.4680748979686447d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6952980021665196d+0
       v=0.4690449806389040d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7041245497695400d+0
       v=0.4699877075860818d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6744033088306065d-1
       v=0.2099942281069176d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1678684485334166d+0
       v=0.3172269150712804d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2793559049539613d+0
       v=0.3832051358546523d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3935264218057639d+0
       v=0.4252193818146985d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5052629268232558d+0
       v=0.4513807963755000d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6107905315437531d+0
       v=0.4657797469114178d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1135081039843524d+0
       b=0.3331954884662588d-1
       v=0.2733362800522836d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1612866626099378d+0
       b=0.7247167465436538d-1
       v=0.3235485368463559d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2100786550168205d+0
       b=0.1151539110849745d+0
       v=0.3624908726013453d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2592282009459942d+0
       b=0.1599491097143677d+0
       v=0.3925540070712828d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3081740561320203d+0
       b=0.2058699956028027d+0
       v=0.4156129781116235d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3564289781578164d+0
       b=0.2521624953502911d+0
       v=0.4330644984623263d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4035587288240703d+0
       b=0.2982090785797674d+0
       v=0.4459677725921312d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4491671196373903d+0
       b=0.3434762087235733d+0
       v=0.4551593004456795d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4928854782917489d+0
       b=0.3874831357203437d+0
       v=0.4613341462749918d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5343646791958988d+0
       b=0.4297814821746926d+0
       v=0.4651019618269806d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5732683216530990d+0
       b=0.4699402260943537d+0
       v=0.4670249536100625d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2214131583218986d+0
       b=0.3873602040643895d-1
       v=0.3549555576441708d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2741796504750071d+0
       b=0.8089496256902013d-1
       v=0.3856108245249010d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3259797439149485d+0
       b=0.1251732177620872d+0
       v=0.4098622845756882d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3765441148826891d+0
       b=0.1706260286403185d+0
       v=0.4286328604268950d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4255773574530558d+0
       b=0.2165115147300408d+0
       v=0.4427802198993945d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4727795117058430d+0
       b=0.2622089812225259d+0
       v=0.4530473511488561d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5178546895819012d+0
       b=0.3071721431296201d+0
       v=0.4600805475703138d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5605141192097460d+0
       b=0.3508998998801138d+0
       v=0.4644599059958017d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6004763319352512d+0
       b=0.3929160876166931d+0
       v=0.4667274455712508d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3352842634946949d+0
       b=0.4202563457288019d-1
       v=0.4069360518020356d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3891971629814670d+0
       b=0.8614309758870850d-1
       v=0.4260442819919195d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4409875565542281d+0
       b=0.1314500879380001d+0
       v=0.4408678508029063d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4904893058592484d+0
       b=0.1772189657383859d+0
       v=0.4518748115548597d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5375056138769549d+0
       b=0.2228277110050294d+0
       v=0.4595564875375116d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5818255708669969d+0
       b=0.2677179935014386d+0
       v=0.4643988774315846d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6232334858144959d+0
       b=0.3113675035544165d+0
       v=0.4668827491646946d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4489485354492058d+0
       b=0.4409162378368174d-1
       v=0.4400541823741973d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5015136875933150d+0
       b=0.8939009917748489d-1
       v=0.4514512890193797d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5511300550512623d+0
       b=0.1351806029383365d+0
       v=0.4596198627347549d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5976720409858000d+0
       b=0.1808370355053196d+0
       v=0.4648659016801781d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6409956378989354d+0
       b=0.2257852192301602d+0
       v=0.4675502017157673d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5581222330827514d+0
       b=0.4532173421637160d-1
       v=0.4598494476455523d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6074705984161695d+0
       b=0.9117488031840314d-1
       v=0.4654916955152048d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6532272537379033d+0
       b=0.1369294213140155d+0
       v=0.4684709779505137d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6594761494500487d+0
       b=0.4589901487275583d-1
       v=0.4691445539106986d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 2702-point angular grid
!
       subroutine ld2702(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(2702)
! y-coordinates of the gridpoints
       real(8) :: y(2702)
! z-coordinates of the gridpoints
       real(8) :: z(2702)
! weight of the grid point
       real(8) :: w(2702)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.2998675149888161d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.4077860529495355d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2065562538818703d-1
       v=0.1185349192520667d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5250918173022379d-1
       v=0.1913408643425751d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8993480082038376d-1
       v=0.2452886577209897d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1306023924436019d+0
       v=0.2862408183288702d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1732060388531418d+0
       v=0.3178032258257357d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2168727084820249d+0
       v=0.3422945667633690d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2609528309173586d+0
       v=0.3612790520235922d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3049252927938952d+0
       v=0.3758638229818521d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3483484138084404d+0
       v=0.3868711798859953d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3908321549106406d+0
       v=0.3949429933189938d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4320210071894814d+0
       v=0.4006068107541156d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4715824795890053d+0
       v=0.4043192149672723d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5091984794078453d+0
       v=0.4064947495808078d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5445580145650803d+0
       v=0.4075245619813152d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6072575796841768d+0
       v=0.4076423540893566d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6339484505755803d+0
       v=0.4074280862251555d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6570718257486958d+0
       v=0.4074163756012244d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6762557330090709d+0
       v=0.4077647795071246d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6911161696923790d+0
       v=0.4084517552782530d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7012841911659961d+0
       v=0.4092468459224052d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7064559272410020d+0
       v=0.4097872687240906d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6123554989894765d-1
       v=0.1738986811745028d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1533070348312393d+0
       v=0.2659616045280191d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2563902605244206d+0
       v=0.3240596008171533d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3629346991663361d+0
       v=0.3621195964432943d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4683949968987538d+0
       v=0.3868838330760539d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5694479240657952d+0
       v=0.4018911532693111d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6634465430993955d+0
       v=0.4089929432983252d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1033958573552305d+0
       b=0.3034544009063584d-1
       v=0.2279907527706409d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1473521412414395d+0
       b=0.6618803044247135d-1
       v=0.2715205490578897d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1924552158705967d+0
       b=0.1054431128987715d+0
       v=0.3057917896703976d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2381094362890328d+0
       b=0.1468263551238858d+0
       v=0.3326913052452555d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2838121707936760d+0
       b=0.1894486108187886d+0
       v=0.3537334711890037d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3291323133373415d+0
       b=0.2326374238761579d+0
       v=0.3700567500783129d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3736896978741460d+0
       b=0.2758485808485768d+0
       v=0.3825245372589122d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4171406040760013d+0
       b=0.3186179331996921d+0
       v=0.3918125171518296d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4591677985256915d+0
       b=0.3605329796303794d+0
       v=0.3984720419937579d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4994733831718418d+0
       b=0.4012147253586509d+0
       v=0.4029746003338211d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5377731830445096d+0
       b=0.4403050025570692d+0
       v=0.4057428632156627d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5737917830001331d+0
       b=0.4774565904277483d+0
       v=0.4071719274114857d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2027323586271389d+0
       b=0.3544122504976147d-1
       v=0.2990236950664119d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2516942375187273d+0
       b=0.7418304388646328d-1
       v=0.3262951734212878d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3000227995257181d+0
       b=0.1150502745727186d+0
       v=0.3482634608242413d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3474806691046342d+0
       b=0.1571963371209364d+0
       v=0.3656596681700892d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3938103180359209d+0
       b=0.1999631877247100d+0
       v=0.3791740467794218d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4387519590455703d+0
       b=0.2428073457846535d+0
       v=0.3894034450156905d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4820503960077787d+0
       b=0.2852575132906155d+0
       v=0.3968600245508371d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5234573778475101d+0
       b=0.3268884208674639d+0
       v=0.4019931351420050d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5627318647235282d+0
       b=0.3673033321675939d+0
       v=0.4052108801278599d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5996390607156954d+0
       b=0.4061211551830290d+0
       v=0.4068978613940934d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3084780753791947d+0
       b=0.3860125523100059d-1
       v=0.3454275351319704d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3589988275920223d+0
       b=0.7928938987104867d-1
       v=0.3629963537007920d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4078628415881973d+0
       b=0.1212614643030087d+0
       v=0.3770187233889873d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4549287258889735d+0
       b=0.1638770827382693d+0
       v=0.3878608613694378d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5000278512957279d+0
       b=0.2065965798260176d+0
       v=0.3959065270221274d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5429785044928199d+0
       b=0.2489436378852235d+0
       v=0.4015286975463570d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5835939850491711d+0
       b=0.2904811368946891d+0
       v=0.4050866785614717d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6216870353444856d+0
       b=0.3307941957666609d+0
       v=0.4069320185051913d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4151104662709091d+0
       b=0.4064829146052554d-1
       v=0.3760120964062763d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4649804275009218d+0
       b=0.8258424547294755d-1
       v=0.3870969564418064d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5124695757009662d+0
       b=0.1251841962027289d+0
       v=0.3955287790534055d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5574711100606224d+0
       b=0.1679107505976331d+0
       v=0.4015361911302668d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5998597333287227d+0
       b=0.2102805057358715d+0
       v=0.4053836986719548d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6395007148516600d+0
       b=0.2518418087774107d+0
       v=0.4073578673299117d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5188456224746252d+0
       b=0.4194321676077518d-1
       v=0.3954628379231406d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5664190707942778d+0
       b=0.8457661551921499d-1
       v=0.4017645508847530d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6110464353283153d+0
       b=0.1273652932519396d+0
       v=0.4059030348651293d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6526430302051563d+0
       b=0.1698173239076354d+0
       v=0.4080565809484880d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6167551880377548d+0
       b=0.4266398851548864d-1
       v=0.4063018753664651d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6607195418355383d+0
       b=0.8551925814238349d-1
       v=0.4087191292799671d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 3074-point angular grid
!
       subroutine ld3074(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(3074)
! y-coordinates of the gridpoints
       real(8) :: y(3074)
! z-coordinates of the gridpoints
       real(8) :: z(3074)
! weight of the grid point
       real(8) :: w(3074)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.2599095953754734d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3603134089687541d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3586067974412447d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1886108518723392d-1
       v=0.9831528474385880d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4800217244625303d-1
       v=0.1605023107954450d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8244922058397242d-1
       v=0.2072200131464099d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1200408362484023d+0
       v=0.2431297618814187d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1595773530809965d+0
       v=0.2711819064496707d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2002635973434064d+0
       v=0.2932762038321116d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2415127590139982d+0
       v=0.3107032514197368d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2828584158458477d+0
       v=0.3243808058921213d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3239091015338138d+0
       v=0.3349899091374030d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3643225097962194d+0
       v=0.3430580688505218d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4037897083691802d+0
       v=0.3490124109290343d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4420247515194127d+0
       v=0.3532148948561955d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4787572538464938d+0
       v=0.3559862669062833d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5137265251275234d+0
       v=0.3576224317551411d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5466764056654611d+0
       v=0.3584050533086076d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6054859420813535d+0
       v=0.3584903581373224d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6308106701764562d+0
       v=0.3582991879040586d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6530369230179584d+0
       v=0.3582371187963125d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6718609524611158d+0
       v=0.3584353631122350d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6869676499894013d+0
       v=0.3589120166517785d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6980467077240748d+0
       v=0.3595445704531601d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7048241721250522d+0
       v=0.3600943557111074d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5591105222058232d-1
       v=0.1456447096742039d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1407384078513916d+0
       v=0.2252370188283782d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2364035438976309d+0
       v=0.2766135443474897d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3360602737818170d+0
       v=0.3110729491500851d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4356292630054665d+0
       v=0.3342506712303391d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5321569415256174d+0
       v=0.3491981834026860d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6232956305040554d+0
       v=0.3576003604348932d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9469870086838469d-1
       b=0.2778748387309470d-1
       v=0.1921921305788564d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1353170300568141d+0
       b=0.6076569878628364d-1
       v=0.2301458216495632d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1771679481726077d+0
       b=0.9703072762711040d-1
       v=0.2604248549522893d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2197066664231751d+0
       b=0.1354112458524762d+0
       v=0.2845275425870697d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2624783557374927d+0
       b=0.1750996479744100d+0
       v=0.3036870897974840d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3050969521214442d+0
       b=0.2154896907449802d+0
       v=0.3188414832298066d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3472252637196021d+0
       b=0.2560954625740152d+0
       v=0.3307046414722089d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3885610219026360d+0
       b=0.2965070050624096d+0
       v=0.3398330969031360d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4288273776062765d+0
       b=0.3363641488734497d+0
       v=0.3466757899705373d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4677662471302948d+0
       b=0.3753400029836788d+0
       v=0.3516095923230054d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5051333589553359d+0
       b=0.4131297522144286d+0
       v=0.3549645184048486d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5406942145810492d+0
       b=0.4494423776081795d+0
       v=0.3570415969441392d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5742204122576457d+0
       b=0.4839938958841502d+0
       v=0.3581251798496118d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1865407027225188d+0
       b=0.3259144851070796d-1
       v=0.2543491329913348d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2321186453689432d+0
       b=0.6835679505297343d-1
       v=0.2786711051330776d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2773159142523882d+0
       b=0.1062284864451989d+0
       v=0.2985552361083679d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3219200192237254d+0
       b=0.1454404409323047d+0
       v=0.3145867929154039d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3657032593944029d+0
       b=0.1854018282582510d+0
       v=0.3273290662067609d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4084376778363622d+0
       b=0.2256297412014750d+0
       v=0.3372705511943501d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4499004945751427d+0
       b=0.2657104425000896d+0
       v=0.3448274437851510d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4898758141326335d+0
       b=0.3052755487631557d+0
       v=0.3503592783048583d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5281547442266309d+0
       b=0.3439863920645423d+0
       v=0.3541854792663162d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5645346989813992d+0
       b=0.3815229456121914d+0
       v=0.3565995517909428d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5988181252159848d+0
       b=0.4175752420966734d+0
       v=0.3578802078302898d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2850425424471603d+0
       b=0.3562149509862536d-1
       v=0.2958644592860982d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3324619433027876d+0
       b=0.7330318886871096d-1
       v=0.3119548129116835d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3785848333076282d+0
       b=0.1123226296008472d+0
       v=0.3250745225005984d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4232891028562115d+0
       b=0.1521084193337708d+0
       v=0.3355153415935208d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4664287050829722d+0
       b=0.1921844459223610d+0
       v=0.3435847568549328d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5078458493735726d+0
       b=0.2321360989678303d+0
       v=0.3495786831622488d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5473779816204180d+0
       b=0.2715886486360520d+0
       v=0.3537767805534621d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5848617133811376d+0
       b=0.3101924707571355d+0
       v=0.3564459815421428d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6201348281584888d+0
       b=0.3476121052890973d+0
       v=0.3578464061225468d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3852191185387871d+0
       b=0.3763224880035108d-1
       v=0.3239748762836212d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4325025061073423d+0
       b=0.7659581935637135d-1
       v=0.3345491784174287d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4778486229734490d+0
       b=0.1163381306083900d+0
       v=0.3429126177301782d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5211663693009000d+0
       b=0.1563890598752899d+0
       v=0.3492420343097421d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5623469504853703d+0
       b=0.1963320810149200d+0
       v=0.3537399050235257d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6012718188659246d+0
       b=0.2357847407258738d+0
       v=0.3566209152659172d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6378179206390117d+0
       b=0.2743846121244060d+0
       v=0.3581084321919782d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4836936460214534d+0
       b=0.3895902610739024d-1
       v=0.3426522117591512d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5293792562683797d+0
       b=0.7871246819312640d-1
       v=0.3491848770121379d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5726281253100033d+0
       b=0.1187963808202981d+0
       v=0.3539318235231476d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6133658776169068d+0
       b=0.1587914708061787d+0
       v=0.3570231438458694d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6515085491865307d+0
       b=0.1983058575227646d+0
       v=0.3586207335051714d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5778692716064976d+0
       b=0.3977209689791542d-1
       v=0.3541196205164025d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6207904288086192d+0
       b=0.7990157592981152d-1
       v=0.3574296911573953d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6608688171046802d+0
       b=0.1199671308754309d+0
       v=0.3591993279818963d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6656263089489130d+0
       b=0.4015955957805969d-1
       v=0.3595855034661997d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 3470-point angular grid
!
       subroutine ld3470(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(3470)
! y-coordinates of the gridpoints
       real(8) :: y(3470)
! z-coordinates of the gridpoints
       real(8) :: z(3470)
! weight of the grid point
       real(8) :: w(3470)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.2040382730826330d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.3178149703889544d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1721420832906233d-1
       v=0.8288115128076110d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4408875374981770d-1
       v=0.1360883192522954d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7594680813878681d-1
       v=0.1766854454542662d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1108335359204799d+0
       v=0.2083153161230153d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1476517054388567d+0
       v=0.2333279544657158d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1856731870860615d+0
       v=0.2532809539930247d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2243634099428821d+0
       v=0.2692472184211158d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2633006881662727d+0
       v=0.2819949946811885d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3021340904916283d+0
       v=0.2920953593973030d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3405594048030089d+0
       v=0.2999889782948352d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3783044434007372d+0
       v=0.3060292120496902d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4151194767407910d+0
       v=0.3105109167522192d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4507705766443257d+0
       v=0.3136902387550312d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4850346056573187d+0
       v=0.3157984652454632d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5176950817792470d+0
       v=0.3170516518425422d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5485384240820989d+0
       v=0.3176568425633755d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6039117238943308d+0
       v=0.3177198411207062d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6279956655573113d+0
       v=0.3175519492394733d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6493636169568952d+0
       v=0.3174654952634756d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6677644117704504d+0
       v=0.3175676415467654d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6829368572115624d+0
       v=0.3178923417835410d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6946195818184121d+0
       v=0.3183788287531909d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7025711542057026d+0
       v=0.3188755151918807d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7066004767140119d+0
       v=0.3191916889313849d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5132537689946062d-1
       v=0.1231779611744508d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1297994661331225d+0
       v=0.1924661373839880d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2188852049401307d+0
       v=0.2380881867403424d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3123174824903457d+0
       v=0.2693100663037885d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4064037620738195d+0
       v=0.2908673382834366d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4984958396944782d+0
       v=0.3053914619381535d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5864975046021365d+0
       v=0.3143916684147777d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6686711634580175d+0
       v=0.3187042244055363d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8715738780835950d-1
       b=0.2557175233367578d-1
       v=0.1635219535869790d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1248383123134007d+0
       b=0.5604823383376681d-1
       v=0.1968109917696070d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1638062693383378d+0
       b=0.8968568601900765d-1
       v=0.2236754342249974d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2035586203373176d+0
       b=0.1254086651976279d+0
       v=0.2453186687017181d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2436798975293774d+0
       b=0.1624780150162012d+0
       v=0.2627551791580541d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2838207507773806d+0
       b=0.2003422342683208d+0
       v=0.2767654860152220d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3236787502217692d+0
       b=0.2385628026255263d+0
       v=0.2879467027765895d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3629849554840691d+0
       b=0.2767731148783578d+0
       v=0.2967639918918702d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4014948081992087d+0
       b=0.3146542308245309d+0
       v=0.3035900684660351d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4389818379260225d+0
       b=0.3519196415895088d+0
       v=0.3087338237298308d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4752331143674377d+0
       b=0.3883050984023654d+0
       v=0.3124608838860167d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5100457318374018d+0
       b=0.4235613423908649d+0
       v=0.3150084294226743d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5432238388954868d+0
       b=0.4574484717196220d+0
       v=0.3165958398598402d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5745758685072442d+0
       b=0.4897311639255524d+0
       v=0.3174320440957372d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1723981437592809d+0
       b=0.3010630597881105d-1
       v=0.2182188909812599d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2149553257844597d+0
       b=0.6326031554204694d-1
       v=0.2399727933921445d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2573256081247422d+0
       b=0.9848566980258631d-1
       v=0.2579796133514652d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2993163751238106d+0
       b=0.1350835952384266d+0
       v=0.2727114052623535d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3407238005148000d+0
       b=0.1725184055442181d+0
       v=0.2846327656281355d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3813454978483264d+0
       b=0.2103559279730725d+0
       v=0.2941491102051334d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4209848104423343d+0
       b=0.2482278774554860d+0
       v=0.3016049492136107d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4594519699996300d+0
       b=0.2858099509982883d+0
       v=0.3072949726175648d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4965640166185930d+0
       b=0.3228075659915428d+0
       v=0.3114768142886460d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5321441655571562d+0
       b=0.3589459907204151d+0
       v=0.3143823673666223d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5660208438582166d+0
       b=0.3939630088864310d+0
       v=0.3162269764661535d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5980264315964364d+0
       b=0.4276029922949089d+0
       v=0.3172164663759821d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2644215852350733d+0
       b=0.3300939429072552d-1
       v=0.2554575398967435d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3090113743443063d+0
       b=0.6803887650078501d-1
       v=0.2701704069135677d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3525871079197808d+0
       b=0.1044326136206709d+0
       v=0.2823693413468940d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3950418005354029d+0
       b=0.1416751597517679d+0
       v=0.2922898463214289d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4362475663430163d+0
       b=0.1793408610504821d+0
       v=0.3001829062162428d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4760661812145854d+0
       b=0.2170630750175722d+0
       v=0.3062890864542953d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5143551042512103d+0
       b=0.2545145157815807d+0
       v=0.3108328279264746d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5509709026935597d+0
       b=0.2913940101706601d+0
       v=0.3140243146201245d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5857711030329428d+0
       b=0.3274169910910705d+0
       v=0.3160638030977130d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6186149917404392d+0
       b=0.3623081329317265d+0
       v=0.3171462882206275d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3586894569557064d+0
       b=0.3497354386450040d-1
       v=0.2812388416031796d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4035266610019441d+0
       b=0.7129736739757095d-1
       v=0.2912137500288045d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4467775312332510d+0
       b=0.1084758620193165d+0
       v=0.2993241256502206d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4883638346608543d+0
       b=0.1460915689241772d+0
       v=0.3057101738983822d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5281908348434601d+0
       b=0.1837790832369980d+0
       v=0.3105319326251432d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5661542687149311d+0
       b=0.2212075390874021d+0
       v=0.3139565514428167d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6021450102031452d+0
       b=0.2580682841160985d+0
       v=0.3161543006806366d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6360520783610050d+0
       b=0.2940656362094121d+0
       v=0.3172985960613294d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4521611065087196d+0
       b=0.3631055365867002d-1
       v=0.2989400336901431d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4959365651560963d+0
       b=0.7348318468484350d-1
       v=0.3054555883947677d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5376815804038283d+0
       b=0.1111087643812648d+0
       v=0.3104764960807702d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5773314480243768d+0
       b=0.1488226085145408d+0
       v=0.3141015825977616d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6148113245575056d+0
       b=0.1862892274135151d+0
       v=0.3164520621159896d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6500407462842380d+0
       b=0.2231909701714456d+0
       v=0.3176652305912204d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5425151448707213d+0
       b=0.3718201306118944d-1
       v=0.3105097161023939d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5841860556907931d+0
       b=0.7483616335067346d-1
       v=0.3143014117890550d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6234632186851500d+0
       b=0.1125990834266120d+0
       v=0.3168172866287200d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6602934551848843d+0
       b=0.1501303813157619d+0
       v=0.3181401865570968d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6278573968375105d+0
       b=0.3767559930245720d-1
       v=0.3170663659156037d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6665611711264577d+0
       b=0.7548443301360158d-1
       v=0.3185447944625510d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 3890-point angular grid
!
       subroutine ld3890(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(3890)
! y-coordinates of the gridpoints
       real(8) :: y(3890)
! z-coordinates of the gridpoints
       real(8) :: z(3890)
! weight of the grid point
       real(8) :: w(3890)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1807395252196920d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2848008782238827d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2836065837530581d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1587876419858352d-1
       v=0.7013149266673816d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4069193593751206d-1
       v=0.1162798021956766d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7025888115257997d-1
       v=0.1518728583972105d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1027495450028704d+0
       v=0.1798796108216934d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1371457730893426d+0
       v=0.2022593385972785d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1727758532671953d+0
       v=0.2203093105575464d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2091492038929037d+0
       v=0.2349294234299855d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2458813281751915d+0
       v=0.2467682058747003d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2826545859450066d+0
       v=0.2563092683572224d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3191957291799622d+0
       v=0.2639253896763318d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3552621469299578d+0
       v=0.2699137479265108d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3906329503406230d+0
       v=0.2745196420166739d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4251028614093031d+0
       v=0.2779529197397593d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4584777520111870d+0
       v=0.2803996086684265d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4905711358710193d+0
       v=0.2820302356715842d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5212011669847385d+0
       v=0.2830056747491068d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5501878488737995d+0
       v=0.2834808950776839d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6025037877479342d+0
       v=0.2835282339078929d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6254572689549016d+0
       v=0.2833819267065800d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6460107179528248d+0
       v=0.2832858336906784d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6639541138154251d+0
       v=0.2833268235451244d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6790688515667495d+0
       v=0.2835432677029253d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6911338580371512d+0
       v=0.2839091722743049d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6999385956126490d+0
       v=0.2843308178875841d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7053037748656896d+0
       v=0.2846703550533846d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4732224387180115d-1
       v=0.1051193406971900d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1202100529326803d+0
       v=0.1657871838796974d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2034304820664855d+0
       v=0.2064648113714232d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2912285643573002d+0
       v=0.2347942745819741d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3802361792726768d+0
       v=0.2547775326597726d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4680598511056146d+0
       v=0.2686876684847025d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5528151052155599d+0
       v=0.2778665755515867d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6329386307803041d+0
       v=0.2830996616782929d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8056516651369069d-1
       b=0.2363454684003124d-1
       v=0.1403063340168372d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1156476077139389d+0
       b=0.5191291632545936d-1
       v=0.1696504125939477d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1520473382760421d+0
       b=0.8322715736994519d-1
       v=0.1935787242745390d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1892986699745931d+0
       b=0.1165855667993712d+0
       v=0.2130614510521968d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2270194446777792d+0
       b=0.1513077167409504d+0
       v=0.2289381265931048d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2648908185093273d+0
       b=0.1868882025807859d+0
       v=0.2418630292816186d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3026389259574136d+0
       b=0.2229277629776224d+0
       v=0.2523400495631193d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3400220296151384d+0
       b=0.2590951840746235d+0
       v=0.2607623973449605d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3768217953335510d+0
       b=0.2951047291750847d+0
       v=0.2674441032689209d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4128372900921884d+0
       b=0.3307019714169930d+0
       v=0.2726432360343356d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4478807131815630d+0
       b=0.3656544101087634d+0
       v=0.2765787685924545d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4817742034089257d+0
       b=0.3997448951939695d+0
       v=0.2794428690642224d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5143472814653344d+0
       b=0.4327667110812024d+0
       v=0.2814099002062895d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5454346213905650d+0
       b=0.4645196123532293d+0
       v=0.2826429531578994d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5748739313170252d+0
       b=0.4948063555703345d+0
       v=0.2832983542550884d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1599598738286342d+0
       b=0.2792357590048985d-1
       v=0.1886695565284976d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1998097412500951d+0
       b=0.5877141038139065d-1
       v=0.2081867882748234d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2396228952566202d+0
       b=0.9164573914691377d-1
       v=0.2245148680600796d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2792228341097746d+0
       b=0.1259049641962687d+0
       v=0.2380370491511872d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3184251107546741d+0
       b=0.1610594823400863d+0
       v=0.2491398041852455d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3570481164426244d+0
       b=0.1967151653460898d+0
       v=0.2581632405881230d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3949164710492144d+0
       b=0.2325404606175168d+0
       v=0.2653965506227417d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4318617293970503d+0
       b=0.2682461141151439d+0
       v=0.2710857216747087d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4677221009931678d+0
       b=0.3035720116011973d+0
       v=0.2754434093903659d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5023417939270955d+0
       b=0.3382781859197439d+0
       v=0.2786579932519380d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5355701836636128d+0
       b=0.3721383065625942d+0
       v=0.2809011080679474d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5672608451328771d+0
       b=0.4049346360466055d+0
       v=0.2823336184560987d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5972704202540162d+0
       b=0.4364538098633802d+0
       v=0.2831101175806309d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2461687022333596d+0
       b=0.3070423166833368d-1
       v=0.2221679970354546d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2881774566286831d+0
       b=0.6338034669281885d-1
       v=0.2356185734270703d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3293963604116978d+0
       b=0.9742862487067941d-1
       v=0.2469228344805590d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3697303822241377d+0
       b=0.1323799532282290d+0
       v=0.2562726348642046d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4090663023135127d+0
       b=0.1678497018129336d+0
       v=0.2638756726753028d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4472819355411712d+0
       b=0.2035095105326114d+0
       v=0.2699311157390862d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4842513377231437d+0
       b=0.2390692566672091d+0
       v=0.2746233268403837d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5198477629962928d+0
       b=0.2742649818076149d+0
       v=0.2781225674454771d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5539453011883145d+0
       b=0.3088503806580094d+0
       v=0.2805881254045684d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5864196762401251d+0
       b=0.3425904245906614d+0
       v=0.2821719877004913d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6171484466668390d+0
       b=0.3752562294789468d+0
       v=0.2830222502333124d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3350337830565727d+0
       b=0.3261589934634747d-1
       v=0.2457995956744870d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3775773224758284d+0
       b=0.6658438928081572d-1
       v=0.2551474407503706d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4188155229848973d+0
       b=0.1014565797157954d+0
       v=0.2629065335195311d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4586805892009344d+0
       b=0.1368573320843822d+0
       v=0.2691900449925075d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4970895714224235d+0
       b=0.1724614851951608d+0
       v=0.2741275485754276d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5339505133960747d+0
       b=0.2079779381416412d+0
       v=0.2778530970122595d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5691665792531440d+0
       b=0.2431385788322288d+0
       v=0.2805010567646741d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6026387682680377d+0
       b=0.2776901883049853d+0
       v=0.2822055834031040d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6342676150163307d+0
       b=0.3113881356386632d+0
       v=0.2831016901243473d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4237951119537067d+0
       b=0.3394877848664351d-1
       v=0.2624474901131803d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4656918683234929d+0
       b=0.6880219556291447d-1
       v=0.2688034163039377d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5058857069185980d+0
       b=0.1041946859721635d+0
       v=0.2738932751287636d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5443204666713996d+0
       b=0.1398039738736393d+0
       v=0.2777944791242523d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5809298813759742d+0
       b=0.1753373381196155d+0
       v=0.2806011661660987d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6156416039447128d+0
       b=0.2105215793514010d+0
       v=0.2824181456597460d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6483801351066604d+0
       b=0.2450953312157051d+0
       v=0.2833585216577828d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5103616577251688d+0
       b=0.3485560643800719d-1
       v=0.2738165236962878d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5506738792580681d+0
       b=0.7026308631512033d-1
       v=0.2778365208203180d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5889573040995292d+0
       b=0.1059035061296403d+0
       v=0.2807852940418966d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6251641589516930d+0
       b=0.1414823925236026d+0
       v=0.2827245949674705d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6592414921570178d+0
       b=0.1767207908214530d+0
       v=0.2837342344829828d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5930314017533384d+0
       b=0.3542189339561672d-1
       v=0.2809233907610981d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6309812253390175d+0
       b=0.7109574040369549d-1
       v=0.2829930809742694d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6666296011353230d+0
       b=0.1067259792282730d+0
       v=0.2841097874111479d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6703715271049922d+0
       b=0.3569455268820809d-1
       v=0.2843455206008783d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 4334-point angular grid
!
       subroutine ld4334(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(4334)
! y-coordinates of the gridpoints
       real(8) :: y(4334)
! z-coordinates of the gridpoints
       real(8) :: z(4334)
! weight of the grid point
       real(8) :: w(4334)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.1449063022537883d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2546377329828424d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1462896151831013d-1
       v=0.6018432961087496d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3769840812493139d-1
       v=0.1002286583263673d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6524701904096891d-1
       v=0.1315222931028093d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9560543416134648d-1
       v=0.1564213746876724d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1278335898929198d+0
       v=0.1765118841507736d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1613096104466031d+0
       v=0.1928737099311080d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1955806225745371d+0
       v=0.2062658534263270d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2302935218498028d+0
       v=0.2172395445953787d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2651584344113027d+0
       v=0.2262076188876047d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2999276825183209d+0
       v=0.2334885699462397d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3343828669718798d+0
       v=0.2393355273179203d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3683265013750518d+0
       v=0.2439559200468863d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4015763206518108d+0
       v=0.2475251866060002d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4339612026399770d+0
       v=0.2501965558158773d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4653180651114582d+0
       v=0.2521081407925925d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4954893331080803d+0
       v=0.2533881002388081d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5243207068924930d+0
       v=0.2541582900848261d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5516590479041704d+0
       v=0.2545365737525860d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6012371927804176d+0
       v=0.2545726993066799d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6231574466449819d+0
       v=0.2544456197465555d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6429416514181271d+0
       v=0.2543481596881064d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6604124272943595d+0
       v=0.2543506451429194d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6753851470408250d+0
       v=0.2544905675493763d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6876717970626160d+0
       v=0.2547611407344429d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6970895061319234d+0
       v=0.2551060375448869d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7034746912553310d+0
       v=0.2554291933816039d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7067017217542295d+0
       v=0.2556255710686343d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4382223501131123d-1
       v=0.9041339695118195d-4
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1117474077400006d+0
       v=0.1438426330079022d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1897153252911440d+0
       v=0.1802523089820518d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2724023009910331d+0
       v=0.2060052290565496d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3567163308709902d+0
       v=0.2245002248967466d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4404784483028087d+0
       v=0.2377059847731150d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5219833154161411d+0
       v=0.2468118955882525d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5998179868977553d+0
       v=0.2525410872966528d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6727803154548222d+0
       v=0.2553101409933397d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7476563943166086d-1
       b=0.2193168509461185d-1
       v=0.1212879733668632d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1075341482001416d+0
       b=0.4826419281533887d-1
       v=0.1472872881270931d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1416344885203259d+0
       b=0.7751191883575742d-1
       v=0.1686846601010828d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1766325315388586d+0
       b=0.1087558139247680d+0
       v=0.1862698414660208d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2121744174481514d+0
       b=0.1413661374253096d+0
       v=0.2007430956991861d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2479669443408145d+0
       b=0.1748768214258880d+0
       v=0.2126568125394796d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2837600452294113d+0
       b=0.2089216406612073d+0
       v=0.2224394603372113d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3193344933193984d+0
       b=0.2431987685545972d+0
       v=0.2304264522673135d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3544935442438745d+0
       b=0.2774497054377770d+0
       v=0.2368854288424087d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3890571932288154d+0
       b=0.3114460356156915d+0
       v=0.2420352089461772d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4228581214259090d+0
       b=0.3449806851913012d+0
       v=0.2460597113081295d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4557387211304052d+0
       b=0.3778618641248256d+0
       v=0.2491181912257687d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4875487950541643d+0
       b=0.4099086391698978d+0
       v=0.2513528194205857d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5181436529962997d+0
       b=0.4409474925853973d+0
       v=0.2528943096693220d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5473824095600661d+0
       b=0.4708094517711291d+0
       v=0.2538660368488136d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5751263398976174d+0
       b=0.4993275140354637d+0
       v=0.2543868648299022d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1489515746840028d+0
       b=0.2599381993267017d-1
       v=0.1642595537825183d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1863656444351767d+0
       b=0.5479286532462190d-1
       v=0.1818246659849308d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2238602880356348d+0
       b=0.8556763251425254d-1
       v=0.1966565649492420d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2612723375728160d+0
       b=0.1177257802267011d+0
       v=0.2090677905657991d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2984332990206190d+0
       b=0.1508168456192700d+0
       v=0.2193820409510504d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3351786584663333d+0
       b=0.1844801892177727d+0
       v=0.2278870827661928d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3713505522209120d+0
       b=0.2184145236087598d+0
       v=0.2348283192282090d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4067981098954663d+0
       b=0.2523590641486229d+0
       v=0.2404139755581477d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4413769993687534d+0
       b=0.2860812976901373d+0
       v=0.2448227407760734d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4749487182516394d+0
       b=0.3193686757808996d+0
       v=0.2482110455592573d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5073798105075426d+0
       b=0.3520226949547602d+0
       v=0.2507192397774103d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5385410448878654d+0
       b=0.3838544395667890d+0
       v=0.2524765968534880d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5683065353670530d+0
       b=0.4146810037640963d+0
       v=0.2536052388539425d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5965527620663510d+0
       b=0.4443224094681121d+0
       v=0.2542230588033068d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2299227700856157d+0
       b=0.2865757664057584d-1
       v=0.1944817013047896d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2695752998553267d+0
       b=0.5923421684485993d-1
       v=0.2067862362746635d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3086178716611389d+0
       b=0.9117817776057715d-1
       v=0.2172440734649114d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3469649871659077d+0
       b=0.1240593814082605d+0
       v=0.2260125991723423d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3845153566319655d+0
       b=0.1575272058259175d+0
       v=0.2332655008689523d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4211600033403215d+0
       b=0.1912845163525413d+0
       v=0.2391699681532458d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4567867834329882d+0
       b=0.2250710177858171d+0
       v=0.2438801528273928d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4912829319232061d+0
       b=0.2586521303440910d+0
       v=0.2475370504260665d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5245364793303812d+0
       b=0.2918112242865407d+0
       v=0.2502707235640574d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5564369788915756d+0
       b=0.3243439239067890d+0
       v=0.2522031701054241d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5868757697775287d+0
       b=0.3560536787835351d+0
       v=0.2534511269978784d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6157458853519617d+0
       b=0.3867480821242581d+0
       v=0.2541284914955151d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3138461110672113d+0
       b=0.3051374637507278d-1
       v=0.2161509250688394d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3542495872050569d+0
       b=0.6237111233730755d-1
       v=0.2248778513437852d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3935751553120181d+0
       b=0.9516223952401907d-1
       v=0.2322388803404617d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4317634668111147d+0
       b=0.1285467341508517d+0
       v=0.2383265471001355d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4687413842250821d+0
       b=0.1622318931656033d+0
       v=0.2432476675019525d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5044274237060283d+0
       b=0.1959581153836453d+0
       v=0.2471122223750674d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5387354077925727d+0
       b=0.2294888081183837d+0
       v=0.2500291752486870d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5715768898356105d+0
       b=0.2626031152713945d+0
       v=0.2521055942764682d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6028627200136111d+0
       b=0.2950904075286713d+0
       v=0.2534472785575503d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6325039812653463d+0
       b=0.3267458451113286d+0
       v=0.2541599713080121d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3981986708423407d+0
       b=0.3183291458749821d-1
       v=0.2317380975862936d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4382791182133300d+0
       b=0.6459548193880908d-1
       v=0.2378550733719775d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4769233057218166d+0
       b=0.9795757037087952d-1
       v=0.2428884456739118d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5140823911194238d+0
       b=0.1316307235126655d+0
       v=0.2469002655757292d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5496977833862983d+0
       b=0.1653556486358704d+0
       v=0.2499657574265851d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5837047306512727d+0
       b=0.1988931724126510d+0
       v=0.2521676168486082d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6160349566926879d+0
       b=0.2320174581438950d+0
       v=0.2535935662645334d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6466185353209440d+0
       b=0.2645106562168662d+0
       v=0.2543356743363214d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4810835158795404d+0
       b=0.3275917807743992d-1
       v=0.2427353285201535d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5199925041324341d+0
       b=0.6612546183967181d-1
       v=0.2468258039744386d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5571717692207494d+0
       b=0.9981498331474143d-1
       v=0.2500060956440310d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5925789250836378d+0
       b=0.1335687001410374d+0
       v=0.2523238365420979d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6261658523859670d+0
       b=0.1671444402896463d+0
       v=0.2538399260252846d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6578811126669331d+0
       b=0.2003106382156076d+0
       v=0.2546255927268069d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5609624612998100d+0
       b=0.3337500940231335d-1
       v=0.2500583360048449d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5979959659984670d+0
       b=0.6708750335901803d-1
       v=0.2524777638260203d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6330523711054002d+0
       b=0.1008792126424850d+0
       v=0.2540951193860656d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6660960998103972d+0
       b=0.1345050343171794d+0
       v=0.2549524085027472d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6365384364585819d+0
       b=0.3372799460737052d-1
       v=0.2542569507009158d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6710994302899275d+0
       b=0.6755249309678028d-1
       v=0.2552114127580376d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 4802-point angular grid
!
       subroutine ld4802(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(4802)
! y-coordinates of the gridpoints
       real(8) :: y(4802)
! z-coordinates of the gridpoints
       real(8) :: z(4802)
! weight of the grid point
       real(8) :: w(4802)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.9687521879420705d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2307897895367918d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2297310852498558d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2335728608887064d-1
       v=0.7386265944001919d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4352987836550653d-1
       v=0.8257977698542210d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6439200521088801d-1
       v=0.9706044762057630d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9003943631993181d-1
       v=0.1302393847117003d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1196706615548473d+0
       v=0.1541957004600968d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1511715412838134d+0
       v=0.1704459770092199d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1835982828503801d+0
       v=0.1827374890942906d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2165081259155405d+0
       v=0.1926360817436107d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2496208720417563d+0
       v=0.2008010239494833d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2827200673567900d+0
       v=0.2075635983209175d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3156190823994346d+0
       v=0.2131306638690909d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3481476793749115d+0
       v=0.2176562329937335d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3801466086947226d+0
       v=0.2212682262991018d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4114652119634011d+0
       v=0.2240799515668565d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4419598786519751d+0
       v=0.2261959816187525d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4714925949329543d+0
       v=0.2277156368808855d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4999293972879466d+0
       v=0.2287351772128336d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5271387221431248d+0
       v=0.2293490814084085d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5529896780837761d+0
       v=0.2296505312376273d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6000856099481712d+0
       v=0.2296793832318756d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6210562192785175d+0
       v=0.2295785443842974d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6401165879934240d+0
       v=0.2295017931529102d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6571144029244334d+0
       v=0.2295059638184868d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6718910821718863d+0
       v=0.2296232343237362d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6842845591099010d+0
       v=0.2298530178740771d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6941353476269816d+0
       v=0.2301579790280501d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7012965242212991d+0
       v=0.2304690404996513d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7056471428242644d+0
       v=0.2307027995907102d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4595557643585895d-1
       v=0.9312274696671092d-4
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1049316742435023d+0
       v=0.1199919385876926d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1773548879549274d+0
       v=0.1598039138877690d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2559071411236127d+0
       v=0.1822253763574900d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3358156837985898d+0
       v=0.1988579593655040d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4155835743763893d+0
       v=0.2112620102533307d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4937894296167472d+0
       v=0.2201594887699007d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5691569694793316d+0
       v=0.2261622590895036d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6405840854894251d+0
       v=0.2296458453435705d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7345133894143348d-1
       b=0.2177844081486067d-1
       v=0.1006006990267000d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1009859834044931d+0
       b=0.4590362185775188d-1
       v=0.1227676689635876d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1324289619748758d+0
       b=0.7255063095690877d-1
       v=0.1467864280270117d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1654272109607127d+0
       b=0.1017825451960684d+0
       v=0.1644178912101232d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1990767186776461d+0
       b=0.1325652320980364d+0
       v=0.1777664890718961d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2330125945523278d+0
       b=0.1642765374496765d+0
       v=0.1884825664516690d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2670080611108287d+0
       b=0.1965360374337889d+0
       v=0.1973269246453848d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3008753376294316d+0
       b=0.2290726770542238d+0
       v=0.2046767775855328d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3344475596167860d+0
       b=0.2616645495370823d+0
       v=0.2107600125918040d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3675709724070786d+0
       b=0.2941150728843141d+0
       v=0.2157416362266829d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4001000887587812d+0
       b=0.3262440400919066d+0
       v=0.2197557816920721d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4318956350436028d+0
       b=0.3578835350611916d+0
       v=0.2229192611835437d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4628239056795531d+0
       b=0.3888751854043678d+0
       v=0.2253385110212775d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4927563229773636d+0
       b=0.4190678003222840d+0
       v=0.2271137107548774d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5215687136707969d+0
       b=0.4483151836883852d+0
       v=0.2283414092917525d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5491402346984905d+0
       b=0.4764740676087880d+0
       v=0.2291161673130077d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5753520160126075d+0
       b=0.5034021310998277d+0
       v=0.2295313908576598d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1388326356417754d+0
       b=0.2435436510372806d-1
       v=0.1438204721359031d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1743686900537244d+0
       b=0.5118897057342652d-1
       v=0.1607738025495257d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2099737037950268d+0
       b=0.8014695048539634d-1
       v=0.1741483853528379d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2454492590908548d+0
       b=0.1105117874155699d+0
       v=0.1851918467519151d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2807219257864278d+0
       b=0.1417950531570966d+0
       v=0.1944628638070613d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3156842271975842d+0
       b=0.1736604945719597d+0
       v=0.2022495446275152d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3502090945177752d+0
       b=0.2058466324693981d+0
       v=0.2087462382438514d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3841684849519686d+0
       b=0.2381284261195919d+0
       v=0.2141074754818308d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4174372367906016d+0
       b=0.2703031270422569d+0
       v=0.2184640913748162d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4498926465011892d+0
       b=0.3021845683091309d+0
       v=0.2219309165220329d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4814146229807701d+0
       b=0.3335993355165720d+0
       v=0.2246123118340624d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5118863625734701d+0
       b=0.3643833735518232d+0
       v=0.2266062766915125d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5411947455119144d+0
       b=0.3943789541958179d+0
       v=0.2280072952230796d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5692301500357246d+0
       b=0.4234320144403542d+0
       v=0.2289082025202583d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5958857204139576d+0
       b=0.4513897947419260d+0
       v=0.2294012695120025d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2156270284785766d+0
       b=0.2681225755444491d-1
       v=0.1722434488736947d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2532385054909710d+0
       b=0.5557495747805614d-1
       v=0.1830237421455091d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2902564617771537d+0
       b=0.8569368062950249d-1
       v=0.1923855349997633d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3266979823143256d+0
       b=0.1167367450324135d+0
       v=0.2004067861936271d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3625039627493614d+0
       b=0.1483861994003304d+0
       v=0.2071817297354263d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3975838937548699d+0
       b=0.1803821503011405d+0
       v=0.2128250834102103d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4318396099009774d+0
       b=0.2124962965666424d+0
       v=0.2174513719440102d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4651706555732742d+0
       b=0.2445221837805913d+0
       v=0.2211661839150214d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4974752649620969d+0
       b=0.2762701224322987d+0
       v=0.2240665257813102d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5286517579627517d+0
       b=0.3075627775211328d+0
       v=0.2262439516632620d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5586001195731895d+0
       b=0.3382311089826877d+0
       v=0.2277874557231869d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5872229902021319d+0
       b=0.3681108834741399d+0
       v=0.2287854314454994d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6144258616235123d+0
       b=0.3970397446872839d+0
       v=0.2293268499615575d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2951676508064861d+0
       b=0.2867499538750441d-1
       v=0.1912628201529828d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3335085485472725d+0
       b=0.5867879341903510d-1
       v=0.1992499672238701d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3709561760636381d+0
       b=0.8961099205022284d-1
       v=0.2061275533454027d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4074722861667498d+0
       b=0.1211627927626297d+0
       v=0.2119318215968572d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4429923648839117d+0
       b=0.1530748903554898d+0
       v=0.2167416581882652d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4774428052721736d+0
       b=0.1851176436721877d+0
       v=0.2206430730516600d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5107446539535904d+0
       b=0.2170829107658179d+0
       v=0.2237186938699523d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5428151370542935d+0
       b=0.2487786689026271d+0
       v=0.2260480075032884d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5735699292556964d+0
       b=0.2800239952795016d+0
       v=0.2277098884558542d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6029253794562866d+0
       b=0.3106445702878119d+0
       v=0.2287845715109671d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6307998987073145d+0
       b=0.3404689500841194d+0
       v=0.2293547268236294d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3752652273692719d+0
       b=0.2997145098184479d-1
       v=0.2056073839852528d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4135383879344028d+0
       b=0.6086725898678011d-1
       v=0.2114235865831876d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4506113885153907d+0
       b=0.9238849548435643d-1
       v=0.2163175629770551d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4864401554606072d+0
       b=0.1242786603851851d+0
       v=0.2203392158111650d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5209708076611709d+0
       b=0.1563086731483386d+0
       v=0.2235473176847839d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5541422135830122d+0
       b=0.1882696509388506d+0
       v=0.2260024141501235d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5858880915113817d+0
       b=0.2199672979126059d+0
       v=0.2277675929329182d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6161399390603444d+0
       b=0.2512165482924867d+0
       v=0.2289102112284834d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6448296482255090d+0
       b=0.2818368701871888d+0
       v=0.2295027954625118d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4544796274917948d+0
       b=0.3088970405060312d-1
       v=0.2161281589879992d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4919389072146628d+0
       b=0.6240947677636835d-1
       v=0.2201980477395102d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5279313026985183d+0
       b=0.9430706144280313d-1
       v=0.2234952066593166d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5624169925571135d+0
       b=0.1263547818770374d+0
       v=0.2260540098520838d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5953484627093287d+0
       b=0.1583430788822594d+0
       v=0.2279157981899988d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6266730715339185d+0
       b=0.1900748462555988d+0
       v=0.2291296918565571d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6563363204278871d+0
       b=0.2213599519592567d+0
       v=0.2297533752536649d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5314574716585696d+0
       b=0.3152508811515374d-1
       v=0.2234927356465995d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5674614932298185d+0
       b=0.6343865291465561d-1
       v=0.2261288012985219d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6017706004970264d+0
       b=0.9551503504223951d-1
       v=0.2280818160923688d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6343471270264178d+0
       b=0.1275440099801196d+0
       v=0.2293773295180159d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6651494599127802d+0
       b=0.1593252037671960d+0
       v=0.2300528767338634d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6050184986005704d+0
       b=0.3192538338496105d-1
       v=0.2281893855065666d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6390163550880400d+0
       b=0.6402824353962306d-1
       v=0.2295720444840727d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6711199107088448d+0
       b=0.9609805077002909d-1
       v=0.2303227649026753d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6741354429572275d+0
       b=0.3211853196273233d-1
       v=0.2304831913227114d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 5294-point angular grid
!
       subroutine ld5294(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(5294)
! y-coordinates of the gridpoints
       real(8) :: y(5294)
! z-coordinates of the gridpoints
       real(8) :: z(5294)
! weight of the grid point
       real(8) :: w(5294)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.9080510764308163d-4
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.2084824361987793d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2303261686261450d-1
       v=0.5011105657239616d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3757208620162394d-1
       v=0.5942520409683854d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5821912033821852d-1
       v=0.9564394826109721d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8403127529194872d-1
       v=0.1185530657126338d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1122927798060578d+0
       v=0.1364510114230331d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1420125319192987d+0
       v=0.1505828825605415d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1726396437341978d+0
       v=0.1619298749867023d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2038170058115696d+0
       v=0.1712450504267789d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2352849892876508d+0
       v=0.1789891098164999d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2668363354312461d+0
       v=0.1854474955629795d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2982941279900452d+0
       v=0.1908148636673661d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3295002922087076d+0
       v=0.1952377405281833d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3603094918363593d+0
       v=0.1988349254282232d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3905857895173920d+0
       v=0.2017079807160050d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4202005758160837d+0
       v=0.2039473082709094d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4490310061597227d+0
       v=0.2056360279288953d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4769586160311491d+0
       v=0.2068525823066865d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5038679887049750d+0
       v=0.2076724877534488d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5296454286519961d+0
       v=0.2081694278237885d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5541776207164850d+0
       v=0.2084157631219326d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5990467321921213d+0
       v=0.2084381531128593d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6191467096294587d+0
       v=0.2083476277129307d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6375251212901849d+0
       v=0.2082686194459732d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6540514381131168d+0
       v=0.2082475686112415d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6685899064391510d+0
       v=0.2083139860289915d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6810013009681648d+0
       v=0.2084745561831237d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6911469578730340d+0
       v=0.2087091313375890d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6988956915141736d+0
       v=0.2089718413297697d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7041335794868720d+0
       v=0.2092003303479793d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7067754398018567d+0
       v=0.2093336148263241d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3840368707853623d-1
       v=0.7591708117365267d-4
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9835485954117399d-1
       v=0.1083383968169186d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1665774947612998d+0
       v=0.1403019395292510d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2405702335362910d+0
       v=0.1615970179286436d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3165270770189046d+0
       v=0.1771144187504911d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3927386145645443d+0
       v=0.1887760022988168d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4678825918374656d+0
       v=0.1973474670768214d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5408022024266935d+0
       v=0.2033787661234659d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6104967445752438d+0
       v=0.2072343626517331d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6760910702685738d+0
       v=0.2091177834226918d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6655644120217392d-1
       b=0.1936508874588424d-1
       v=0.9316684484675566d-4
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9446246161270182d-1
       b=0.4252442002115869d-1
       v=0.1116193688682976d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1242651925452509d+0
       b=0.6806529315354374d-1
       v=0.1298623551559414d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1553438064846751d+0
       b=0.9560957491205369d-1
       v=0.1450236832456426d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1871137110542670d+0
       b=0.1245931657452888d+0
       v=0.1572719958149914d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2192612628836257d+0
       b=0.1545385828778978d+0
       v=0.1673234785867195d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2515682807206955d+0
       b=0.1851004249723368d+0
       v=0.1756860118725188d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2838535866287290d+0
       b=0.2160182608272384d+0
       v=0.1826776290439367d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3159578817528521d+0
       b=0.2470799012277111d+0
       v=0.1885116347992865d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3477370882791392d+0
       b=0.2781014208986402d+0
       v=0.1933457860170574d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3790576960890540d+0
       b=0.3089172523515731d+0
       v=0.1973060671902064d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4097938317810200d+0
       b=0.3393750055472244d+0
       v=0.2004987099616311d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4398256572859637d+0
       b=0.3693322470987730d+0
       v=0.2030170909281499d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4690384114718480d+0
       b=0.3986541005609877d+0
       v=0.2049461460119080d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4973216048301053d+0
       b=0.4272112491408562d+0
       v=0.2063653565200186d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5245681526132446d+0
       b=0.4548781735309936d+0
       v=0.2073507927381027d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5506733911803888d+0
       b=0.4815315355023251d+0
       v=0.2079764593256122d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5755339829522475d+0
       b=0.5070486445801855d+0
       v=0.2083150534968778d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1305472386056362d+0
       b=0.2284970375722366d-1
       v=0.1262715121590664d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1637327908216477d+0
       b=0.4812254338288384d-1
       v=0.1414386128545972d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1972734634149637d+0
       b=0.7531734457511935d-1
       v=0.1538740401313898d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2308694653110130d+0
       b=0.1039043639882017d+0
       v=0.1642434942331432d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2643899218338160d+0
       b=0.1334526587117626d+0
       v=0.1729790609237496d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2977171599622171d+0
       b=0.1636414868936382d+0
       v=0.1803505190260828d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3307293903032310d+0
       b=0.1942195406166568d+0
       v=0.1865475350079657d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3633069198219073d+0
       b=0.2249752879943753d+0
       v=0.1917182669679069d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3953346955922727d+0
       b=0.2557218821820032d+0
       v=0.1959851709034382d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4267018394184914d+0
       b=0.2862897925213193d+0
       v=0.1994529548117882d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4573009622571704d+0
       b=0.3165224536636518d+0
       v=0.2022138911146548d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4870279559856109d+0
       b=0.3462730221636496d+0
       v=0.2043518024208592d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5157819581450322d+0
       b=0.3754016870282835d+0
       v=0.2059450313018110d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5434651666465393d+0
       b=0.4037733784993613d+0
       v=0.2070685715318472d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5699823887764627d+0
       b=0.4312557784139123d+0
       v=0.2077955310694373d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5952403350947741d+0
       b=0.4577175367122110d+0
       v=0.2081980387824712d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2025152599210369d+0
       b=0.2520253617719557d-1
       v=0.1521318610377956d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2381066653274425d+0
       b=0.5223254506119000d-1
       v=0.1622772720185755d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2732823383651612d+0
       b=0.8060669688588620d-1
       v=0.1710498139420709d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3080137692611118d+0
       b=0.1099335754081255d+0
       v=0.1785911149448736d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3422405614587601d+0
       b=0.1399120955959857d+0
       v=0.1850125313687736d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3758808773890420d+0
       b=0.1702977801651705d+0
       v=0.1904229703933298d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4088458383438932d+0
       b=0.2008799256601680d+0
       v=0.1949259956121987d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4410450550841152d+0
       b=0.2314703052180836d+0
       v=0.1986161545363960d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4723879420561312d+0
       b=0.2618972111375892d+0
       v=0.2015790585641370d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5027843561874343d+0
       b=0.2920013195600270d+0
       v=0.2038934198707418d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5321453674452458d+0
       b=0.3216322555190551d+0
       v=0.2056334060538251d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5603839113834030d+0
       b=0.3506456615934198d+0
       v=0.2068705959462289d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5874150706875146d+0
       b=0.3789007181306267d+0
       v=0.2076753906106002d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6131559381660038d+0
       b=0.4062580170572782d+0
       v=0.2081179391734803d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2778497016394506d+0
       b=0.2696271276876226d-1
       v=0.1700345216228943d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3143733562261912d+0
       b=0.5523469316960465d-1
       v=0.1774906779990410d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3501485810261827d+0
       b=0.8445193201626464d-1
       v=0.1839659377002642d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3851430322303653d+0
       b=0.1143263119336083d+0
       v=0.1894987462975169d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4193013979470415d+0
       b=0.1446177898344475d+0
       v=0.1941548809452595d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4525585960458567d+0
       b=0.1751165438438091d+0
       v=0.1980078427252384d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4848447779622947d+0
       b=0.2056338306745660d+0
       v=0.2011296284744488d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5160871208276894d+0
       b=0.2359965487229226d+0
       v=0.2035888456966776d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5462112185696926d+0
       b=0.2660430223139146d+0
       v=0.2054516325352142d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5751425068101757d+0
       b=0.2956193664498032d+0
       v=0.2067831033092635d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6028073872853596d+0
       b=0.3245763905312779d+0
       v=0.2076485320284876d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6291338275278409d+0
       b=0.3527670026206972d+0
       v=0.2081141439525255d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3541797528439391d+0
       b=0.2823853479435550d-1
       v=0.1834383015469222d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3908234972074657d+0
       b=0.5741296374713106d-1
       v=0.1889540591777677d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4264408450107590d+0
       b=0.8724646633650199d-1
       v=0.1936677023597375d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4609949666553286d+0
       b=0.1175034422915616d+0
       v=0.1976176495066504d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4944389496536006d+0
       b=0.1479755652628428d+0
       v=0.2008536004560983d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5267194884346086d+0
       b=0.1784740659484352d+0
       v=0.2034280351712291d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5577787810220990d+0
       b=0.2088245700431244d+0
       v=0.2053944466027758d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5875563763536670d+0
       b=0.2388628136570763d+0
       v=0.2068077642882360d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6159910016391269d+0
       b=0.2684308928769185d+0
       v=0.2077250949661599d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6430219602956268d+0
       b=0.2973740761960252d+0
       v=0.2082062440705320d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4300647036213646d+0
       b=0.2916399920493977d-1
       v=0.1934374486546626d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4661486308935531d+0
       b=0.5898803024755659d-1
       v=0.1974107010484300d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5009658555287261d+0
       b=0.8924162698525409d-1
       v=0.2007129290388658d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5344824270447704d+0
       b=0.1197185199637321d+0
       v=0.2033736947471293d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5666575997416371d+0
       b=0.1502300756161382d+0
       v=0.2054287125902493d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5974457471404752d+0
       b=0.1806004191913564d+0
       v=0.2069184936818894d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6267984444116886d+0
       b=0.2106621764786252d+0
       v=0.2078883689808782d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6546664713575417d+0
       b=0.2402526932671914d+0
       v=0.2083886366116359d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5042711004437253d+0
       b=0.2982529203607657d-1
       v=0.2006593275470817d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5392127456774380d+0
       b=0.6008728062339922d-1
       v=0.2033728426135397d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5726819437668618d+0
       b=0.9058227674571398d-1
       v=0.2055008781377608d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6046469254207278d+0
       b=0.1211219235803400d+0
       v=0.2070651783518502d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6350716157434952d+0
       b=0.1515286404791580d+0
       v=0.2080953335094320d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6639177679185454d+0
       b=0.1816314681255552d+0
       v=0.2086284998988521d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5757276040972253d+0
       b=0.3026991752575440d-1
       v=0.2055549387644668d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6090265823139755d+0
       b=0.6078402297870770d-1
       v=0.2071871850267654d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6406735344387661d+0
       b=0.9135459984176636d-1
       v=0.2082856600431965d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6706397927793709d+0
       b=0.1218024155966590d+0
       v=0.2088705858819358d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6435019674426665d+0
       b=0.3052608357660639d-1
       v=0.2083995867536322d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6747218676375681d+0
       b=0.6112185773983089d-1
       v=0.2090509712889637d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end
!
!    Lebedev 5810-point angular grid
!
       subroutine ld5810(x,y,z,w,n)
       implicit none
! x-coordinates of the gridpoints
       real(8) :: x(5810)
! y-coordinates of the gridpoints
       real(8) :: y(5810)
! z-coordinates of the gridpoints
       real(8) :: z(5810)
! weight of the grid point
       real(8) :: w(5810)
! number of grid points
       integer(4) :: n
! internal variables
       real(8) :: a,b,v
       n=1
       v=0.9735347946175486d-5
       call gen_oh( 1, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1907581241803167d-3
       call gen_oh( 2, n, x(n), y(n), z(n), w(n), a, b, v)
       v=0.1901059546737578d-3
       call gen_oh( 3, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1182361662400277d-1
       v=0.3926424538919212d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3062145009138958d-1
       v=0.6667905467294382d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5329794036834243d-1
       v=0.8868891315019135d-4
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7848165532862220d-1
       v=0.1066306000958872d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1054038157636201d+0
       v=0.1214506743336128d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1335577797766211d+0
       v=0.1338054681640871d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1625769955502252d+0
       v=0.1441677023628504d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1921787193412792d+0
       v=0.1528880200826557d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2221340534690548d+0
       v=0.1602330623773609d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2522504912791132d+0
       v=0.1664102653445244d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2823610860679697d+0
       v=0.1715845854011323d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3123173966267560d+0
       v=0.1758901000133069d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3419847036953789d+0
       v=0.1794382485256736d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3712386456999758d+0
       v=0.1823238106757407d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3999627649876828d+0
       v=0.1846293252959976d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4280466458648093d+0
       v=0.1864284079323098d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4553844360185711d+0
       v=0.1877882694626914d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4818736094437834d+0
       v=0.1887716321852025d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5074138709260629d+0
       v=0.1894381638175673d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5319061304570707d+0
       v=0.1898454899533629d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5552514978677286d+0
       v=0.1900497929577815d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5981009025246183d+0
       v=0.1900671501924092d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6173990192228116d+0
       v=0.1899837555533510d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6351365239411131d+0
       v=0.1899014113156229d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6512010228227200d+0
       v=0.1898581257705106d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6654758363948120d+0
       v=0.1898804756095753d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6778410414853370d+0
       v=0.1899793610426402d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6881760887484110d+0
       v=0.1901464554844117d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6963645267094598d+0
       v=0.1903533246259542d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7023010617153579d+0
       v=0.1905556158463228d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.7059004636628753d+0
       v=0.1907037155663528d-3
       call gen_oh( 4, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3552470312472575d-1
       v=0.5992997844249967d-4
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.9151176620841283d-1
       v=0.9749059382456978d-4
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1566197930068980d+0
       v=0.1241680804599158d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2265467599271907d+0
       v=0.1437626154299360d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2988242318581361d+0
       v=0.1584200054793902d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3717482419703886d+0
       v=0.1694436550982744d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4440094491758889d+0
       v=0.1776617014018108d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5145337096756642d+0
       v=0.1836132434440077d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5824053672860230d+0
       v=0.1876494727075983d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6468283961043370d+0
       v=0.1899906535336482d-3
       call gen_oh( 5, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6095964259104373d-1
       b=0.1787828275342931d-1
       v=0.8143252820767350d-4
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.8811962270959388d-1
       b=0.3953888740792096d-1
       v=0.9998859890887728d-4
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1165936722428831d+0
       b=0.6378121797722990d-1
       v=0.1156199403068359d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1460232857031785d+0
       b=0.8985890813745037d-1
       v=0.1287632092635513d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1761197110181755d+0
       b=0.1172606510576162d+0
       v=0.1398378643365139d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2066471190463718d+0
       b=0.1456102876970995d+0
       v=0.1491876468417391d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2374076026328152d+0
       b=0.1746153823011775d+0
       v=0.1570855679175456d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2682305474337051d+0
       b=0.2040383070295584d+0
       v=0.1637483948103775d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2989653312142369d+0
       b=0.2336788634003698d+0
       v=0.1693500566632843d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3294762752772209d+0
       b=0.2633632752654219d+0
       v=0.1740322769393633d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3596390887276086d+0
       b=0.2929369098051601d+0
       v=0.1779126637278296d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3893383046398812d+0
       b=0.3222592785275512d+0
       v=0.1810908108835412d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4184653789358347d+0
       b=0.3512004791195743d+0
       v=0.1836529132600190d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4469172319076166d+0
       b=0.3796385677684537d+0
       v=0.1856752841777379d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4745950813276976d+0
       b=0.4074575378263879d+0
       v=0.1872270566606832d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5014034601410262d+0
       b=0.4345456906027828d+0
       v=0.1883722645591307d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5272493404551239d+0
       b=0.4607942515205134d+0
       v=0.1891714324525297d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5520413051846366d+0
       b=0.4860961284181720d+0
       v=0.1896827480450146d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5756887237503077d+0
       b=0.5103447395342790d+0
       v=0.1899628417059528d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1225039430588352d+0
       b=0.2136455922655793d-1
       v=0.1123301829001669d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1539113217321372d+0
       b=0.4520926166137188d-1
       v=0.1253698826711277d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1856213098637712d+0
       b=0.7086468177864818d-1
       v=0.1366266117678531d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2174998728035131d+0
       b=0.9785239488772918d-1
       v=0.1462736856106918d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2494128336938330d+0
       b=0.1258106396267210d+0
       v=0.1545076466685412d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2812321562143480d+0
       b=0.1544529125047001d+0
       v=0.1615096280814007d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3128372276456111d+0
       b=0.1835433512202753d+0
       v=0.1674366639741759d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3441145160177973d+0
       b=0.2128813258619585d+0
       v=0.1724225002437900d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3749567714853510d+0
       b=0.2422913734880829d+0
       v=0.1765810822987288d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4052621732015610d+0
       b=0.2716163748391453d+0
       v=0.1800104126010751d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4349335453522385d+0
       b=0.3007127671240280d+0
       v=0.1827960437331284d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4638776641524965d+0
       b=0.3294470677216479d+0
       v=0.1850140300716308d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4920046410462687d+0
       b=0.3576932543699155d+0
       v=0.1867333507394938d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5192273554861704d+0
       b=0.3853307059757764d+0
       v=0.1880178688638289d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5454609081136522d+0
       b=0.4122425044452694d+0
       v=0.1889278925654758d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5706220661424140d+0
       b=0.4383139587781027d+0
       v=0.1895213832507346d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5946286755181518d+0
       b=0.4634312536300553d+0
       v=0.1898548277397420d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.1905370790924295d+0
       b=0.2371311537781979d-1
       v=0.1349105935937341d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2242518717748009d+0
       b=0.4917878059254806d-1
       v=0.1444060068369326d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2577190808025936d+0
       b=0.7595498960495142d-1
       v=0.1526797390930008d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2908724534927187d+0
       b=0.1036991083191100d+0
       v=0.1598208771406474d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3236354020056219d+0
       b=0.1321348584450234d+0
       v=0.1659354368615331d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3559267359304543d+0
       b=0.1610316571314789d+0
       v=0.1711279910946440d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3876637123676956d+0
       b=0.1901912080395707d+0
       v=0.1754952725601440d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4187636705218842d+0
       b=0.2194384950137950d+0
       v=0.1791247850802529d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4491449019883107d+0
       b=0.2486155334763858d+0
       v=0.1820954300877716d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4787270932425445d+0
       b=0.2775768931812335d+0
       v=0.1844788524548449d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5074315153055574d+0
       b=0.3061863786591120d+0
       v=0.1863409481706220d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5351810507738336d+0
       b=0.3343144718152556d+0
       v=0.1877433008795068d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5619001025975381d+0
       b=0.3618362729028427d+0
       v=0.1887444543705232d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5875144035268046d+0
       b=0.3886297583620408d+0
       v=0.1894009829375006d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6119507308734495d+0
       b=0.4145742277792031d+0
       v=0.1897683345035198d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2619733870119463d+0
       b=0.2540047186389353d-1
       v=0.1517327037467653d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.2968149743237949d+0
       b=0.5208107018543989d-1
       v=0.1587740557483543d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3310451504860488d+0
       b=0.7971828470885599d-1
       v=0.1649093382274097d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3646215567376676d+0
       b=0.1080465999177927d+0
       v=0.1701915216193265d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3974916785279360d+0
       b=0.1368413849366629d+0
       v=0.1746847753144065d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4295967403772029d+0
       b=0.1659073184763559d+0
       v=0.1784555512007570d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4608742854473447d+0
       b=0.1950703730454614d+0
       v=0.1815687562112174d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4912598858949903d+0
       b=0.2241721144376724d+0
       v=0.1840864370663302d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5206882758945558d+0
       b=0.2530655255406489d+0
       v=0.1860676785390006d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5490940914019819d+0
       b=0.2816118409731066d+0
       v=0.1875690583743703d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5764123302025542d+0
       b=0.3096780504593238d+0
       v=0.1886453236347225d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6025786004213506d+0
       b=0.3371348366394987d+0
       v=0.1893501123329645d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6275291964794956d+0
       b=0.3638547827694396d+0
       v=0.1897366184519868d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3348189479861771d+0
       b=0.2664841935537443d-1
       v=0.1643908815152736d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.3699515545855295d+0
       b=0.5424000066843495d-1
       v=0.1696300350907768d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4042003071474669d+0
       b=0.8251992715430854d-1
       v=0.1741553103844483d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4375320100182624d+0
       b=0.1112695182483710d+0
       v=0.1780015282386092d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4699054490335947d+0
       b=0.1402964116467816d+0
       v=0.1812116787077125d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5012739879431952d+0
       b=0.1694275117584291d+0
       v=0.1838323158085421d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5315874883754966d+0
       b=0.1985038235312689d+0
       v=0.1859113119837737d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5607937109622117d+0
       b=0.2273765660020893d+0
       v=0.1874969220221698d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5888393223495521d+0
       b=0.2559041492849764d+0
       v=0.1886375612681076d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6156705979160163d+0
       b=0.2839497251976899d+0
       v=0.1893819575809276d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6412338809078123d+0
       b=0.3113791060500690d+0
       v=0.1897794748256767d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4076051259257167d+0
       b=0.2757792290858463d-1
       v=0.1738963926584846d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4423788125791520d+0
       b=0.5584136834984293d-1
       v=0.1777442359873466d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4760480917328258d+0
       b=0.8457772087727143d-1
       v=0.1810010815068719d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5085838725946297d+0
       b=0.1135975846359248d+0
       v=0.1836920318248129d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5399513637391218d+0
       b=0.1427286904765053d+0
       v=0.1858489473214328d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5701118433636380d+0
       b=0.1718112740057635d+0
       v=0.1875079342496592d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5990240530606021d+0
       b=0.2006944855985351d+0
       v=0.1887080239102310d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6266452685139695d+0
       b=0.2292335090598907d+0
       v=0.1894905752176822d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6529320971415942d+0
       b=0.2572871512353714d+0
       v=0.1898991061200695d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.4791583834610126d+0
       b=0.2826094197735932d-1
       v=0.1809065016458791d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5130373952796940d+0
       b=0.5699871359683649d-1
       v=0.1836297121596799d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5456252429628476d+0
       b=0.8602712528554394d-1
       v=0.1858426916241869d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5768956329682385d+0
       b=0.1151748137221281d+0
       v=0.1875654101134641d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6068186944699046d+0
       b=0.1442811654136362d+0
       v=0.1888240751833503d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6353622248024907d+0
       b=0.1731930321657680d+0
       v=0.1896497383866979d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6624927035731797d+0
       b=0.2017619958756061d+0
       v=0.1900775530219121d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5484933508028488d+0
       b=0.2874219755907391d-1
       v=0.1858525041478814d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.5810207682142106d+0
       b=0.5778312123713695d-1
       v=0.1876248690077947d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6120955197181352d+0
       b=0.8695262371439526d-1
       v=0.1889404439064607d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6416944284294319d+0
       b=0.1160893767057166d+0
       v=0.1898168539265290d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6697926391731260d+0
       b=0.1450378826743251d+0
       v=0.1902779940661772d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6147594390585488d+0
       b=0.2904957622341456d-1
       v=0.1890125641731815d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6455390026356783d+0
       b=0.5823809152617197d-1
       v=0.1899434637795751d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6747258588365477d+0
       b=0.8740384899884715d-1
       v=0.1904520856831751d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       a=0.6772135750395347d+0
       b=0.2919946135808105d-1
       v=0.1905534498734563d-3
       call gen_oh( 6, n, x(n), y(n), z(n), w(n), a, b, v)
       n=n-1
       return
       end

