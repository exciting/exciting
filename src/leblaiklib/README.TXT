Lebedev grids of orders n=6m+5 where m=0,1,...,21 in 16 digit precision
=======================================================================

The file Lebedev-Laikov.F implements a set of subroutines providing 
Lebedev-Laikov grids of order n=2m+1, where m=1,2,...,15, and additionally
grids of order n=6m+5, where m=5,6,...,21. The parameters ensure 
that angular integration of polynomials x**k * y**l * z**m, where k+l+m <= 131 
can be performed with a relative accuracy of 2e-14 [1]. Note that the weights
are normalised to add up to 1.0.

For each order n a separate subroutine is provided named 
LD. The parameters X, Y, Z are arrays for the 
cartesian components of each point, and the parameter W is an array for the
weights. The subroutines increase the integer parameter N by number of grid
points generated. All these routines use the subroutine gen_oh which takes care 
of the octahedral symmetry of the grids.

Christoph van Wuellen (Ruhr-Universitaet, Bochum, Germany) generated the 
routines in Lebedev-Laikov.F by translating the original C-routines kindly 
provided by Dmitri Laikov (Moscow State University, Moscow, Russia). We 
are in debt to Dmitri Laikov for giving us permission to make these routines
publically available.

Huub van Dam
Daresbury Laboratory, Daresbury, United Kingdom
April, 2000

References
==========

[1] V.I. Lebedev, and D.N. Laikov
    "A quadrature formula for the sphere of the 131st
     algebraic order of accuracy"
    Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
