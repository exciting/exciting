      FUNCTION DLAMCH ( CMACH ) RESULT(RMACH)
!
!  -- LAPACK auxiliary routine Replacement for DLAMCH.f
!     use Fortran 90 Machine Parameter built-in functions.
!
      character(len=1) :: CMACH
      real(kind(1.d0)) :: xdbl, ydbl, RMACH
!
!  Purpose
!  =======
!
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!
! =====================================================================
      xdbl=1.d0
      IF( CMACH == 'E' .or. CMACH == 'e' ) THEN
!          eps   = relative machine precision
         RMACH = Epsilon(xdbl)
      ELSE IF( CMACH == 'S' .or. CMACH == 's' ) THEN
!          sfmin = safe minimum, such that 1/sfmin does not overflow
         RMACH = Tiny(xdbl)
      ELSE IF( CMACH == 'B' .or. CMACH == 'b' ) THEN
!          base  = base of the machine
         RMACH = Radix(xdbl)
      ELSE IF( CMACH == 'P' .or. CMACH == 'p' ) THEN
!          prec  = eps*base
         RMACH = Radix(xdbl)*Epsilon(xdbl)
      ELSE IF( CMACH == 'N' .or. CMACH == 'n' ) THEN
!          t     = number of (base) digits in the mantissa
         RMACH = Digits(xdbl)
      ELSE IF( CMACH == 'R' .or. CMACH == 'r' ) THEN
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          Assume rounding (IEEE).
         RMACH = 1.0
      ELSE IF( CMACH == 'M' .or. CMACH == 'm' ) THEN
!          emin  = minimum exponent before (gradual) underflow
         RMACH = Minexponent(xdbl)
      ELSE IF( CMACH == 'U' .or. CMACH == 'u' ) THEN
!          rmin  = underflow threshold - base**(emin-1)
         RMACH = Tiny(xdbl)
      ELSE IF( CMACH == 'L' .or. CMACH == 'l' ) THEN
!          emax  = largest exponent before overflow
         RMACH = Maxexponent(xdbl)
      ELSE IF( CMACH == 'O' .or. CMACH == 'o' ) THEN
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
         RMACH = Huge(xdbl)
      END IF
!
      END

*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*
*     .. Local variables ..
      INTEGER I
      DOUBLE PRECISION X(10),Y(10)
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
*      DLAMC3 = A + B

*     Modification by JKD to ensure variables are flushed to memory
      DO I=1,10
        X(I)=A
      END DO
      DO I=1,10
        Y(I)=X(I)+B
      END DO
      DLAMC3=Y(10)
*
      RETURN
*
*     End of DLAMC3
*
      END

