!
!
!
!
!
      Subroutine gencoulcut (cctype, rccut, ngpmax, ngp, vgpc, gpc, ccut_)
      Use modmain
      Use modmpi
      Use modxs
      Implicit None
! arguments
      Character (*), Intent (In) :: cctype
      Integer, Intent (In) :: ngpmax, ngp
      Real (8), Intent (In) :: vgpc (3, ngpmax), gpc (ngpmax)
      Real (8), Intent (Out) :: ccut_ (ngpmax)
 ! local variables
      Integer :: igp
      Real (8) :: rccut, b0, maxgpc, gmax
!      Real (8), Allocatable :: samplej0(:,:)
      Integer :: n,n0
! external functions
      Real (8), External :: ccut



!     sample bessel_J0 function for numerical integrations in the 1D case.
!     If(cctype .eq. '1d') then
!     n0 = 20
!     !location of first zero of bessel_j0(x)
!     b0 = 2.4048
!     n = nint( n0 * max(1, rccut * maxval())
!     Allocate(samplej0(n,2))
!     End If
      gmax = maxval(gpc)
        write(*,*) "cctype:",cctype
        write(*,*) "len cctype:",len(cctype)

      Do igp = 1, ngp
        write(*,*) "calculating cutoff for (rank, ig) = ",rank, igp
         ccut_ (igp) = ccut (cctype, rccut, vgpc(:, igp), gmax, gpc(igp))
      End Do
      End Subroutine gencoulcut
!
!
      Real (8) Function ccut (cctype, rccut, vgpc, gmax, gpc)
      Use modmain, Only: pi, twopi, fourpi
      Use modmpi
      Use modinput
      Implicit None
! arguments
      Character (*), Intent (In) :: cctype
      Real (8), Intent (In) :: vgpc (3), gpc, rccut, gmax
! local variables
      real(8), parameter :: eps=1.d-8
      Real (8) :: t1, gx, agx, gper, gpar, gz, gr
      Real (8) :: x, dx, b0, h
      Integer ::  n0, n, i
! external functions
      Real (8), External :: besk0, besk1
      
      !write(*,*) "cctype in function:", cctype
      Select Case (trim(cctype))
      Case ('nocutoff')
     !        set up the square root of the Coulomb potential from analytical
     !        expression (no cutoff)
         ccut = 0.d0
         if (gpc .gt. eps) then
            ccut = Sqrt (fourpi) / gpc
         else
! set sqrt of Coulomb potential to zero for |G+q| = 0
            ccut = 0.d0
         end if
      Case ('0d')
     !        0D spherical cutoff

         if (gpc .gt. eps) then
            ccut = fourpi / gpc**2
            ccut = ccut * (1.d0 - Cos(gpc * rccut))
         else
            ccut = twopi * rccut**2
         end if

      Case ('1d')
     !        1D infinite cylinder
         Write (*,*)
         Write (*, '("Warning(genccut): 1D cutoff applied. Assumes &
         & periodicity along x direction")')
         Write (*,*)
         gx = abs(vgpc(1))
         gper = Sqrt( vgpc(2)**2 + vgpc(3)**2 )

         if (gx .gt. eps) then
            ccut = fourpi / gpc**2
            t1 = 1.d0
            t1 = t1 + gper*rccut * bessel_jn(1,gper*rccut) * besk0(gx*rccut)
            t1 = t1 - gx*rccut * bessel_jn(0,gper*rccut) * besk1(gx*rccut)
            ccut = ccut * t1
         else if ( (gx .le.eps) .and. (gper .gt. eps) ) then
!     location of first zero of bessel_J0(x)
            b0 = 2.4048
!     desired number of intervals in each bessel_j0 cycle
            n0 = 60
!     maximum gper

            !gmax = gqmax

!     total number of integration intervals.
            n = nint( n0 * max(1.0, rccut * gmax/b0))
!     length of the finite cylinder (assumes periodicity of the wire
!     along x direction, defined by the first basevect of form (lx,0,0))
            h = 2 * input%groundstate%ngridk(1) * &
            & input%structure%crystal%basevect(1,1)

            write(*,*) "max g = ", gmax
            write(*,*) "gqmax = ", input%xs%gqmax

            dx = min(rccut/(n0*1.d0) , b0/(n0 * gmax))

            ccut = 0.d0
            do i = 1, n
               x = i*dx
               ccut = ccut + x * bessel_jn(0,x*gper) * log(x)
            end do
            ccut = -fourpi * ccut * dx

            ccut = ccut + fourpi * rccut * log(2*h) * &
            & bessel_jn(1,gper * rccut)/gper

         else if( (gx .le.eps) .and. (gper .le. eps) ) then
            ccut = - pi * rccut**2 * (2 * log(rccut) - 1)
         end if

            write (*,*) "gx = ",gx
            write (*,*) "gper = ",gper
            write (*,*) "ccut = ",ccut
      Case ('2d')
     !        2D infinite slab
         !Write (*,*)
         !Write (*, '("Warning(genccut): 2D cutoff applied. Assumes &
         !& periodicity along xy direction")')
         !Write (*,*)

         gpar = Sqrt( vgpc(1)**2 + vgpc(2)**2 )
         gper = abs(vgpc(3))
         gr = gper * rccut

         if (gpar .gt. eps) then
            t1 = gper/gpar * sin(gr) - cos(gr)
            ccut = fourpi / gpc**2
            ccut = ccut * ( 1.d0 + exp (-gpar * rccut) * t1 )
         else if ((gpar .le. eps) .and. (gper .gt. eps)) then
            t1 = -cos(gr) - gr * sin(gr)
            ccut = fourpi / gper**2
            ccut = ccut * ( 1.d0 + t1 )
         else
            ccut = -fourpi * rccut **2 / 2.d0
         end if

      Case Default
         Write (*,*)
         Write (*, '("Error(gencoulcut): unknown cutoff type for Coulomb &
         &potential: ", a)') cctype
         Write (*,*)
         Call terminate
      End Select
!
      End Function ccut

!
!
!
! This code is based on a FORTRAN90 library SPECFUN which evaluates certain special
! functions
! Original FORTRAN77 version by William Cody, Laura Stoltz.
! FORTRAN90 version by John Burkardt.
!
! The original (FORTRAN77) version of SPECFUN is available
! through NETLIB: http://www.netlib.org/specfun/index.html".
!
! This piece of code is available under the GNU LGPL license.
!
!BOP
! !ROUTINE: besk0
! !INTERFACE:
Real (8) Function besk0 (x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   besk0 evaluates the Bessel K0(X) function.
!   Returns the error function ${\rm erf}(x)$ using a rational function
!   approximation. This procedure is numerically stable and accurate to near
!   machine precision.
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order zero
!    for arguments 0.0 < ARG <= XMAX.
!    See comments heading CALCK0.
!
! !REVISION HISTORY:
!EOP
!BOC

  implicit none

  integer (4) jint
  Real (8) result
  Real (8) x

  jint = 1
  call calck0 ( x, result, jint )
  besk0 = result

  return
End Function
!EOC

!
!
!
! This code is based on a FORTRAN90 library SPECFUN which evaluates certain special
! functions
! Original FORTRAN77 version by William Cody, Laura Stoltz.
! FORTRAN90 version by John Burkardt.
!
! The original (FORTRAN77) version of SPECFUN is available
! through NETLIB: http://www.netlib.org/specfun/index.html".
!
! This piece of code is available under the GNU LGPL license.
!
!BOP
! !ROUTINE: besk1
! !INTERFACE:
Real (8) Function besk1 ( x )
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   besk1 evaluates the Bessel K1(X) function.
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order one
!    for arguments XLEAST <= ARG <= XMAX.
!
! !REVISION HISTORY:
!EOP
!BOC
  implicit none

  integer (4) jint
  real (8) result
  real (8) x

  jint = 1
  call calck1 ( x, result, jint )
  besk1 = result

  return
End Function
!EOC

!
!
!
! This code is based on a FORTRAN90 library SPECFUN which evaluates certain special
! functions
! Original FORTRAN77 version by William Cody, Laura Stoltz.
! FORTRAN90 version by John Burkardt.
!
! The original (FORTRAN77) version of SPECFUN is available
! through NETLIB: http://www.netlib.org/specfun/index.html".
!
! This piece of code is available under the GNU LGPL license.
!
!BOP
! !ROUTINE: calck0
! !INTERFACE:
subroutine calck0 ( arg, result, jint )
! !INPUT/OUTPUT PARAMETERS:
!   arg : real argument (in,real).
!    0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!   result : real argument (in,real)
!    Input, real ( kind = 8 ) ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K0(x);
!    2, RESULT = exp(x) * K0(x);
!
!   jint : real argument (in,real)
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
!
!
! !DESCRIPTION:
!
!    CALCK0 computes various K0 Bessel functions.
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
! !REVISION HISTORY:
!EOP
!BOC

  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) arg
  real ( kind = 8 ) f(4)
  real ( kind = 8 ) g(3)
  real ( kind = 8 ) p(6)
  real ( kind = 8 ) pp(10)
  real ( kind = 8 ) q(2)
  real ( kind = 8 ) qq(10)
  real ( kind = 8 ) result
  real ( kind = 8 ) sumf
  real ( kind = 8 ) sumg
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xx
!
!  Machine-dependent constants
!
  data xsmall /1.11d-16/
  data xinf /1.79d+308/
  data xmax /705.342d0/
!
!  Coefficients for XSMALL <= ARG <= 1.0
!
  data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01, &
            1.1999463724910714109d+01, 4.6850901201934832188d+02, &
            5.9169059852270512312d+03, 2.4708152720399552679d+03/
  data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
  data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02, &
           -1.7733784684952985886d+04,-4.0320340761145482298d+05/
  data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
           -1.6128136304458193998d+06/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03, &
            3.1075408980684392399d+04, 1.0577068948034021957d+05, &
            1.7398867902565686251d+05, 1.5097646353289914539d+05, &
            7.1557062783764037541d+04, 1.8321525870183537725d+04, &
            2.3444738764199315021d+03, 1.1600249425076035558d+02/
  data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03, &
            3.1474655750295278825d+04, 9.7418829762268075784d+04, &
            1.5144644673520157801d+05, 1.2689839587977598727d+05, &
            5.8824616785857027752d+04, 1.4847228371802360957d+04, &
            1.8821890840982713696d+03, 9.2556599177304839811d+01/

  x = arg
!
!  0.0 < ARG <= 1.0.
!
  if ( 0.0D+00 < x ) then

    if ( x <= 1.0D+00 ) then

      temp = log ( x )

      if ( x < xsmall ) then
!
!  Return for small ARG.
!
        result = p(6) / q(2) - temp

      else

        xx = x * x

        sump = (((( &
                 p(1) &
          * xx + p(2) ) &
          * xx + p(3) ) &
          * xx + p(4) ) &
          * xx + p(5) ) &
          * xx + p(6)

        sumq = ( xx + q(1) ) * xx + q(2)
        sumf = ( ( &
                 f(1) &
          * xx + f(2) ) &
          * xx + f(3) ) &
          * xx + f(4)

        sumg = ( ( xx + g(1) ) * xx + g(2) ) * xx + g(3)

        result = sump / sumq - xx * sumf * temp / sumg - temp

        if ( jint == 2 ) then
          result = result * exp ( x )
        end if

      end if

    else if ( jint == 1 .and. xmax < x ) then
!
!  Error return for XMAX < ARG.
!
      result = 0.0D+00

    else
!
!  1.0 < ARG.
!
      xx = 1.0D+00 / x
      sump = pp(1)
      do i = 2, 10
        sump = sump * xx + pp(i)
      end do

      sumq = xx
      do i = 1, 9
        sumq = ( sumq + qq(i) ) * xx
      end do
      sumq = sumq + qq(10)
      result = sump / sumq / sqrt ( x )

      if ( jint == 1 ) then
        result = result * exp ( -x )
      end if

    end if

  else
!
!  Error return for ARG <= 0.0.
!
    result = xinf

  end if

  return
End Subroutine
!EOC

!
!
!
! This code is based on a FORTRAN90 library SPECFUN which evaluates certain special
! functions
! Original FORTRAN77 version by William Cody, Laura Stoltz.
! FORTRAN90 version by John Burkardt.
!
! The original (FORTRAN77) version of SPECFUN is available
! through NETLIB: http://www.netlib.org/specfun/index.html".
!
! This piece of code is available under the GNU LGPL license.
!
!BOP
! !ROUTINE: calck0
! !INTERFACE:
subroutine calck1 ( arg, result, jint )
! !INPUT/OUTPUT PARAMETERS:
!   arg : real argument (in,real).
!    0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!   result : real argument (in,real)
!    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K1(x);
!    2, RESULT = exp(x) * K1(x);
!
!   jint : real argument (in,real)
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);
!
! !DESCRIPTION:
!
!    CALCK1 computes various K1 Bessel functions.
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
! !REVISION HISTORY:
!EOP
!BOC

  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) f(5)
  real ( kind = 8 ) g(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) p(5)
  real ( kind = 8 ) pp(11)
  real ( kind = 8 ) q(3)
  real ( kind = 8 ) qq(9)
  real ( kind = 8 ) result
  real ( kind = 8 ) sumf
  real ( kind = 8 ) sumg
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xleast
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xx
!
!  Machine-dependent constants
!
  data xleast /2.23d-308/
  data xsmall /1.11d-16/
  data xinf /1.79d+308/
  data xmax /705.343d+0/
!
!  Coefficients for  XLEAST <=  ARG  <= 1.0
!
  data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1, &
            7.1885382604084798576d+3, 1.7733324035147015630d+5, &
            7.1938920065420586101d+5/
  data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4, &
           -2.2149374878243304548d+6/
  data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1, &
           -4.5051623763436087023d+3,-1.4758069205414222471d+5, &
           -1.3531161492785421328d+6/
  data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4, &
           -2.7062322985570842656d+6/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0, &
            1.3182609918569941308d+2, 8.1094256146537402173d+2, &
            2.3123742209168871550d+3, 3.4540675585544584407d+3, &
            2.8590657697910288226d+3, 1.3319486433183221990d+3, &
            3.4122953486801312910d+2, 4.4137176114230414036d+1, &
            2.2196792496874548962d+0/
  data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2, &
            1.2082692316002348638d+3, 2.1181000487171943810d+3, &
            1.9448440788918006154d+3, 9.6929165726802648634d+2, &
            2.5951223655579051357d+2, 3.4552228452758912848d+1, &
            1.7710478032601086579d+0/

  x = arg
!
!  Error return for ARG < XLEAST.
!
  if ( x < xleast ) then

    result = xinf
!
!  XLEAST <= ARG <= 1.0.
!
  else if ( x <= 1.0D+00 ) then

    if ( x < xsmall ) then
!
!  Return for small ARG.
!
      result = 1.0D+00 / x

    else

      xx = x * x

      sump = (((( &
               p(1) &
        * xx + p(2) ) &
        * xx + p(3) ) &
        * xx + p(4) ) &
        * xx + p(5) ) &
        * xx + q(3)

      sumq = (( &
          xx + q(1) ) &
        * xx + q(2) ) &
        * xx + q(3)

      sumf = ((( &
               f(1) &
        * xx + f(2) ) &
        * xx + f(3) ) &
        * xx + f(4) ) &
        * xx + f(5)

      sumg = (( &
          xx + g(1) ) &
        * xx + g(2) ) &
        * xx + g(3)

      result = ( xx * log ( x ) * sumf / sumg + sump / sumq ) / x

      if ( jint == 2 ) then
        result = result * exp ( x )
      end if

    end if

  else if ( jint == 1 .and. xmax < x ) then
!
!  Error return for XMAX < ARG.
!
    result = 0.0D+00

  else
!
!  1.0 < ARG.
!
    xx = 1.0D+00 / x

    sump = pp(1)
    do i = 2, 11
      sump = sump * xx + pp(i)
    end do

    sumq = xx
    do i = 1, 8
      sumq = ( sumq + qq(i) ) * xx
    end do
    sumq = sumq + qq(9)

    result = sump / sumq / sqrt ( x )

    if ( jint == 1 ) then
      result = result * exp ( -x )
    end if

  end if

  return
end

