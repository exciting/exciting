      SUBROUTINE gga_x_wpbeh(n,RHO,GRHO,sx,V1X, OMEGA)
!      SUBROUTINE gga_x_wpbeh(RHO,GRHO,sx,V1X,V2X,OMEGA)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER, INTENT(in) :: n
      REAL*8, INTENT(in) :: rho(n)
      REAL*8, INTENT(in) :: grho(n)
      REAL*8, INTENT(out) :: sx(n)
      REAL*8, INTENT(out) :: v1x(n)
      REAL*8 :: v2x(n)
      REAL*8, INTENT(in) :: omega
! NOTE, here sx is the total energy density,
! not just the gradient correction energy density as e.g. in pbex()
! And the same goes for the potentials V1X, V2X

      real*8, PARAMETER :: SMALL=1.D-20,SMAL2=1.D-08
      real*8, PARAMETER :: US=0.161620459673995492D0,AX=-0.738558766382022406D0, &
                UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK
      REAL*8, PARAMETER :: f1 = -1.10783814957303361d0, alpha = 2.0d0/3.0d0
      !local variables
      integer i    
      REAL*8 :: RS, VX, AA, RR, EX, S2, S, DSDN, DSDG, FX, D1X, D2X
!     ==--------------------------------------------------------------==

sx=0.0d0
v1x=0.0d0
v2x=0.0d0
      do i=1,n
!!      CALL XC(RHO,EX,EC,VX,VC)
        RS = RHO(i)**(1.0d0/3.0d0)
        VX = (4.0d0/3.0d0)*f1*alpha*RS

!!      AA    = DMAX1(GRHO,SMAL2)
        AA    = GRHO(i)
!!      RR    = RHO**(-4.0_DP/3.0_DP)
        RR    = 1.0d0/(RHO(i)*RS)
        EX    = AX/RR  !rho(i)*E_x^{LDA}
        S2    = AA*RR*RR*US*US !US=1/(2*(3*pi**2)**(1/3))

        S = SQRT(S2)
        IF(S.GT.8.3d0) THEN
          S = 8.572844D0 - 18.796223D0/S2
        ENDIF
        CALL wpbe_analy_erfc_approx_grad(RHO(i),S,OMEGA,FX,D1X,D2X)
        !sx(i) = EX*FX        ! qunatum esprsso
        sx(i) = (EX*FX)/rho(i)        ! real energy density
        DSDN = -4.D0/3.D0*S/RHO(i)
        !D1X = D FX / D r_alpha, D2X = D FX / D s
        V1X(i) = VX*FX + (DSDN*D2X+D1X)*EX   ! - VX
        write(1002,*) "sr: i, sx, v1x", i, sx(i), v1x(i) 
        call flushifc(1002)
        DSDG = US*RR
        V2X(i) = EX*1.D0/SQRT(AA)*DSDG*D2X  
      enddo

!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE gga_x_wpbeh

!      subroutine gcxc(rho, grho, sx, v1x, v2x)
!      !-----------------------------------------------------------------------
!      !     gradient corrections for exchange and correlation - Hartree a.u.
!      !     See comments at the beginning of module for implemented cases
!      !
!      !     input:  rho, grho=|\nabla rho|^2
!      !     definition:  E_x = \int E_x(rho,grho) dr
!      !     output: sx = E_x(rho,grho)
!      !             v1x= D(E_x)/D(rho)
!      !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
!      !
!      implicit none
!
!      real*8 :: rho, grho, sx, v1x, v2x
!      real*8 :: sx__,v1x__, v2x__
!      real*8 :: sxsr, v1xsr, v2xsr
!      real*8 , parameter:: small = 1.E-10d0
!
!      ! exchange
!      if (rho <= small) then
!         sx = 0.0d0
!         v1x = 0.0d0
!         v2x = 0.0d0
!      else  ! 'pbexsr'
!     call pbex (rho, grho, 1, sx, v1x, v2x)
!     if(exx_started) then
!         call pbexsr (rho, grho, sxsr, v1xsr, v2xsr, screening_parameter)
!       sx = sx - exx_fraction * sxsr
!       v1x = v1x - exx_fraction * v1xsr
!       v2x = v2x - exx_fraction * v2xsr
!      endif
!
!      return
!      end subroutine gxcx

