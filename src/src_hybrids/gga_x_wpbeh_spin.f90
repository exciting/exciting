      SUBROUTINE gga_x_wpbeh_spin(n,RHOin,GRHOin,sx,V1X, v2x,OMEGA)
!      SUBROUTINE gga_x_wpbeh(RHO,GRHO,sx,V1X,V2X,OMEGA)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER, INTENT(in) :: n
      REAL*8, INTENT(in) :: rhoin(n)
      REAL*8, INTENT(in) :: grhoin(n)
      REAL*8, allocatable :: rho(:)
      REAL*8, allocatable :: grho(:)
      REAL*8, INTENT(out) :: sx(n)
      REAL*8, INTENT(out) :: v1x(n)
      REAL*8, intent(out) :: v2x(n)
      REAL*8, INTENT(in) :: omega
! NOTE, here sx is the total energy density,
! not just the gradient correction energy density as e.g. in pbex()
! And the same goes for the potentials V1X, V2X
      real*8, parameter :: tolerance = 1.D-14
!> Tolerance below which the charge density (rho) is considered zero 
!> This choice of cutoff was determined empirically and avoids divergence
!of the quantities that are divided by the density for 1D and 2D materials 
      !real*8, parameter :: SMAL2 = 1.D-08

      real*8, PARAMETER :: US=0.161620459673995492D0,AX=-0.738558766382022406D0, &
                UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK
      REAL*8, PARAMETER :: f1 = -1.10783814957303361d0, alpha = 2.0d0/3.0d0
      !local variables
      integer i    
      REAL*8 :: RS, VX, AA, RR, EX, S2, S, DSDN, DSDG, FX, D1X, D2X
!     ==--------------------------------------------------------------==
      If (allocated(rho)) deallocate(rho)
      Allocate(rho(n))
      If (allocated(grho)) deallocate(grho)
      Allocate(grho(n))
      sx(:) = 0.0d0
      v1x(:) = 0.0d0
      v2x(:) = 0.0d0
        rho(:)=2.d0*rhoin(:)
        grho(:)=2.d0*grhoin(:) 
      do i=1,n
        if (rho(i)<tolerance) cycle
!!      CALL XC(RHO,EX,EC,VX,VC)
        RS = RHO(i)**(1.0d0/3.0d0)
        VX = (4.0d0/3.0d0)*f1*alpha*RS

!!      AA    = DMAX1(GRHO,SMAL2)
!!        AA    = GRHO(i)
        AA    = GRHO(i)**2
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
        DSDG = US*RR
        V2X(i) = EX*1.D0/SQRT(AA)*DSDG*D2X  
        !V2X(i) = EX*DSDG*D2X  
      enddo

!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE gga_x_wpbeh_spin
