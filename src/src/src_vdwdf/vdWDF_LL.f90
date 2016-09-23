! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!============================================================
subroutine vdWDF_LL(dis,dens1,graddens1,dens2,graddens2,func)
!------------------------------------------------------------
use param
!------------------------------------------------------------

implicit none
real*8  :: dis, phi, func
real*8  :: dens1, dens2
real*8  :: graddens1, graddens2

real*8  :: Zab  
real*8  :: kf1, rs1, ec1, q01,  kf2, rs2, ec2, q02
real*8  :: d1, d2, D, delta

integer :: iD,idelta
real*8  :: D0,delta0
real*8  :: phi11,phi12,phi21,phi22,phi1,phi2

!-----------------------------------------
! Constants
real*8, parameter :: pi34 = 0.6203504908994d0, & ! pi34=(3/4pi)^(1/3)
                     third = 1.d0 / 3.d0,      &
                     c2 = 3.093667726280136d0    ! c2=(3 pi^2)^(1/3)

real*8, parameter :: CvdW = 32.6650486837726d0   ! CvdW = 12 (4*pi/9)^3

!-----------------------------------------
! New: Interpolation for phi in the limit D->0 
real*8, parameter :: p0 = -0.166068,   &
                     p1 =  0.317261,   &
                     p2 = -0.0262737,  &
                     p3 =  0.000927087,&
                     p4 = -6.93206e-06

!----------------------------------------------!
!----------------------------------------------!
!----------------------------------------------!
  
  func = 0.0d0
  if (dis<eps) return

  if (vdWDF_version=='vdW-DF')  Zab = Zab_04
  if (vdWDF_version=='vdW-DF2') Zab = Zab_10
  
! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  rs1 = pi34 / dens1**third
  call pw (rs1, ec1)
  kf1 = c2 * dens1**third
  q01 = kf1 - 4.0d0*pi/3.0d0*ec1 - &
        Zab/(36.0d0*kf1)*(graddens1/dens1)**2.0d0

! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  rs2 = pi34 / dens2**third
  call pw (rs2, ec2)
  kf2 = c2 * dens2**third
  q02 = kf2 - 4.0d0*pi/3.0d0*ec2 - &
        Zab/(36.0d0*kf2)*(graddens2/dens2)**2.0d0

  !write(*,*) q01, q02 
  !write(*,*) q01/kf1, q02/kf2
  
  d1    = dis*q01
  d2    = dis*q02
  D     = 0.5d0*(d1 + d2)
  if ((d1.lt.eps).and.(d2.lt.eps)) then
      delta=0d0
  else
      delta = (d1 - d2)/(d1 + d2)
  end if

  !write(6,*) 'd1   =', d1
  !write(6,*) 'd2   =', d2
  !write(6,*) 'D    =', D
  !write(6,*) 'delta=', delta

  delta = dabs(delta)
  if (delta.lt.deltamin .or. delta.gt.deltamax) then
      write(21,*) 'delta out of range',delta
      stop
  end if 
   
! use empirical relation to determine whether the asymptotic for phi can be used
! (compare the plots for the exact phi and its asymptotic form)
  
  if ( ( (0.167*D - 0.1296*D*delta).gt.1.0d0 ) .or. ( D.gt.Dmax ) ) then
       ! use asymptotic form of Phi(D,delta)     
       phi = - CvdW/( (d1*d2)**2d0*(d1**2d0+d2**2d0) )

  else if (D .le. Dcutoff) then
       ! phi = kernel(idelta,1)
       phi = p0+p1/D+p2/D/D+p3/D/D/D+p4/D/D/D/D

  else
       ! bilinear interpolation between tabulated values
       
       ! to avoid some numerical instabilities while rounding,
       ! an additional small number eps is added
       iD     = 1 + int((D     - Dmin     + eps)/Dstep)
       idelta = 1 + int((delta - deltamin + eps)/deltastep)
      
       D0     = Dmin + dble(iD-1)*Dstep
       delta0 = deltamin + dble(idelta-1)*deltastep
       
       phi11  = kernel(idelta    ,iD)
       phi12  = kernel(idelta + 1,iD)
       phi21  = kernel(idelta    ,iD + 1)
       phi22  = kernel(idelta + 1,iD + 1)                  
       
       phi1   = phi11 + (D - D0)*(phi21 - phi11)/Dstep
       phi2   = phi12 + (D - D0)*(phi22 - phi12)/Dstep     
       
       phi    = phi1  + (delta - delta0)*(phi2 - phi1)/deltastep
       
  end if

!-------------------------------------------
! The integrand valuse (const is predefined in the main program)
  func = 0.5d0*const*dens1*phi*dens2

return
end

!-----------------------------------------------------------------------
subroutine pw (rs, ec)
!-----------------------------------------------------------------------
  implicit none
  real(8) :: rs, ec
  !
  real(8) :: a, a1, b1, b2, b3, b4
  real(8) :: rs12, rs32, rs2, om, dom, olog
  
  parameter (a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = 3.5876d0, &
         b3 = 1.6382d0, b4 = 0.49294d0 )
  
     ! interpolation formula
     rs12 = sqrt (rs)
     rs32 = rs * rs12
     rs2 = rs**2.0d0
     om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
     dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 + &
           2.d0 * b4 * rs2)
     olog = log (1.d0 + 1.0d0 / om)
     ec = - 2.d0 * a * (1.d0 + a1 * rs) * olog
!     vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) &
!          * olog - 2.d0 / 3.d0 * a * (1.d0 + a1 * rs) * dom / &
!          (om * (om + 1.d0) )

  return
end subroutine pw
