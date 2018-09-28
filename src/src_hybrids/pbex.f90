subroutine pbex (rho, grho, iflag, sx, v1x, v2x)
  !---------------------------------------------------------------
  !
  ! PBE exchange (without Slater exchange):
  ! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  ! iflag=2  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
  ! iflag=3  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
  ! iflag=4  PBEQ2D: L. Chiodo et al., PRL 108, 126402 (2012)
  ! iflag=5  optB88: Klimes et al., J. Phys. Cond. Matter, 22, 022201 (2010)
  ! iflag=6  optB86b: Klimes et al., Phys. Rev. B 83, 195131 (2011)
  ! iflag=7  ev: Engel and Vosko, PRB 47, 13164 (1991)
  !
  USE kinds, ONLY : DP
    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
    Use mod_hybrids
!  USE constants, ONLY : pi
#if defined(__LIBXC)
  use xc_f90_types_m
  use xc_f90_lib_m
#endif
  implicit none
  real(dp), intent(in) :: rho, grho
  ! input: charge and squared gradient
  real(dp), intent(out):: sx, v1x, v2x
  ! output: energy, potential
  integer, intent(in) :: iflag
!  REAL(DP), PARAMETER :: pi= 3.14159265358979323846_DP
#if defined(__LIBXC) 
  ! local variables
  integer :: func_id = -1 ! not set
  integer :: size = 1
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  real(dp) :: exc, ex_lda = 0.0d0 , vx_lda = 0.0d0
   
  
  if (iflag.eq.1)  func_id = 101
  if (iflag.eq.2)  func_id = 102
  if (iflag.eq.3)  func_id = 116
  if (iflag.eq.5)  func_id = 141
  if (func_id==-1) call errore('pbex','case not implemented with libxc',iflag)

  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)    
  call xc_f90_gga_exc_vxc(xc_func, size, rho, grho, exc, v1x, v2x)  
  call xc_f90_func_end(xc_func)  

  ! remove Slater term for compatibility with QE  
  call xc_f90_func_init(xc_func, xc_info, 1, XC_UNPOLARIZED)       
  call xc_f90_lda_exc_vxc(xc_func, size, rho, ex_lda, vx_lda)  
  call xc_f90_func_end(xc_func) 
  exc = exc - ex_lda 
  v1x = v1x - vx_lda 
  
  !sx = exc * rho  ! e_x = rho * \epsilon_x
  sx = exc   ! e_x = rho * \epsilon_x
  !v2x = v2x*2.0_dp
  v2x = v2x*2.0_dp*sqrt(grho)

#else
  ! local variables
  real(DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(DP) :: dxunif, dfx, f1, f2, f3, dfx1
  real(DP) :: p,  amu, ab, c, dfxdp, dfxds, upbe, uge, s, ak, aa
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(DP), parameter :: third = 1._DP / 3._DP, c1 = 0.75_DP / pi , &
       c2 = 3.093667726280136_DP, c5 = 4._DP * third, &
       c6 = c2*2.51984210, c7=5._DP/6._DP, c8=0.8_DP ! (3pi^2)^(1/3)*2^(4/3)
  ! parameters of the functional
  real(DP) :: k (6), mu(6), ev(6)
  !           pbe        rpbe        pbesol   pbeq2d      optB88  optB86b
  data k / 0.804_DP,   1.2450_DP,   0.804_DP , 0.804_DP,  0.0_dp,  0.0_dp/, &
       mu/ 0.2195149727645171_DP, 0.2195149727645171_DP, 0.12345679012345679_DP, &
           0.12345679012345679_DP,  0.22_dp, 0.1234_dp/, &
       ev / 1.647127_DP, 0.980118_DP, 0.017399_DP, 1.523671_DP, 0.367229_DP, &
                                   0.011282_DP /  ! a and b parameters of Engel and Vosko
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5_DP / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  !   Energy
  !
  if ( iflag == 4) then
     p=s1*s1
     s=s1
     ak=0.804_DP
     amu=10._DP/81._DP
     ab=0.5217_DP
     c=2._DP
     fx =  ak - ak / (1.0_dp + amu * p / ak)  + p**2 * (1 + p) &
           /(10**c + p**3) * (-1.0_dp - ak + ak / (1.0_dp + amu * p / ak) &
           + ab * p ** (-0.1d1/ 0.4D1))
  elseif ( iflag == 5) then
     ab=mu(iflag)*c7 ! mu/ab=1.2
     p=s1*c6
     c=log(p+sqrt(p*p+1)) ! asinh(p)
     dfx1=1+ab*s1*c
     fx =  mu(iflag)*s1*s1/dfx1
  elseif ( iflag == 6) then
     p=mu(iflag)*s1*s1
     fx =  p / ( 1 + p )**c8
  elseif ( iflag == 7) then
     s=s2*s2
     f1 =  1 + ev(1)*s2 + ev(2)*s + ev(3)*s*s2
     f2 =  1 + ev(4)*s2 + ev(5)*s + ev(6)*s*s2
     fx = f1 / f2 - 1
  else
     f1 = s2 * mu(iflag) / k (iflag)
     f2 = 1._DP + f1
     f3 = k (iflag) / f2
     fx = 1 +  k (iflag) - f3
  end if
  exunif = - c1 * kf
  sx = exunif * fx
  !
  !   Potential
  !
  dxunif = exunif * third
  if ( iflag == 4) then
      dfxdp = dble(1 / (1 + amu * p / ak) ** 2 * amu) + dble(2 * p * (1 &
     + p) / (10 ** c + p ** 3) * (-1 - ak + ak / (1 + amu * p / ak) + ab &
      * p ** (-0.1d1 / 0.4D1))) + dble(p ** 2 / (10 ** c + p ** 3) * ( &
     -1 - ak + ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) - &
      dble(3 * p ** 4 * (1 + p) / (10 ** c + p ** 3) ** 2 * (-1 - ak + &
     ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) + dble(p ** &
      2) * dble(1 + p) / dble(10 ** c + p ** 3) * (-dble(1 / (1 + amu * &
      p / ak) ** 2 * amu) - dble(ab * p ** (-0.5d1 / 0.4D1)) / 0.4D1)

      dfxds=dfxdp*2._DP*s
      dfx=dfxds
  elseif (iflag == 5) then
     dfx=2*fx/s1-fx/dfx1*(ab*c+ab*s1/sqrt(p*p+1)*c6)
  elseif (iflag == 6) then
     dfx=2*mu(iflag)*s1*fx*(1+(1-c8)*p)/(p*(1+p))
  elseif (iflag == 7) then
    dfx  =  ev(1) + 2*ev(2)*s2 + 3*ev(3)*s  
    dfx1 =  ev(4) + 2*ev(5)*s2 + 3*ev(6)*s 
    dfx  = 2 * s1 * ( dfx - f1*dfx1/f2 ) / f2
  else
     dfx1 = f2 * f2
     dfx = 2._DP * mu(iflag) * s1 / dfx1
  end if
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho
  !v2x = exunif * dfx * dsg 
  !sx = sx * rho
  sx = sx 
#endif
  return
end subroutine pbex

