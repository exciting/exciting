!> Index parameters for XC functionals implemented by exciting (rather than
!> libXC). Also see xml/schema/groundstate.xmd for the definitions. 
module vx_enums
  implicit none
  private

  !> Null entry for when exciting does not use an internal functional
  integer, public, parameter :: no_exciting_xc = 1

  !> LDA, Perdew-Zunger/Ceperley-Alder Phys.Rev.B 23, 5048 (1981)
  integer, public, parameter :: LDA_PZ = 2
  
  !> LSDA, Perdew-Wang/Ceperley-Alder, Phys.Rev.B, 13244 (1992) 
  integer, public, parameter :: LDA_PW = 3

  !> LDA, X-alpha approximation, J. C. Slater, Phys. Rev. 81, 385 (1951) 
  integer, public, parameter :: LDA_XALPHA = 4

  !> LSDA, von Barth-Hedin, J. Phys. C 5, 1629 (1972)
  integer, public, parameter :: LDA_vBH = 5

  !> GGA, Perdew-Burke-Ernzerhof (PBE), Phys. Rev. Lett. 77, 3865 (1996)               
  integer, public, parameter :: GGA_PBE = 20

  !> GGA, Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998) 
  integer, public, parameter :: GGA_PBE_R = 21

  !> GGA, PBEsol, arXiv:0707.2088v1 (2007)
  integer, public, parameter :: GGA_PBE_SOL = 22

  !> GGA, wPBEh (w=0), J. Chem. Phys. 118, 8207 (2003) - 
  !  PBE short-range part used in the hybrid functional HSE. 
  !  Useful for future hybrids implementations.  
  integer, public, parameter :: GGA_PBE_SR = 23

  !> GGA, Wu-Cohen exchange (WC06) with PBE correlation, Phys. Rev. B 73, 235116 (2006) 
  integer, public, parameter :: GGA_WC = 26

  !> GGA, Armiento-Mattsson (AM05) spin-unpolarised functional, Phys. Rev. B 72, 085108 (2005)
  integer, public, parameter :: GGA_AM05 = 30

  !> GGA, asymptotically corrected PBE  (acPBE), arXiv:1409.4834 (2014) 
  integer, public, parameter :: GGA_AC_PBE = 300

  !> Hybrid, PBE0, J. Chem. Phys. 110, 6158, (1999), J. Chem. Phys. 110, 5029 (1999)
  integer, public, parameter :: HYB_PBE0 = 406

  !> Unknown (Implemented for testing by Dima)
  integer, public, parameter :: HYB_LDA0 = 407

  !> Hybrid, HSE06, J. Chem. Phys. 125, 224106 (2006)
  integer, public, parameter :: HYB_HSE = 408
  
end module vx_enums
