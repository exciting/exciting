module DFT_D2_parameters
  implicit none
  real(8), parameter :: au_to_ang = 0.52917726d0
  real(8), parameter :: J_to_au = 4.35974417d-18
  real(8), parameter :: N_A = 6.02214129d23 !1/mol                                                                                                                                                                                                                             
  real(8), parameter :: damping_const = 20d0
  integer, parameter :: max_elem = 86
  real(8), parameter :: cutoff = 95 !au                                                                                                                                                                                                                                        
  real(8), parameter :: s6 = 0.75 !for pbe: s6 = 0.75   ?                                                                                                                                                                                                                      
  real(8), parameter :: rs6 = 1.1 !for pbe: rs6 = 1.1   ?                                                                                                                                                                                                                     
end module DFT_D2_parameters
