
module modgw
    use gw_info, only: fgw
    use mod_core_states
    use mod_dielectric_function
    use mod_frequency
    use mod_gaunt_coefficients
    use mod_kpointset
    use mod_kqpts       ! original definitions, should be completely replaced by mod_kpointset
    use mod_misc_gw
    use mod_product_basis
    use mod_selfenergy
    use precision, only: dp, i32, str_128
    
    implicit none

    ! debug info
    integer(i32) :: fdebug
    ! gw results file name
    character(str_128) :: fgwh5
    character(str_128) :: path, cik

    !-----------------------------------
    ! variables for e-ph coupling calc  
    !-----------------------------------
    integer(i32) :: ibeph 
    integer(i32) :: nbeph 
    integer(i32) :: ibsumeph 
    integer(i32) :: nbsumeph 
    integer(i32) :: nomegeph 
    integer(i32) :: ngridkqtot
    real(dp) :: efnew, cbm
    real(dp), allocatable :: g2eph (:) 

    ! Lower band index for GW output
    integer(i32) :: ibgw
    ! Upper band index for GW output
    integer(i32) :: nbgw
    ! Number of bands for gw output      
    integer(i32) :: nbandsgw
    ! Number of electrons used in GW
    real(dp)     :: nvelgw   
    
    !----------------------------!
    ! frequency grid         !
    !----------------------------!
    type(frequency):: freq
    
    !--------------------------------!
    ! reciprocal space variables !
    !--------------------------------!
    type(k_set)  :: kset
    type(k_set)  :: qsetd
    type(kkqmt_set) :: kqsetd
    type(kq_set) :: kqset
    type(G_set)  :: Gset
    type(Gk_set) :: Gkset, Gkqset, Gqset, Gqbarc
     
    !--------------------------------!
    ! BZ integration weights
    !--------------------------------!
    real(dp), allocatable :: ciw(:, :)
    real(dp), allocatable :: kiw(:, :)
    real(dp), allocatable :: kwfer(:, :)
    
    complex(dp), allocatable, target :: fnm_sum(:, :, :) 
    complex(dp), allocatable, target :: fnm_tet(:, :, :, :)
    complex(dp), pointer, contiguous :: fnm(:, :, :)
    complex(dp), allocatable :: unw(:, :, :, :, :)   
    complex(dp), allocatable :: kcw(:, :, :, :)
    
    
    !------------------------------------
    ! Memory size of large global arrays
    !------------------------------------
    real(dp) :: msize
    ! byte to Mbyte converter factor 
    real(dp), parameter :: b2mb = 1._dp/1024/1024
    
    !--------------------
    ! Matrix block size 
    !--------------------
    integer(i32) :: mblksiz
    
    !---------!
    ! Timing
    !---------!
    real(dp) :: time_total
    real(dp) :: time_io
    
    real(dp) :: time_initgw
    real(dp) :: time_initscf
    real(dp) :: time_initkpt
    real(dp) :: time_initfreq
    real(dp) :: time_initeval
    real(dp) :: time_initmb

    real(dp) :: time_pmat
    real(dp) :: time_barcmb
    real(dp) :: time_vxc
    real(dp) :: time_bzinit
    
    real(dp) :: time_df
    real(dp) :: time_dfhead
    real(dp) :: time_dfwing
    real(dp) :: time_dfbody
    real(dp) :: time_dfinv
    
    real(dp) :: time_eprod
    real(dp) :: time_diagsgi
    real(dp) :: time_mpwipw
    real(dp) :: time_micm
    real(dp) :: time_minc
    real(dp) :: time_minm
    
    real(dp) :: time_selfx
    real(dp) :: time_selfc
    
    real(dp) :: time_rotmb
      
contains

!-------------------------------------------------------------------------------      

    subroutine init_timing()
        time_total = 0._dp
        
        time_io = 0._dp
        
        time_initgw = 0._dp
        time_initscf = 0._dp
        time_initkpt = 0._dp
        time_initeval = 0._dp
        time_initfreq = 0._dp
        time_initmb = 0._dp
        
        time_pmat = 0._dp
        time_barcmb = 0._dp
        time_vxc = 0._dp
        time_bzinit = 0._dp
        
        time_df = 0._dp
        time_dfhead = 0._dp
        time_dfwing = 0._dp
        time_dfbody = 0._dp
        time_dfinv = 0._dp
        
        time_eprod = 0._dp
        time_diagsgi = 0._dp
        time_mpwipw = 0._dp
        time_micm = 0._dp
        time_minc = 0._dp
        time_minm = 0._dp
        
        time_selfx = 0._dp
        time_selfc = 0._dp
        
        time_rotmb = 0._dp
        
    end subroutine

!-------------------------------------------------------------------------------      
    subroutine print_timing
      use modmpi, only: rank

      if (rank==0) then
        call boxmsg(fgw,'=','GW timing info (seconds)') 

        write(fgw, '(" Initialization", T45, ": ", F12.2)') time_initgw
        write(fgw, '("     - init_scf", T45,": ", F12.2)') time_initscf
        write(fgw, '("     - init_kpt", T45,": ", F12.2)') time_initkpt
        write(fgw, '("     - init_eval", T45,": ", F12.2)') time_initeval
        write(fgw, '("     - init_freq", T45,": ", F12.2)') time_initfreq
        write(fgw, '("     - init_mb", T45,": ", F12.2)') time_initmb

        write(fgw, '(" Subroutines", T45, ": ", F12.2)')
        
        write(fgw, '("     - calcpmat", T45,": ", F12.2)') time_pmat
        write(fgw, '("     - calcbarcmb", T45,": ", F12.2)') time_barcmb
        write(fgw, '("     - BZ integration weights", T45,": ", F12.2)') time_bzinit
        
        write(fgw, '("     Dielectric function", T45,": ", F12.2)') time_df
        write(fgw, '("     - head", T45,": ", F12.2)') time_dfhead
        write(fgw, '("     - wings", T45,": ", F12.2)') time_dfwing
        write(fgw, '("     - body", T45,": ", F12.2)') time_dfbody
        write(fgw, '("     - inversion", T45,": ", F12.2)') time_dfinv
        
        write(fgw, '("     WF products expansion", T45,": ", F12.2)') &
        &  time_diagsgi+time_mpwipw+time_eprod
        write(fgw, '("     - diagsgi", T45,": ", F12.2)') time_diagsgi
        write(fgw, '("     - calcmpwipw", T45,": ", F12.2)') time_mpwipw
        write(fgw, '("     - calcmicm", T45,": ", F12.2)') time_micm
        write(fgw, '("     - calcminc", T45,": ", F12.2)') time_minc
        write(fgw, '("     - calcminm", T45,": ", F12.2)') time_minm
        
        write(fgw, '("     - symmetry", T45,": ", F12.2)') time_rotmb
        
        write(fgw, '("     Self-energy", T45,": ", F12.2)') &
        &  time_selfx+time_selfc
        write(fgw, '("     - calcselfx", T45,": ", F12.2)') time_selfx
        write(fgw, '("     - calcselfc", T45,": ", F12.2)') time_selfc
        
        write(fgw, '("     - calcvxcnn", T45,": ", F12.2)') time_vxc
       
        write(fgw, '("     - input/output", T45,": ", F12.2)') time_io
        write(fgw,*)'_________________________________________________________'
        write(fgw, '(" Total", T45, ": ", F12.2)') time_total
        write(fgw,*)
      end if        
    end subroutine
    
end module modgw

