
module modgw

    use mod_core_states
    use mod_product_basis
    use mod_gaunt_coefficients
    use mod_vxc
    use mod_misc_gw
    use mod_bands
    use mod_kpointset
    use mod_frequency
    use mod_coulomb_potential
    use mod_dielectric_function
    use mod_selfenergy


    use mod_kqpts       ! original definitions, should be completely replaced by 
                        ! mod_kpointset
    
    !use mod_bzintw

    ! general gw output file
    integer(4) :: fgw
    ! debug info
    integer(4) :: fdebug
    ! gw results file name
    character(128) :: fgwh5
    character(128) :: path, cik

    !-----------------------------------
    ! variables for e-ph coupling calc  
    !-----------------------------------
    integer(4) :: ibeph 
    integer(4) :: nbeph 
    integer(4) :: ibsumeph 
    integer(4) :: nbsumeph 
    integer(4) :: nomegeph 
    integer(4) :: ngridkqtot
    real   (8) :: efnew, cbm
    real   (8), allocatable :: g2eph (:) 

    ! Lower band index for GW output
    integer(4) :: ibgw
    ! Upper band index for GW output
    integer(4) :: nbgw
    ! Number of bands for gw output      
    integer(4) :: nbandsgw
    ! Number of electrons used in GW
    real(8)    :: nvelgw   
    
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
    type(Gk_set) :: Gkset, Gqset, Gqbarc
     
    !--------------------------------!
    ! BZ integration weights
    !--------------------------------!
    real(8), allocatable :: ciw(:,:)
    real(8), allocatable :: kiw(:,:)
    real(8), allocatable :: kwfer(:,:)
    
    complex(8), allocatable :: fnm(:,:,:,:)
    complex(8), allocatable :: unw(:,:,:,:,:)   
    complex(8), allocatable :: kcw(:,:,:,:)
    
    
    !------------------------------------
    ! Memory size of large global arrays
    !------------------------------------
    real(8) :: msize
    ! byte to Mbyte converter factor 
    real(8), parameter :: b2mb = 1.d0/1024/1024
    
    !--------------------
    ! Matrix block size 
    !--------------------
    integer(4) :: mblksiz
    
    !---------!
    ! Timing
    !---------!
    real(8) :: time_total
    real(8) :: time_io
    
    real(8) :: time_initgw
    real(8) :: time_initscf
    real(8) :: time_initkpt
    real(8) :: time_initfreq
    real(8) :: time_initeval
    real(8) :: time_initmb

    real(8) :: time_pmat
    real(8) :: time_barcmb
    real(8) :: time_vxc
    real(8) :: time_bzinit
    
    real(8) :: time_df
    real(8) :: time_dfhead
    real(8) :: time_dfwing
    real(8) :: time_dfbody
    real(8) :: time_dfinv
    
    real(8) :: time_eprod
    real(8) :: time_diagsgi
    real(8) :: time_mpwipw
    real(8) :: time_micm
    real(8) :: time_minc
    real(8) :: time_minm
    
    real(8) :: time_selfx
    real(8) :: time_selfc
    
    real(8) :: time_rotmb
      
contains

!-------------------------------------------------------------------------------      

    subroutine init_timing
        implicit none
        time_total = 0.d0
        
        time_io = 0.d0
        
        time_initgw = 0.d0
        time_initscf = 0.d0
        time_initkpt = 0.d0
        time_initeval = 0.d0
        time_initfreq = 0.d0
        time_initmb = 0.d0
        
        time_pmat = 0.d0
        time_barcmb = 0.d0
        time_vxc = 0.d0
        time_bzinit = 0.d0
        
        time_df = 0.d0
        time_dfhead = 0.d0
        time_dfwing = 0.d0
        time_dfbody = 0.d0
        time_dfinv = 0.d0
        
        time_eprod = 0.d0
        time_diagsgi = 0.d0
        time_mpwipw = 0.d0
        time_micm = 0.d0
        time_minc = 0.d0
        time_minm = 0.d0
        
        time_selfx = 0.d0
        time_selfc = 0.d0
        
        time_rotmb = 0.d0
        
    end subroutine

!-------------------------------------------------------------------------------      
    subroutine print_timing
      use modmpi, only: rank
      implicit none
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
        
        if (input%gw%reduceq) then
          write(fgw, '("     - genmbrotmat", T45,": ", F12.2)') time_rotmb
        end if
        
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

