
subroutine init_hybrids()
    use modinput
    use modgw
    use modmpi, only : rank
    use mod_mpi_gw
    implicit none
    integer :: lmax, ik

    ! initialize GW MPI environment
    call init_mpi_gw

!---------------------------------------
! MB parameters are taken from GW
!---------------------------------------
    if (.not.associated(input%gw)) &
    &  input%gw => getstructgw(emptynode)

!---------------------------------------
! Options for the mixed basis functions
!---------------------------------------
    if (.not.associated(input%gw%MixBasis)) &
    &  input%gw%MixBasis => getstructmixbasis(emptynode)

!---------------------------------------
! Parameters for the bare coulomb potential
!---------------------------------------
    if (.not.associated(input%gw%BareCoul)) &
    &  input%gw%BareCoul => getstructbarecoul(emptynode)

!---------------------------------------
! print out information on MB to INFO.OUT
!--------------------------------------- 

    ! Hartree-Fock related debugging info
    fgw = 600
    open(fgw, File='HYBRIDS.OUT', Action='WRITE', Form='FORMATTED')
    if (rank == 0) then
        write(fgw,*)
        write(fgw,*) 'Mixed basis parameters:'
        write(fgw,*) '- Interstitial:'
        write(fgw,*) '  -- maximum |G| of IPW in gmaxvr units (gmb):', input%gw%MixBasis%gmb
        write(fgw,*) '- MT-Spheres:'
        write(fgw,*) '  -- l_max (lmaxmb): ', input%gw%MixBasis%lmaxmb
        write(fgw,*) '  -- linear dependence tolerance (epsmb): ', input%gw%MixBasis%epsmb
        write(fgw,*)
        write(fgw,*) 'Bare Coulomb parameters:'
        write(fgw,*) 'Maximum |G| in gmaxvr*gmb units:', input%gw%BareCoul%pwm
        write(fgw,*) 'Error tolerance for struct. const.:', input%gw%BareCoul%stctol
        write(fgw,*) 'Tolerance to choose basis functions from bare Coulomb &
        &  matrix eigenvectors: ', input%gw%BareCoul%barcevtol
        call linmsg(fgw,'=','')
        call flushifc(fgw)
    end if

!---------------------------------------------------------
! Intialize auxiliary arrays used further for convenience    
!---------------------------------------------------------
    call init_misc_gw

!---------------------------------------
!   Initialize k/q grids
!---------------------------------------
    input%gw%ngridq  = input%groundstate%ngridk
    input%gw%vqloff  = input%groundstate%vkloff
    input%gw%reduceq = input%groundstate%reducek

    call init_kqpoint_set
    ! write(*,*) nkpt, kset%nkpt
    ! do ik = 1, nkpt
    !   write(*,'(3f8.4,4x,3f8.4)') vkl(:,ik), kset%vkl(:,ik)
    ! end do

!--------------------------------------------------------------
! Calculate the integrals to treat the singularities at G+q->0
!--------------------------------------------------------------
    call setsingc

! Gaunt coefficients
    lmax = max(input%groundstate%lmaxapw+1, &
    &          2*(input%gw%mixbasis%lmaxmb+1))
    call calcgauntcoef(lmax)

    return
end subroutine
