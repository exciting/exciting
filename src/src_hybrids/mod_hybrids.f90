
module mod_hybrids

    use modinput
    use modmain
    use modgw
    use modmpi, only: rank
    implicit none

! number of HF cycles
    integer :: ihyb

! non-local exchange energy
    real(8) :: exnl

! non-local exchange potential
    complex(8), allocatable :: vxnl(:,:,:)

! APW matrix elements of the non-local potential
    complex(8), allocatable :: vnlmat(:,:,:)

!*******************************************************************************
contains

    subroutine init_hybrids
        implicit none
        integer(4) :: i, ia, is, ias, ic, il, ist, m
        integer(4) :: ir, l, ie1, ikp
        real(8)    :: norm

! print debugging information    
        debug = .false.

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
        if (rank .Eq. 0) call boxmsg(fgw,'=','Mixed Product Basis')
!---------------------------------------
! treatment of core electrons
!---------------------------------------
        select case (input%gw%coreflag)
            case('all','ALL')
                iopcore=0
                if (rank .Eq. 0) write(fgw,*)' all: All electron calculation'
            case('xal','XAL')
                iopcore=1
                if (rank .Eq. 0) write(fgw,*)' xal: all electron for exchange, valence only for correlation'
            case('val')
                iopcore=2
                if (rank .Eq. 0) write(fgw,*)' val: Valence electrons only'
            case('vab')
                iopcore=3
                if (rank .Eq. 0) write(fgw,*)' vab: Valence-only without core states in mixbasis'
        end select         
        if (rank .Eq. 0) then
!---------------------------------------
! MB related parameters
!---------------------------------------
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

!---------------------------------------
! Determine the number of core states for each species (auxiliary arrays)
!---------------------------------------    
        if (allocated(ncore)) deallocate(ncore)
        allocate(ncore(nspecies))
        ncmax=0
        nclm=0
        ncg=0
        lcoremax=0
        do is=1,nspecies
            ncore(is)=0
            ic = 0
            do ist=1,spnst(is)
                if (spcore(ist,is)) then
                    ncore(is)=ncore(is)+1
                    il=spl(ist,is)
                    do m=-spk(ist,is),spk(ist,is)-1
                        ic=ic+1
                    end do
                end if
            end do
            ncmax=max(ncmax,ncore(is))
            nclm=max(nclm,ic)
            lcoremax=max(lcoremax,il)
            ncg=ncg+ic*natoms(is)
        end do
! setting a unique index for all the core states of all atoms
        call setcorind

! reciprocal cell volume
        vi=1.0d0/omega

! shortcut for basis vectors 
        avec(:,1)=input%structure%crystal%basevect(:,1)
        avec(:,2)=input%structure%crystal%basevect(:,2)
        avec(:,3)=input%structure%crystal%basevect(:,3)

! reciprocal lattice basis lengths
        do i=1,3
            alat(i)=dsqrt(avec(1,i)*avec(1,i)+avec(2,i)*avec(2,i)+avec(3,i)*avec(3,i))
            pia(i)=2.0d0*pi/alat(i)
        end do

! additional arrays used simply for convenience
        do is=1,nspecies
            do ia=1,natoms(is)
                ias=idxas(ia,is)
! shortcut for atomic positions
                atposl(:,ia,is)=input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            end do
! calculate the muffin-tin volume
            vmt(is)=4.0d0*pi*rmt(is)*rmt(is)*rmt(is)/(3.0d0*omega)
        end do

        if (allocated(ipwint)) deallocate(ipwint)
        allocate(ipwint(ngrtot))
        ipwint(:) = conjg(cfunig(:))

        return
    end subroutine

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! deallocate hybrids related data
    subroutine exit_hybrids

! deallocate global
        deallocate(vxnl)
        deallocate(vnlmat)

! deallocate mixed-basis stuff
        deallocate(ncore)
        deallocate(corind)
        deallocate(ipwint)
        deallocate(nup)
        deallocate(nmix)
        deallocate(umix)
        deallocate(bigl)
        deallocate(mbl)
        deallocate(rtl)
        deallocate(rrint)
        deallocate(bradketc)
        deallocate(bradketa)
        deallocate(bradketlo)
        deallocate(locmixind)
      
! avoid running gw
        if (associated(input%gw)) input%gw%taskname="skip"

    end subroutine

end module
