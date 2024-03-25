!----------------------------!
!   Dielectric function      !
!----------------------------!

module mod_dielectric_function
    use precision, only: wp

    ! dielectric function \epsilon(q)
    complex(8), allocatable :: epsilon(:,:,:)
    !
    !----------------------- Manoar Hossain ----------------------------------------------
    ! Variables for cubic GW implementation
    !> Green's function in R and MT and for occupied part
    !!complex(wp), allocatable :: greenr_pp_occ(:,:,:,:,:)
    !> Green's function in R and MT and for unoccupied part
    !!complex(wp), allocatable :: greenr_pp_uno(:,:,:,:,:)
    !> Polarizability in MT-MT region
    !!complex(wp), allocatable :: polarizabilityR_mtmt(:,:,:,:,:)
    !> Polarizability in q-space in MT-MT region
    complex(wp), allocatable :: pola_q_mtmt(:,:,:,:,:)
    !> Polarizability in q-space in MT-I region
    complex(wp), allocatable :: pola_q_mti(:,:,:,:)
    !> Polarizability in q-space in MT-I region
    complex(wp), allocatable :: pola_q_ii(:,:,:)
    !------------------------------------------------------------------------

    !-------------------------------------------------
    ! Analytical treatment of q=0 singularity
    !-------------------------------------------------

    ! valence-valence momentum matrix elements
    complex(8), allocatable :: pmatvv(:,:,:)

    ! core-valence momentum matrix elements
    complex(8), allocatable :: pmatcv(:,:,:)

    ! head of the dielectric function (tensor)
    complex(8), allocatable :: epsh(:,:,:)

    ! the vertical wing of the dielectric matrix (vector)
    complex(8), allocatable :: epsw1(:,:,:)

    ! the horizontal wing of the dielectric matrix (vector)
    complex(8), allocatable :: epsw2(:,:,:)

    !----------------------------------------------------------------------
    ! Used for calculating the macroscopic dielectric function (task_emac)
    !----------------------------------------------------------------------

    complex(8), allocatable :: eps00(:,:,:)

    !--------------------------------------
    complex(8), allocatable :: vPv(:,:,:)
    complex(8), allocatable :: vPvh(:)
    complex(8), allocatable :: vPvw1(:,:)
    complex(8), allocatable :: vPvw2(:,:)

    !----------------------------------------------------------------------
    ! files containing data on PMAT and PMATCOR
    !----------------------------------------------------------------------
    integer :: fid_pmatvv=300
    character(24) :: fname_pmatvv='PMATVV.OUT'

    integer :: fid_pmatcv=310
    character(24) :: fname_pmatcv='PMATCV.OUT'

    !----------------------------------------------------------------------
    ! file to store the dielectric function
    !----------------------------------------------------------------------
    integer :: fid_eps
    character(24) :: fname_eps='EPSILON'
    character(24) :: fname_head='HEAD.OUT'
    character(24) :: fname_wings='WINGS.OUT'

contains

    subroutine init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
        use modinput
        implicit none
        integer, intent(in) :: mbsiz
        integer, intent(in) :: iomstart, iomend
        logical, intent(in) :: Gamma
        ! q-dependent dielectric function
        if (allocated(epsilon)) deallocate(epsilon)
        allocate(epsilon(mbsiz,mbsiz,iomstart:iomend))
        epsilon(:,:,:) = 0.d0
        ! head and wings of the dielectric function when q->0
        if (Gamma) then
          if (allocated(epsh)) deallocate(epsh)
          allocate(epsh(3,3,iomstart:iomend))
          epsh(:,:,:) = 0.d0
          if (allocated(epsw1)) deallocate(epsw1)
          allocate(epsw1(mbsiz,3,iomstart:iomend))
          epsw1(:,:,:) = 0.d0
          if (allocated(epsw2)) deallocate(epsw2)
          allocate(epsw2(mbsiz,3,iomstart:iomend))
          epsw2(:,:,:) = 0.d0
          ! macroscopic dielectric tensor
          if (allocated(eps00)) deallocate(eps00)
          allocate(eps00(3,3,iomstart:iomend))
          eps00(:,:,:) = 0.d0
        end if ! Gamma
    end subroutine

    subroutine delete_dielectric_function(Gamma)
        use modinput
        implicit none
        logical, intent(in) :: Gamma
        if (allocated(epsilon)) deallocate(epsilon)
        if (Gamma) then
          if (allocated(epsh)) deallocate(epsh)
          if (allocated(epsw1)) deallocate(epsw1)
          if (allocated(epsw2)) deallocate(epsw2)
          if (allocated(eps00)) deallocate(eps00)
        end if
    end subroutine

end module
