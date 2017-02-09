!----------------------------!
!   Dielectric function      !
!----------------------------!

module mod_dielectric_function

    ! dielectric function \epsilon(q)
    complex(8), allocatable :: epsilon(:,:,:)
    
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
          allocate(epsh(iomstart:iomend,3,3))
          epsh(:,:,:) = 0.d0    
          if (allocated(epsw1)) deallocate(epsw1)
          allocate(epsw1(mbsiz,iomstart:iomend,3))
          epsw1(:,:,:) = 0.d0
          if (allocated(epsw2)) deallocate(epsw2)
          allocate(epsw2(mbsiz,iomstart:iomend,3))
          epsw2(:,:,:) = 0.d0
          ! macroscopic dielectric tensor
          if (allocated(eps00)) deallocate(eps00)
          allocate(eps00(iomstart:iomend,3,3))
          eps00(:,:,:) = 0.d0
        end if ! Gamma
        ! Second order screened potential W^{(2)}
        if (input%gw%selfenergy%secordw) then
          if (allocated(vPv)) deallocate(vPv)
          allocate(vPv(mbsiz,mbsiz,iomstart:iomend))
          if (Gamma) then
            if (allocated(vPvh)) deallocate(vPvh)
            allocate(vPvh(iomstart:iomend))
            vPvh(:) = 0.d0    
            if (allocated(vPvw1)) deallocate(vPvw1)
            allocate(vPvw1(mbsiz,iomstart:iomend))
            vPvw1(:,:) = 0.d0
            if (allocated(vPvw2)) deallocate(vPvw2)
            allocate(vPvw2(mbsiz,iomstart:iomend))
            vPvw2(:,:) = 0.d0
          end if
        end if
    
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
        if (input%gw%selfenergy%secordw) then
          if (allocated(vPv)) deallocate(vPv)
          if (Gamma) then
            if (allocated(vPvh)) deallocate(vPvh)
            if (allocated(vPvw1)) deallocate(vPvw1)
            if (allocated(vPvw2)) deallocate(vPvw2)
          end if
        end if
    end subroutine
                                         
end module
