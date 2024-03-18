!----------------------------!
!   Dielectric function      !
!----------------------------!

module mod_dielectric_function
    use gw_io, only: build_file_name, write_to_file, read_from_file
    use precision, only: i32

    implicit none

    private

    ! dielectric function \epsilon(q)
    complex(8), public, allocatable :: epsilon(:,:,:)

    !-------------------------------------------------
    ! Analytical treatment of q=0 singularity
    !-------------------------------------------------

    ! valence-valence momentum matrix elements
    complex(8), public, allocatable :: pmatvv(:,:,:)

    ! core-valence momentum matrix elements
    complex(8), public, allocatable :: pmatcv(:,:,:)

    ! head of the dielectric function (tensor)
    complex(8), public, allocatable :: epsh(:,:,:)

    ! the vertical wing of the dielectric matrix (vector)
    complex(8), public, allocatable :: epsw1(:,:,:)

    ! the horizontal wing of the dielectric matrix (vector)
    complex(8), public, allocatable :: epsw2(:,:,:)

    !----------------------------------------------------------------------
    ! Used for calculating the macroscopic dielectric function (task_emac)
    !----------------------------------------------------------------------

    complex(8), public, allocatable :: eps00(:,:,:)

    !----------------------------------------------------------------------
    ! files containing data on PMAT and PMATCOR
    !----------------------------------------------------------------------
    character(24), parameter, public :: fname_pmatvv='PMATVV.OUT'
    character(24), parameter, public :: fname_pmatcv='PMATCV.OUT'
    
    !----------------------------------------------------------------------
    ! files to store the dielectric function
    !----------------------------------------------------------------------
    character(len=*), parameter, private :: file_name_epsilon = 'EPSILON-GW_'
    character(len=*), parameter, private :: file_name_epsilon_head = 'EPSH'
    character(len=*), parameter, private :: file_name_epsilon_wings1 = 'EPSW1'
    character(len=*), parameter, private :: file_name_epsilon_wings2 = 'EPSW2'

    integer(i32), parameter   :: max_string_length = 40

    public :: write_epsilon_to_file, init_dielectric_function, delete_dielectric_function

    
contains

    subroutine init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
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
        logical, intent(in) :: Gamma
        if (allocated(epsilon)) deallocate(epsilon)
        if (Gamma) then
          if (allocated(epsh)) deallocate(epsh)
          if (allocated(epsw1)) deallocate(epsw1)
          if (allocated(epsw2)) deallocate(epsw2)
          if (allocated(eps00)) deallocate(eps00)
        end if
    end subroutine


    subroutine write_epsilon_to_file( iq, is_Gamma_point, file_format )
      integer(i32), intent(in)  :: iq 
      logical, intent(in)       :: is_Gamma_point
      character(len=*), intent(in) :: file_format
    
      character(len=max_string_length) :: file_name
    
      call build_file_name( file_name_epsilon, iq, file_name )
      call write_to_file( file_name, epsilon, file_format )
      if( is_Gamma_point ) then
          call build_file_name( file_name_epsilon_head, file_name )
          call write_to_file( file_name, epsh, file_format )
          call build_file_name( file_name_epsilon_wings1, file_name )
          call write_to_file( file_name, epsw1, file_format )
          call build_file_name( file_name_epsilon_wings2, file_name )
          call write_to_file( file_name, epsw2, file_format )
      end if
      
  end subroutine


  subroutine read_epsilon_from_file( iq, is_Gamma_point, file_format )
      integer(i32), intent(in)  :: iq 
      logical, intent(in)       :: is_Gamma_point
      character(len=*), intent(in) :: file_format
    
      character(len=max_string_length) :: file_name

      call build_file_name( file_name_epsilon, iq, file_name )
      call read_from_file( file_name, epsilon, file_format )
      if( is_Gamma_point ) then
          call build_file_name( file_name_epsilon_head, file_name )
          call read_from_file( file_name, epsh, file_format )
          call build_file_name( file_name_epsilon_wings1, file_name )
          call read_from_file( file_name, epsw1, file_format )
          call build_file_name( file_name_epsilon_wings2, file_name )
          call read_from_file( file_name, epsw2, file_format )
      end if
      
  end subroutine

end module
