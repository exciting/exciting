!----------------------------!
!   Dielectric function      !
!----------------------------!

module mod_dielectric_function
    use gw_io, only: build_file_name, write_to_file, read_from_file
    use precision, only: i32, dp
#include "offload.fpp"

    implicit none

    private

    ! dielectric function \epsilon(q)
    complex(dp), public, allocatable :: epsilon(:,:,:)

    !-------------------------------------------------
    ! Analytical treatment of q=0 singularity
    !-------------------------------------------------

    ! valence-valence momentum matrix elements
    complex(dp), public, allocatable :: pmatvv(:,:,:)

    ! core-valence momentum matrix elements
    complex(dp), public, allocatable :: pmatcv(:,:,:)

    ! head of the dielectric function (tensor)
    complex(dp), public, allocatable :: epsh(:,:,:)

    ! the vertical wing of the dielectric matrix (vector)
    complex(dp), public, allocatable :: epsw1(:,:,:)

    ! the horizontal wing of the dielectric matrix (vector)
    complex(dp), public, allocatable :: epsw2(:,:,:)

    !----------------------------------------------------------------------
    ! Used for calculating the macroscopic dielectric function (task_emac)
    !----------------------------------------------------------------------

    complex(dp), public, allocatable :: eps00(:,:,:)

    !----------------------------------------------------------------------
    ! files containing data on PMAT and PMATCOR
    !----------------------------------------------------------------------
    character(len=*), parameter, public :: fname_pmatvv='PMATVV.OUT'
    character(len=*), parameter, public :: fname_pmatcv='PMATCV.OUT'
    
    !----------------------------------------------------------------------
    ! files to store the dielectric function
    !----------------------------------------------------------------------
    character(len=*), parameter, private :: file_name_epsilon = 'EPSILON-GW_Q'
    character(len=*), parameter, private :: file_name_epsilon_irreducible = 'EPSILON-GW_IQ'
    character(len=*), parameter, private :: file_name_epsilon_head = 'EPSH'
    character(len=*), parameter, private  :: file_name_epsilon_wings1 = 'EPSW1'
    character(len=*), parameter, private  :: file_name_epsilon_wings2 = 'EPSW2'

    !----------------------------------------------------------------------
    ! files to store the inverse of the dielectric function
    !----------------------------------------------------------------------
    character(len=*), parameter, private :: file_name_inverse_epsilon = 'INVERSE-EPSILON_Q'
    character(len=*), parameter, private :: file_name_inverse_epsilon_irreducible = 'INVERSE-EPSILON_IQ'
    character(len=*), parameter, private :: file_name_inverse_epsilon_head = 'INVERSE-EPSH'
    character(len=*), parameter, private :: file_name_inverse_epsilon_wings1 = 'INVERSE-EPSW1'
    character(len=*), parameter, private :: file_name_inverse_epsilon_wings2 = 'INVERSE-EPSW2'
    


    integer(i32), parameter   :: max_string_length = 40

    public :: write_epsilon_to_file, init_dielectric_function, delete_dielectric_function, &
              write_inverse_epsilon_to_file, read_epsilon_from_file, read_inverse_epsilon_from_file

    
contains

    subroutine init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
        use constants, only: zzero
        implicit none
        integer, intent(in) :: mbsiz
        integer, intent(in) :: iomstart, iomend
        logical, intent(in) :: Gamma
        ! q-dependent dielectric function
        if (allocated(epsilon)) then 
          OMP_OFFLOAD target exit data map(delete: epsilon)
          deallocate(epsilon)
        end if
        allocate(epsilon(mbsiz,mbsiz,iomstart:iomend), source=zzero)
        OMP_OFFLOAD target enter data map(always, to: epsilon)
        ! head and wings of the dielectric function when q->0
        if (Gamma) then
          if (allocated(epsh)) deallocate(epsh)
          allocate(epsh(3,3,iomstart:iomend), source=zzero)
          
          if (allocated(epsw1)) deallocate(epsw1)
          allocate(epsw1(mbsiz,3,iomstart:iomend), source=zzero)
          
          if (allocated(epsw2)) deallocate(epsw2)
          allocate(epsw2(mbsiz,3,iomstart:iomend), source=zzero)
          
          ! macroscopic dielectric tensor
          if (allocated(eps00)) deallocate(eps00)
          allocate(eps00(3,3,iomstart:iomend), source=zzero)

        end if ! Gamma
    end subroutine

    subroutine delete_dielectric_function(Gamma)
        logical, intent(in) :: Gamma
        if (allocated(epsilon)) then 
          OMP_OFFLOAD target exit data map(delete: epsilon)
          deallocate(epsilon)
        end if
        if (Gamma) then
          if (allocated(epsh)) deallocate(epsh)
          if (allocated(epsw1)) deallocate(epsw1)
          if (allocated(epsw2)) deallocate(epsw2)
          if (allocated(eps00)) deallocate(eps00)
        end if
    end subroutine


    subroutine write_epsilon_to_file( iq, is_Gamma_point, file_format, irreducible )
      integer(i32), intent(in)  :: iq 
      logical, intent(in)       :: is_Gamma_point
      character(len=*), intent(in) :: file_format
      logical, intent(in), optional :: irreducible
    
      character(len=max_string_length) :: file_name

      logical :: irreducible_local 

      if (present(irreducible)) then
        irreducible_local = irreducible
      else
        irreducible_local = .false.
      end if

      if (irreducible_local) then 
        call build_file_name( file_name_epsilon_irreducible, iq, file_name )
      else
        call build_file_name( file_name_epsilon, iq, file_name )
      end if

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


  subroutine read_epsilon_from_file( iq, is_Gamma_point, file_format, irreducible )
      integer(i32), intent(in)  :: iq 
      logical, intent(in)       :: is_Gamma_point
      character(len=*), intent(in) :: file_format
      logical, intent(in), optional :: irreducible
    
      character(len=max_string_length) :: file_name

      logical :: irreducible_local 

      if (present(irreducible)) then
        irreducible_local = irreducible
      else
        irreducible_local = .false.
      end if

      if (irreducible_local) then
        call build_file_name( file_name_epsilon_irreducible, iq, file_name )
      else 
        call build_file_name( file_name_epsilon, iq, file_name )
      end if

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

  !> Write the inverse of epsilon into files
  !> Attention: actually exciting stores the inverse of epsilon in
  !> the same matrices as the dielectric matrix
  subroutine write_inverse_epsilon_to_file( iq, is_Gamma_point, file_format, irreducible )
    !> q-point index
    integer(i32), intent(in)  :: iq 
    !> When true, it is the \(\Gamma\) point
    logical, intent(in)       :: is_Gamma_point
    !> Format of output file
    character(len=*), intent(in) :: file_format
    !> Using irreducible points indexing
    logical, intent(in), optional :: irreducible
  
    integer(i32), parameter   :: max_length = 40
    character(len=max_length) :: file_name
    logical :: irreducible_local 

    if (present(irreducible)) then
      irreducible_local = irreducible
    else
      irreducible_local = .false.
    end if
  
    if (irreducible_local) then
      call build_file_name( file_name_inverse_epsilon_irreducible, iq, file_name )
    else
      call build_file_name( file_name_inverse_epsilon, iq, file_name )
    end if

    call write_to_file( file_name, epsilon, file_format )
    if( is_Gamma_point ) then
        call build_file_name( file_name_inverse_epsilon_head, file_name )
        call write_to_file( file_name, epsh, file_format )
        call build_file_name( file_name_inverse_epsilon_wings1, file_name )
        call write_to_file( file_name, epsw1, file_format )
        call build_file_name( file_name_inverse_epsilon_wings2, file_name )
        call write_to_file( file_name, epsw2, file_format )
    end if
    
  end subroutine


  !> Read the inverse of epsilon from files
  !> Attention: actually exciting stores the inverse of epsilon in
  !> the same matrices as the dielectric matrix
  subroutine read_inverse_epsilon_from_file( iq, is_Gamma_point, file_format, irreducible )
    !> q-point index
    integer(i32), intent(in)  :: iq 
    !> If true, the actual q-point refers to the \( \Gamma \) point
    logical, intent(in)       :: is_Gamma_point
    !> File format used in the files where the inverse of epsilon is stored
    character(len=*), intent(in) :: file_format
    !> Using irreducible points indexing
    logical, intent(in), optional :: irreducible
  
    integer(i32), parameter   :: max_length = 40
    character(len=max_length) :: file_name
    logical :: irreducible_local 

    if (present(irreducible)) then
      irreducible_local = irreducible
    else
      irreducible_local = .false.
    end if
  
    if (irreducible_local) then
      call build_file_name( file_name_inverse_epsilon_irreducible, iq, file_name )
    else
      call build_file_name( file_name_inverse_epsilon, iq, file_name )
    end if

    call read_from_file( file_name, epsilon, file_format )
    if( is_Gamma_point ) then
        call build_file_name( file_name_inverse_epsilon_head, file_name )
        call read_from_file( file_name, epsh, file_format )
        call build_file_name( file_name_inverse_epsilon_wings1, file_name )
        call read_from_file( file_name, epsw1, file_format )
        call build_file_name( file_name_inverse_epsilon_wings2, file_name )
        call read_from_file( file_name, epsw2, file_format )
    end if
    
  end subroutine

end module
