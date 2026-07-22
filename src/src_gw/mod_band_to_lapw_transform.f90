!> Module containing procedures to transform an operator from a
!> band representation to LAPW+LO basis, and to load and dump
!> the quantities required for such a transformation
module mod_band_to_lapw_transform
    use constants,                only: zzero, zone
    use precision,                only: dp, i32, str_512
    use modgw,                    only: Gkset, kset
    use mod_eigenvalue_occupancy, only: nstfv
    use modinput,                 only: input
    use mod_eigensystem,          only: nmatmax, nmat

    implicit none
    private

    public :: transform_from_band_representation_to_lapwlo, read_overlap_from_a_file, write_overlap_to_a_file, overlap_rootname

    !> Rootname for the overlap files
    character(len=*), parameter :: overlap_rootname = "SOVERLAP_K"

contains

    !------------------------------------------------------------------
    !> Transform an operator from band representation to LAPW+LO basis
    !> for a single k-point
    !>
    !> M^{LAPW} = adjoint(adjoint(C) S) M^{band} (adjoint(C) S)
    !------------------------------------------------------------------
    subroutine transform_from_band_representation_to_lapwlo(ik, ikir, use_full_bz, matrix_band, matrix_lapwlo, file_format)

        !> k-point index (full BZ)
        integer(i32), intent(in) :: ik
        !> k-point index (irreducible)
        integer(i32), intent(in) :: ikir
        !> Should we use ik or ikir
        logical, intent(in) :: use_full_bz
        !> Matrix in band representation
        complex(dp), intent(in)  :: matrix_band(nstfv, nstfv)
        !> Matrix in LAPW+LO representation
        complex(dp), allocatable, intent(inout) :: matrix_lapwlo(:,:)
        !> Format of the overlap files
        character(len=*), intent(in) :: file_format

        integer(i32) :: nmatp, i, j, ik_internal

        complex(dp), allocatable :: evec(:,:)
        complex(dp), allocatable :: temp(:,:)
        complex(dp), allocatable :: soverlap(:,:)

        if (allocated(matrix_lapwlo)) deallocate(matrix_lapwlo)

        ik_internal = merge(ik, ikir, use_full_bz)

        !--------------------------------------------------------------
        ! Compute LAPW(+LO) basis dimension
        !--------------------------------------------------------------
        nmatp = nmat(1,ik_internal)

        !--------------------------------------------------------------
        ! Overlap matrix S_{GG'}
        !--------------------------------------------------------------
        call read_overlap_from_a_file(soverlap, ik_internal, file_format) 

        !--------------------------------------------------------------
        ! temp = M^{band} * \adjoint(C) S
        !--------------------------------------------------------------
        allocate(temp(nstfv, nmatp))
        call zgemm('n', 'n', nstfv, nmatp, nstfv, &
                   zone, matrix_band, nstfv, &
                   soverlap, nstfv, &
                   zzero, temp, nstfv)

        !--------------------------------------------------------------
        ! M^{LAPW} = adjoint(\adjoint(C) S) * temp1
        !--------------------------------------------------------------
        allocate(matrix_lapwlo(nmatp, nmatp), source=zzero)
        call zgemm('c', 'n', nmatp, nmatp, nstfv, &
                   zone, soverlap, nstfv, &
                   temp, nstfv, &
                   zzero, matrix_lapwlo, nmatp)

        !--------------------------------------------------------------
        ! Cleanup
        !--------------------------------------------------------------
        deallocate(temp, soverlap)
        
    end subroutine transform_from_band_representation_to_lapwlo

    subroutine write_overlap_to_a_file(overlap, ik, file_format)
        use gw_io, only: build_file_name, write_to_file
        !> Overlap
        complex(dp),  intent(in) :: overlap(:,:)
        !> Index of the current k-point
        integer(i32), intent(in) :: ik
        !> Format of the output file
        character(len=*), optional, intent(in) :: file_format
        
        character(len=str_512) :: file_name
        character(len=str_512) :: file_format_local

        if (present(file_format)) then
            file_format_local = file_format
        else
            file_format_local = "binary"
        end if

        call build_file_name( overlap_rootname, ik, file_name )
        call write_to_file( file_name, overlap(:,:), file_format_local)

    end subroutine write_overlap_to_a_file

    subroutine read_overlap_from_a_file(overlap, ik, file_format)
        use gw_io, only: build_file_name, read_from_file
        !> Overlap
        complex(dp),  allocatable, intent(inout) :: overlap(:,:)
        !> Index of the current k-point
        integer(i32), intent(in) :: ik
        !> Format of the output file
        character(len=*), optional, intent(in) :: file_format
        
        character(len=str_512) :: file_name
        character(len=str_512) :: file_format_local

        if (present(file_format)) then
            file_format_local = file_format
        else
            file_format_local = "binary"
        end if

        call build_file_name( overlap_rootname, ik, file_name )
        call read_from_file( file_name, overlap, file_format_local)

    end subroutine read_overlap_from_a_file
      

end module mod_band_to_lapw_transform

