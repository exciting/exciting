! This "template" implements function enforces the symmetry on the input vector
! So code repetition is avoided
! TEMPLATE
! TYPE1 : type of the input_vector
! PRECISION1 : precision of TYPE1
#if !defined(TYPE1) || !defined(PRECISION1)
#error "TYPE1 and PRECISION1 must be defined before using this template (enfore_degeneracy_template.inc)"
#endif

#define PROCEDURNAME1 enforce_degeneracy
#define PASTE(X) X
#define PASTE2(X) PASTE(X)_
#define CONCATENATE(X,Y) PASTE2(X)Y
#define EXTENSION CONCATENATE(TYPE1,PRECISION1)
#define GENERATENAME(name) CONCATENATE(name,EXTENSION) 

!> That function enforces the symmetry on the input vector
subroutine GENERATENAME(PROCEDURNAME1)(ikp, input_vector)
    !> irreducible k-point index
    integer(i32), intent(in) :: ikp
    !> Data at ikp in which the degeneracy must be strictly enforced (TYPE1 of PRECISION1)
    TYPE1(PRECISION1), intent(inout)  :: input_vector(:)

    ! Index to iterate over degenerate spaces
    integer(i32) :: ispace
    ! Upper and lower limits
    integer(i32) :: ispace_init, ispace_final

    ! Fail if degenerate_subspaces is not allocated
    if (.not. allocated(degenerate_subspaces)) then
        call terminate("Error(GENERATENAME(PROCEDURNAME1)): degenerate_subspaces is not allocated. " // &
        "Please call initialize_degeneracy_module before.")
    end if

    do ispace = 1_i32, size(degenerate_subspaces,2)
        if (degenerate_subspaces(3,ispace,ikp) == -1) cycle
        input_vector(degenerate_subspaces(1,ispace,ikp):degenerate_subspaces(2,ispace,ikp)) =  &
            sum(input_vector(degenerate_subspaces(1,ispace,ikp):degenerate_subspaces(2,ispace,ikp))) / degenerate_subspaces(3,ispace,ikp)
    end do

end subroutine GENERATENAME(PROCEDURNAME1)

! Clean functions
#undef PASTE
#undef PASTE2 
#undef CONCATENATE
#undef EXTENSION
#undef GENERATENAME

! Clean the name of the procedure
#undef PROCEDURNAME1

! Clean templated type and precision
#undef TYPE1 
#undef PRECISION1
