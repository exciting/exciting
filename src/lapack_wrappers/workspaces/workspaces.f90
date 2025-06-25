!> This module provides a common workspace for LAPACK subroutines
module lapack_workspaces

    use precision, only: i32, dp
    use modmpi, only: terminate_if_false

    private
    public  :: lapack_workspace_complex_dp_t

    !> Workspace base class, note this is an invalid object
    !> that cannot be directly used
    type, abstract :: lapack_workspace_base_t
        integer(i32), private :: lrwork = -1
        integer(i32), private :: liwork = -1
        integer(i32), private :: lwork  = -1
    contains
        procedure, public   :: initialize, computed
        !> The actual implementation is left to the derived types
        procedure(allocate_workspace), deferred :: allocate_workspace
        procedure(deallocate_workspace), deferred :: deallocate_workspace
        procedure(reset), deferred :: reset
    end type lapack_workspace_base_t

    !> Create the interfaces for the deferred types
    interface
        subroutine allocate_workspace(this)
            import :: lapack_workspace_base_t
            class(lapack_workspace_base_t), intent(inout) :: this
        end subroutine allocate_workspace

        subroutine deallocate_workspace(this)
            import :: lapack_workspace_base_t
            class(lapack_workspace_base_t), intent(inout) :: this
        end subroutine deallocate_workspace

        subroutine reset(this, from_query)
            import :: lapack_workspace_base_t
            class(lapack_workspace_base_t), intent(inout) :: this
            logical, optional, intent(in) :: from_query
        end subroutine reset
    end interface

    !> Define the derived types
    type, extends(lapack_workspace_base_t) :: lapack_workspace_complex_dp_t
        real(dp), allocatable     :: rwork(:)
        integer(i32), allocatable :: iwork(:)
        complex(dp), allocatable  :: work(:)
    contains
        procedure :: reset => reset_complex_dp
        procedure :: allocate_workspace => allocate_workspace_complex_dp
        procedure :: deallocate_workspace => deallocate_workspace_complex_dp
    end type lapack_workspace_complex_dp_t

contains

    !> Init the workspace
    subroutine initialize(this, lrwork, liwork, lwork)

        !> lapack_workspace_t class to reset
        class(lapack_workspace_base_t), intent(out) :: this
        !> LWORK to init. If not present this is set 0, so computed works
        integer(i32), intent(in), optional :: lrwork
        !> LIWORK to init. If not present this is set 0, so computed works
        integer(i32), intent(in), optional :: liwork
        !> LWORK to init. If not present this is set 0, so computed works
        integer(i32), intent(in), optional :: lwork

        if (present(lrwork)) then
            this%lrwork = lrwork
        else 
            this%lrwork = 0
        end if

        if (present(liwork)) then
            this%liwork = liwork
        else 
            this%liwork = 0
        end if

        if (present(lwork)) then
            this%lwork = lwork
        else 
            this%lwork = 0
        end if

    end subroutine initialize

    !> Return true if the workspace has been already computed
    pure function computed(this) result(answer)
        
        !> lapack_workspace_t class to check if inited
        class(lapack_workspace_base_t), intent(in) :: this
        
        logical :: answer

        answer = .not. ((this%lrwork == -1) .and. (this%liwork == -1) .and. (this%lwork == -1))
        
    end function computed

    !> allocates all workspace
    subroutine allocate_workspace_complex_dp(this)
        class(lapack_workspace_complex_dp_t), intent(inout) :: this
        call terminate_if_false(this%computed(), "allocate_workspace: trying to allocate the workspace before setting its size.")
        call terminate_if_false(.not. allocated(this%work), "allocate_workspace: work array was already allocated.")
        call terminate_if_false(.not. allocated(this%iwork), "allocate_workspace: iwork array was already allocated.")
        call terminate_if_false(.not. allocated(this%rwork), "allocate_workspace: rwork array was already allocated.")

        allocate(this%work(this%lwork), this%iwork(this%liwork), this%rwork(this%lrwork))

    end subroutine allocate_workspace_complex_dp

    subroutine deallocate_workspace_complex_dp(this)
        class(lapack_workspace_complex_dp_t), intent(inout) :: this
        if(allocated(this%work)) deallocate(this%work)
        if(allocated(this%iwork)) deallocate(this%iwork)
        if(allocated(this%rwork)) deallocate(this%rwork)
    end subroutine deallocate_workspace_complex_dp

    !> Reset a workspace
    subroutine reset_complex_dp(this, from_query)

        !> lapack_workspace_t class to reset
        class(lapack_workspace_complex_dp_t), intent(inout) :: this
        !> If we are reseting the workspace from a query
        logical, optional, intent(in) :: from_query

        logical :: local_from_query

        if (present(from_query)) then
            local_from_query = from_query
        else
            local_from_query = .false.
        end if

        if (.not. local_from_query) then
            this%lrwork = -1
            this%liwork = -1
            this%lwork  = -1
            call this%deallocate_workspace()
        else
            this%lrwork = merge(0_i32, int(this%rwork(1), kind=i32), size(this%rwork) == 0)
            this%liwork = merge(0_i32, int(this%iwork(1), kind=i32), size(this%iwork) == 0)
            this%lwork  = merge(0_i32, int(this%work(1),  kind=i32), size(this%work)  == 0)
            call this%deallocate_workspace()
            call this%allocate_workspace()
        end if

    end subroutine reset_complex_dp

end module lapack_workspaces
