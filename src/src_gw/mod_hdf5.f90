
!===============================================================================
!
! Author: Anton Kozhevnikov
! Added:  10.06.2014 by DIN
! 
! Modifications:
!   Added character type, 11.06.2014 by DIN
!
!===============================================================================

module mod_hdf5

    interface hdf5_write
        module procedure hdf5_write_i4, &
        &                hdf5_write_d,  &
        &                hdf5_write_z,  &
        &                hdf5_write_c,  &
        &                hdf5_write_l
    end interface

    interface hdf5_read
        module procedure hdf5_read_i4, &
        &                hdf5_read_d,  &
        &                hdf5_read_z,  &
        &                hdf5_read_c,  &
        &                hdf5_read_l
    end interface

    public hdf5_initialize
    public hdf5_finalize
    public hdf5_create_file
    public hdf5_create_group

contains

!-------------------------------------------------------------------------------
    subroutine hdf5_initialize
#ifdef _HDF5_
        use hdf5
        implicit none
        integer ierr
        call h5open_f(ierr)
#endif
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_finalize
#ifdef _HDF5_
        use hdf5
        implicit none
        integer ierr
        call h5close_f(ierr)
#endif
    end subroutine

!-------------------------------------------------------------------------------
    subroutine hdf5_create_file(fname)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
#ifdef _HDF5_
        integer ierr
        integer(hid_t) h5_root_id
        call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_file): h5fcreate_f returned ",I6)')ierr
          goto 10
        endif
        call h5fclose_f(h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_file): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  fname: ",A)')trim(fname)
        stop
#endif
    end subroutine
    
!-------------------------------------------------------------------------------
    logical function hdf5_exist_group(fname,path,gname)
        ! Check if group with the given name exists
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: gname
#ifdef _HDF5_
        integer(hid_t):: h5_root_id, h5_group_id
        integer ierr
 
        call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5gopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5lexists_f(h5_group_id,trim(gname),hdf5_exist_group,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5lexists_f returned ",I6)')ierr
          goto 10
        endif
        call h5fclose_f(h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  fname : ",A)')trim(fname)
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  gname : ",A)')trim(gname)  
        stop
#endif
    end function hdf5_exist_group

!-------------------------------------------------------------------------------    
    subroutine hdf5_create_group(fname,path,gname)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: gname
#ifdef _HDF5_
        integer(hid_t) h5_root_id,h5_group_id,h5_new_group_id
        integer ierr

        call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5gopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5gcreate_f(h5_group_id,trim(gname),h5_new_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5gcreate_f returned ",I6)')ierr
          goto 10
        endif
        call h5gclose_f(h5_new_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5gclose_f for the new group returned ",I6)')ierr
          goto 10
        endif
        call h5gclose_f(h5_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5gclose_f for the existing path returned ",I6)')ierr
          goto 10
        endif
        call h5fclose_f(h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  fname : ",A)')trim(fname)
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  gname : ",A)')trim(gname)  
        stop
#endif
    end subroutine

!-------------------------------------------------------------------------------
    subroutine hdf5_write_i4(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        integer(4), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_write_array_i4(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_write_d(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        real(8), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_write_z(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        complex(8), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)+1
          allocate(dims_(ndims))
          dims_(1)=2
          dims_(2:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=2
        endif
        call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine


!-------------------------------------------------------------------------------    
    subroutine hdf5_read_i4(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        integer(4), intent(out) :: val
        integer(4), optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_read_array_i4(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_read_d(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        real(8), intent(out) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_read_z(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        complex(8), intent(out) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)+1
          allocate(dims_(ndims))
          dims_(1)=2
          dims_(2:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=2
        endif
        call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine hdf5_write_c(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        character*(*), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_write_array_c(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine

!-------------------------------------------------------------------------------
    subroutine hdf5_write_l(fname,path,dname,val)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        logical, intent(in) :: val
#ifdef _HDF5_
        if (val) then
          call hdf5_write_array_c("true",1,1,fname,path,dname)
        else
          call hdf5_write_array_c("false",1,1,fname,path,dname)
        end if
#endif
    end subroutine
!-------------------------------------------------------------------------------
    subroutine hdf5_read_c(fname,path,dname,val,dims)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        character*(*), intent(out) :: val
        integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_read_array_c(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
#endif
    end subroutine    
!-------------------------------------------------------------------------------
    subroutine hdf5_read_l(fname,path,dname,val)
#ifdef _HDF5_
        use hdf5
#endif
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        logical, intent(out) :: val
#ifdef _HDF5_
        character(100) :: val_
        call hdf5_read_array_c(val_,1,1,fname,path,dname)
        if (trim(val_) == "true") then
          val=.TRUE.
        else
          val=.FALSE.
        end if
#endif
    end subroutine    
end module

!===============================================================================

#ifdef _HDF5_

!-------------------------------------------------------------------------------
subroutine hdf5_write_array_i4(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    integer(4), intent(in) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
    integer ierr,i
    integer(hsize_t), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_INTEGER,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5dcreate_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5sclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims  : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  nm    : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_write_array_d(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    real(8), intent(in) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
    integer ierr,i
    integer(hsize_t), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5dcreate_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5sclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims  : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  nm    : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_i4(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    integer(4), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataset_id,group_id
    integer ierr,i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo

    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_d(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    real(8), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataset_id,group_id
    integer ierr,i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_write_array_c(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    character(*), intent(in) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm
    integer(HID_T) :: h5_root_id,dataspace_id,dataset_id,group_id
    integer(HID_T) :: filetype
    integer(SIZE_T) :: sdim
    integer :: ierr, i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    integer(HSIZE_T), dimension(ndims+1) :: data_dims
    integer(SIZE_T), dimension(ndims) :: s_len
    character*100 :: errmsg

    do i = 1, ndims
      h_dims(i) = dims(i)
    enddo
    do i = 1, product(h_dims)
      s_len(i) = len_trim(a(i))
    end do
    data_dims(1) = len(a(1))
    do i = 1, ndims
      data_dims(i+1) = h_dims(i)
    enddo
    
    call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    !------------------------------------------------------
    ! Create variable size datatype
    !call H5Tcopy_f(H5T_FORTRAN_S1,filetype,ierr)
    !sdim = maxval(s_len)
    !call H5Tset_size_f(filetype,sdim,ierr)
    CALL H5Tcopy_f(H5T_STRING,filetype,ierr)
    CALL H5Tset_strpad_f(filetype,H5T_STR_NULLPAD_F,ierr)
    !------------------------------------------------------
    call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dcreate_f(group_id,trim(nm),filetype,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5dcreate_f returned ",I6)')ierr
      goto 10
    endif
    !------------------------------------------------------
    !call h5dwrite_f(dataset_id,filetype,a,h_dims,ierr)
    CALL h5dwrite_vl_f(dataset_id,filetype,a,data_dims,s_len,ierr,dataspace_id)
    !------------------------------------------------------
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5sclose_f returned ",I6)')ierr
      goto 10
    endif
    call H5Tclose_f(filetype,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5tclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims  : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  nm    : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_c(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    character(*), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) :: h5_root_id,dataset_id,dataspace_id,group_id
    integer(HID_T) :: filetype
    integer :: ierr, i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    integer(HSIZE_T), dimension(ndims+1) :: data_dims
    integer(SIZE_T), dimension(ndims) :: s_len
    integer(HSIZE_T), dimension(1:2) :: maxdims
    character*100 :: errmsg

    do i = 1, ndims
      h_dims(i) = dims(i)
    enddo
    do i = 1, product(h_dims)
      s_len(i) = len_trim(a(i))
    end do
    data_dims(1) = len(a(1))
    do i = 1, ndims
      data_dims(i+1) = h_dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    !------------------------------------------------------
    ! Get the datatype
    CALL H5Dget_type_f(dataset_id,filetype,ierr)
    ! Get dataspace
    CALL H5Dget_space_f(dataset_id,dataspace_id,ierr)
    CALL H5Sget_simple_extent_dims_f(dataspace_id,h_dims,maxdims,ierr)
    !call h5dread_f(dataset_id,H5T_FORTRAN_S1,a,h_dims,ierr)
    CALL h5dread_vl_f(dataset_id,filetype,a,data_dims,s_len,ierr,dataspace_id)
    !------------------------------------------------------    
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    CALL H5Tclose_f(filetype,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5tclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

#endif

