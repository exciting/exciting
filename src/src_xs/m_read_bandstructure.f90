module m_read_bandstructure

  implicit none

  public

  ! Banstructure arrays
  integer(4), allocatable :: brange_(:)
  real(8), allocatable :: vklpath_(:,:,:), kpathlength_(:,:), energyval_(:,:)

  contains 

    subroutine clear_bandstructure()
      if(allocated(brange_)) deallocate(brange_)
      if(allocated(vklpath_)) deallocate(vklpath_)
      if(allocated(kpathlength_)) deallocate(kpathlength_)
      if(allocated(energyval_)) deallocate(energyval_)
    end subroutine clear_bandstructure

    subroutine read_bandstructure(fname)
      use modmpi
      use m_getunit

      character(*), intent(in) :: fname
      
      logical :: fex
      character(2) :: commentchar
      integer(4) :: un 
      integer(4) :: ib1, ib2, nsteps, istep, ib, ibread, istepread

      call clear_bandstructure()

      ! Check if bandstructure.dat is in working directory
      inquire(file=trim(fname), exist=fex)
      if(.not. fex) then 
        write(*, '("Error(read_bandstructure): Bandstructure file not found.")')
        call terminate
      end if

      ! Read in Bandstructure

      call getunit(un)
      open(unit=un, file=trim(fname), form='formatted', action='read')

      !   Read in index of lowest and highest band and the number of steps along the path
      allocate(brange_(3))
      read(un,*) commentchar, brange_(1), brange_(2), brange_(3)
      ib1=brange_(1)
      ib2=brange_(2)
      nsteps=brange_(3)

      allocate(vklpath_(3, nsteps, ib1:ib2))
      allocate(kpathlength_(nsteps, ib1:ib2))
      allocate(energyval_(nsteps, ib1:ib2))

      ! Read in kpoints, pathlengths and energy values
      do ib = ib1, ib2
        do istep = 1, nsteps
          read(un,*)&
            & ibread, istepread, vklpath_(:,istep,ib), kpathlength_(istep, ib),&
            & energyval_(istep, ib)
        end do
        read(un,*)
      end do
      close(un)

      write(*, '("Info(read_bandstructure):&
        & Bandstructure read for band interval:", 2i8)') ib1, ib2

    end subroutine read_bandstructure

end module m_read_bandstructure
