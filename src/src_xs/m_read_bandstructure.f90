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
      use modmpi, only: terminate
      use m_getunit
      use mod_hdf5, only: hdf5_read

      character(*), intent(in) :: fname
      
      logical :: fex
      character(2) :: commentchar
      integer(4) :: un 
      integer(4) :: ib1, ib2, nsteps, istep, ib, ibread, istepread
#ifdef _HDF5_
      ! bandstructure path and corresponding DFT energy values
      real(8), allocatable :: vkl_(:,:), kpath2_(:), energy_(:,:)
      character(256) :: fhdf5_, gname_, params_
#endif
      call clear_bandstructure()
      
      ! Check if bandstructure.dat is in working directory
      inquire(file=trim(fname), exist=fex)
      if(.not. fex) then 
        write(*, '("Error(read_bandstructure): Bandstructure file not found.")')
        call terminate
      end if

      ! Read in Bandstructure
# ifndef _HDF5_
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
#else
      ! define HDF5 file
      fhdf5_="property.h5"
      ! band ranges are hard-coded to be the whole range of Kohn-Sham states
      allocate(brange_(3))
      ib1=1
      params_='/bandstructure/parameters/'
      gname_='/bandstructure/'
      call hdf5_read(fhdf5_, params_, "nkpt", nsteps)
      call hdf5_read(fhdf5_, params_, "nstsv", ib2)
      if (allocated(vkl_)) deallocate(vkl_)
      allocate(vkl_(3,nsteps))
      call hdf5_read(fhdf5_, params_, "vkl", vkl_(1,1), shape(vkl_))
      brange_(1)=ib1
      brange_(2)=ib2
      brange_(3)=nsteps
      ! get actual data
      if (allocated(kpath2_)) deallocate(kpath2_)
      allocate(kpath2_(nsteps))
      call hdf5_read(fhdf5_, gname_, "points",kpath2_(1), shape(kpath2_))
      if (allocated(energy_)) deallocate(energy_)
      allocate(energy_(ib2,nsteps))
      call hdf5_read(fhdf5_, gname_, "evalsv", energy_(1,1), shape(energy_))
      
      ! allocate output
      allocate(vklpath_(3, nsteps, ib1:ib2))
      allocate(kpathlength_(nsteps, ib1:ib2))
      allocate(energyval_(nsteps, ib1:ib2))

      ! get data in the right format
      ! the output data contains some redundancy, i.e. the kpath and kpathlength
      !is stored for each band in memory, but only once in file.
      do istep=1,nsteps
        do ib =ib1, ib2
          vklpath_(:,istep,ib)=vkl_(:,istep)
          kpathlength_(istep,ib)=kpath2_(istep)
          energyval_(istep,ib)=energy_(ib,istep)
        end do
      end do
      
      ! deallocate intermediate arrays
      deallocate(energy_, kpath2_, vkl_)
#endif
    end subroutine read_bandstructure

end module m_read_bandstructure
