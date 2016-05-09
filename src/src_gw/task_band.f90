
subroutine task_band()

  use modinput
  use modmain
  use modgw
  use modmpi

  implicit none

  integer(4)    :: ik, ib
  real(8)       :: tstart, tend
  integer       :: i, j
  character(80) :: s, fname
  logical       :: exist

  !------------------------
  ! Read KS bandstructure
  !------------------------

  call init0()
  ! call init1

  call readfermi()

  fname = 'bandstructure.dat'
  inquire(File=fname, Exist=exist)
  if (.not.exist) then
    write(*,*) 'ERROR(task_band): bandstructure.dat file is not found!'
    write(*,*) '    Run properties/bandstructure first to produce KS spectrum.'
    stop
  end if

  open(70, File='bandstructure.dat', Action='Read', Status='Old')
  read(70,*) s, nstsv, nkpt
  write(*,*) trim(s), nstsv, nkpt
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(dpp1d)) deallocate(dpp1d)
  allocate(dpp1d(nkpt))
  if (allocated(evalsv)) deallocate(evalsv)
  allocate(evalsv(nstsv,nkpt))
  do ib = 1, nstsv
    do ik = 1, nkpt
      ! Note: evalsv are already shifted to E_f = 0
      read(70,*) i, j, vkl(:,ik), dpp1d(ik), evalsv(ib,ik)
    end do
    read(70,*) ! skip line
  end do
  close(70)

  !--------------------------------------------------------------
  ! read QP energies from file and perform Fourier interpolation 
  !--------------------------------------------------------------
  call getevalqp(nkpt,vkl,evalsv)

  !----------------------------------
  ! write QP bandstructure to disk
  !----------------------------------
  open(50,file='BAND-QP.OUT',action='WRITE',form='FORMATTED')
  do ib = ibgw, min(nbgw,nstsv)
    do ik = 1, nkpt
      write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)+efermi-eferqp
    end do !ik
    write(50,*)
  end do !ib
  close(50)

  return
end subroutine
