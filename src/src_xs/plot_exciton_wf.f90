
subroutine plot_exciton_wf(lambda,rh)

  use modinput
  use modmain
  use modxs
  use mod_rgrid
  use mod_xsf_format
  use modmpi, only : rank
  implicit none

  ! input/output
  integer, intent(in) :: lambda
  real(8), intent(in) :: rh(3)

  ! local
  logical :: exist
  integer :: iostat
  integer :: nsta1, nrnst1, nrnst3, hamsiz, ievec
  character(256) :: locext

  real(8),    allocatable :: beval(:)
  complex(8), allocatable :: bevec(:,:)

  if (.not.associated(input%properties%wfplot)) then
    write(*,*)
    write(*,*) 'WARNING(plot_exciton_wf): wfplot parameters are not specified!'
    write(*,*) '                          Use defaults.'
    input%properties%wfplot => getstructwfplot(emptynode)
  end if

  inquire(file='EXCCOEFF.bin', exist=exist)
  if ( (.not. exist) .and. (rank==0) ) then
    write(*,*)
    write(*,'("Error(bse): file EXCCOEFF.bin does not exist!")')
    write(*,*)
    stop
  else
    open(50,File='EXCCOEFF.bin', & 
            Action='READ',Form='UNFORMATTED', IOstat=iostat)
    if ( (iostat.ne.0) .and. (rank==0) ) then
      write(*,*) iostat
      write(*,'("Error(bse): error reading EXCCOEFF.bin")')
      write(*,*)
      stop
    end if
    read(50) input%xs%storeexcitons%MinNumberExcitons, input%xs%storeexcitons%MaxNumberExcitons, &
        & nkptnr, nsta1, nrnst1, nrnst3, hamsiz
    beval=0.0
    bevec=0.0
    do ievec = input%xs%storeexcitons%MinNumberExcitons, input%xs%storeexcitons%MaxNumberExcitons
       read(50) beval(ievec), bevec(1:hamsiz,ievec)
    end do
    close(50)
  end if
         
  locext = filext 
  filext = 'EVECFV_QMT000.OUT'
  
  inquire(file=trim(filext), exist=exist)
  if ( (.not. exist) .and. (rank==0) ) then
    write(*,*)
    write(*,'("Error(bse): file EVECFV_QMT000.OUT does not exist!")')
    write(*,*)
    stop
  else
    open(50,File=trim(filext), & 
         Action='READ',Form='UNFORMATTED', IOstat=iostat)
    if ( (iostat.ne.0) .and. (rank==0) ) then
      write(*,*) iostat
      write(*,'("Error(bse): error reading EVECFV_QMT000.OUT")')
      write(*,*)
      stop
    end if
  end if
   
  close(50)
      
  filext = locext

  return
end subroutine