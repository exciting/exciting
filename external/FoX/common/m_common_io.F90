module m_common_io

#ifndef DUMMYLIB
  use m_common_error, only : FoX_error

  implicit none
  private

  ! Basic  I/O tools

  integer, save :: io_eor
  integer, save :: io_eof
  integer, save :: io_err
  logical, save :: io_init=.false.

  interface  setup_io
      module procedure setup_io_scratch, setup_io_with_data
  end interface setup_io

  public :: io_eor
  public :: io_eof
  public :: io_err
  public :: get_unit
  public :: setup_io

contains

  subroutine setup_io_scratch()
    if (io_init) return 
    call find_eor_eof(io_eor, io_eof)
    io_init = .true. 
  end subroutine setup_io_scratch

  subroutine setup_io_with_data(err_code, eor_code, eof_code)
     implicit none
     integer,intent(in)  :: err_code, eor_code, eof_code
     if (io_init) return 
     io_err = err_code
     io_eor  = eor_code
     io_eof  = eof_code
     io_init = .true.
   end subroutine setup_io_with_data 

  


  subroutine get_unit(lun,iostat)
    ! Get an available Fortran unit number
    integer, intent(out)  :: lun
    integer, intent(out)  :: iostat

    integer :: i
    logical :: unit_used

    do i = 10, 99
      lun = i
      inquire(unit=lun,opened=unit_used)
      if (.not. unit_used) then
        iostat = 0
        return
      endif
    enddo
    iostat = -1
    lun = -1
  end subroutine get_unit


  subroutine find_eor_eof(io_eor,io_eof)
    ! Determines the values of the iostat values for End of File and 
    ! End of Record (in non-advancing I/O)
#ifdef __NAG__
    use f90_iostat
#endif
    integer, intent(out)           :: io_eor
    integer, intent(out)           :: io_eof

#ifdef __NAG__
    io_eor = ioerr_eor
    io_eof = ioerr_eof
#else
    integer           :: lun, iostat
    character(len=1)  :: c

    call get_unit(lun,iostat)

    if (iostat /= 0) call FoX_error("Out of unit numbers")

    open(unit=lun,status="scratch",form="formatted", &
      action="readwrite",position="rewind",iostat=iostat)
    if (iostat /= 0) call FoX_error("Cannot open test file")

    write(unit=lun,fmt=*)  "a"
    write(unit=lun,fmt=*)  "b"

    rewind(unit=lun)

    io_eor = 0
    do
      read(unit=lun,fmt="(a1)",advance="no",iostat=io_eor) c
      if (io_eor /= 0) exit
    enddo

    io_eof = 0
    do
      read(unit=lun,fmt=*,iostat=io_eof)
      if (io_eof /= 0) exit
    enddo

    close(unit=lun,status="delete")
#endif
    
    ! Invent an io_err ...
    io_err = 1
    do
      if (io_err/=0.and.io_err/=io_eor.and.io_err/=io_eof) exit
      io_err = io_err + 1
    end do
  end subroutine find_eor_eof

#endif
end module m_common_io
