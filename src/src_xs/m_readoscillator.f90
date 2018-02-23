module m_readoscillator

  implicit none

  contains

    subroutine readoscillator(iqmt, io1, evals, bindevals, oscir)
      use modinput
      use m_getunit
      use m_genfilname

      integer(4), intent(in) :: iqmt
      integer(4), intent(in) :: io1
      real(8), allocatable, intent(inout), optional :: evals(:), bindevals(:)
      complex(8), allocatable, intent(inout), optional :: oscir(:)

      logical :: fex, fcoup
      integer(4) :: un, i, nexc, idummy, ncommentlines, nlines
      real(8), allocatable :: reevals(:), rebindevals(:)
      real(8), allocatable :: absoscir(:), reoscir(:), imoscir(:)
      character(256) :: frmt, tdastring, bsetypestring, scrtypestring
      character(256) :: excitondir, fnexc

      character(*), parameter :: thisname = "readoscillator"

      if(iqmt /= 1 .and. io1/=1 ) then 
        write(*,*) "iqmt /= 1 but io1 /= 1"
        stop
      end if

      excitondir='EXCITON'

      fcoup = input%xs%bse%coupling
      if(fcoup) then
        tdastring=''
      else
        if(input%xs%bse%chibarq) then 
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if
      if(input%xs%bse%bsetype == "IP") then
        tdastring=''
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      frmt='(I8,5(1x,E23.16))'


      if(iqmt /= 1) then 
        call genfilname(dirname=trim(excitondir), basename='EXCITON', iqmt=iqmt,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
      else
        call genfilname(dirname=trim(excitondir), basename='EXCITON', tq0=.true., oc1=io1, oc2=io1,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
      end if

      inquire(file=trim(fnexc), exist=fex)
      if(.not. fex) then 
        write(*,'("Error(",a,"): File ",a," does not exist")')&
         & trim(thisname), trim(fnexc)
        stop
      end if
      write(*,'("Info(",a,"): Opening file ",a)') trim(thisname), trim(fnexc)

      ! Get number of excitons
      ncommentlines=14
      call getlines(trim(fnexc), nlines)
      nexc = nlines - ncommentlines

      if(allocated(reevals)) deallocate(reevals)
      allocate(reevals(nexc))
      if(allocated(rebindevals)) deallocate(rebindevals)
      allocate(rebindevals(nexc))
      if(allocated(absoscir)) deallocate(absoscir)
      allocate(absoscir(nexc))
      if(allocated(reoscir)) deallocate(reoscir)
      allocate(reoscir(nexc))
      if(allocated(imoscir)) deallocate(imoscir)
      allocate(imoscir(nexc))

      call getunit(un) 

      open(un, file=trim(fnexc), action="read", form="formatted", status="old")

      do i=1, ncommentlines
        read(un,*)
      end do
      do i=1, nexc
        read(un, fmt=frmt) idummy, reevals(i), rebindevals(i),&
          & absoscir(i), reoscir(i), imoscir(i)
      end do

      close(un)

      if(present(evals)) then 
        if(allocated(evals)) deallocate(evals)
        allocate(evals(nexc))
        evals = reevals
      end if

      if(present(bindevals)) then 
        if(allocated(bindevals)) deallocate(bindevals)
        allocate(bindevals(nexc))
        bindevals = rebindevals
      end if

      if(present(oscir)) then 
        if(allocated(oscir)) deallocate(oscir)
        allocate(oscir(nexc))
        oscir = cmplx(reoscir,imoscir,8)
      end if

      contains

        subroutine getlines(fname, nlines)

          character(*), intent(in) :: fname
          integer(4), intent(out) :: nlines

          integer(4) :: un, stat

          call getunit(un) 

          open(un, file=trim(fname), action="read", form="formatted", status="old")

          nlines = 0
          do
            read(un,*,iostat=stat)
            if(stat < 0) then 
              exit
            else if(stat > 0) then 
              stop 'io error'
            end if
            nlines = nlines + 1
          end do

          close(un)

        end subroutine getlines

    end subroutine readoscillator

end module m_readoscillator
