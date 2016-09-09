subroutine storeexcitons(hamsize, nexc, nkptnr, iuref, bcou, smap, beval, bevec)
  use modxs, only: bcbs
  use modinput, only: input
  use modmpi, only: rank
  use m_getunit

  implicit none

  ! I/O
  integer(4), intent(in) :: hamsize, nexc, nkptnr, iuref, smap(hamsize,3)
  type(bcbs), intent(in) :: bcou
  real(8), intent(in) :: beval(hamsize)
  complex(8), intent(in) :: bevec(hamsize,nexc)

  ! Local
  integer(4) :: stat, ievec, unexc

  !-----------------------------------------------------
  ! Upon request, store array with exciton coefficients
  !-----------------------------------------------------
  if( (input%xs%storeexcitons%minnumberexcitons .lt. 1) .or. &
    & (input%xs%storeexcitons%minnumberexcitons .gt. nexc) .or. &
    & (input%xs%storeexcitons%maxnumberexcitons .lt. 1) .or. &
    & (input%xs%storeexcitons%maxnumberexcitons .gt. nexc) .or. &
    & (input%xs%storeexcitons%minnumberexcitons .gt. &
    & input%xs%storeexcitons%maxnumberexcitons) ) then

    write(*,*)
    write(*,'("Error(bse): Wrong range of exciton indices: ", 2i5)') &
      & input%xs%storeexcitons%minnumberexcitons,&
      & input%xs%storeexcitons%maxnumberexcitons
    write(*,*)

    stop
  end if  

  ! Write bin
  call getunit(unexc)
  open(unexc, file='EXCCOEFF.bin', action='write',form='unformatted', iostat=stat)

  if((stat/=0) .and. (rank==0)) then
    write(*,*) stat
    write(*,'("Error(bse): Error creating EXCCOEFF.bin")')
    write(*,*)
    stop
  end if

!! Discuss with author
  ! Write
  write(unexc) input%xs%storeexcitons%minnumberexcitons,&
    & input%xs%storeexcitons%maxnumberexcitons,&
    & nkptnr, iuref, bcou%il1, bcou%il2, bcou%n1, bcou%n2, hamsize, smap

  do ievec = input%xs%storeexcitons%minnumberexcitons,&
    & input%xs%storeexcitons%maxnumberexcitons

    write(unexc) beval(ievec), bevec(1:hamsize,ievec)
  end do

  close(unexc)

end subroutine storeexcitons
