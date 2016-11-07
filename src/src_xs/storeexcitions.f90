module m_storeexcitons
  use modinput, only: input
  use mod_kpoint, only: nkptnr
  use modmpi
  use modbse
  use m_getunit

  implicit none

  contains

    subroutine storeexcitons(nexc, beval, bevec, fcoup)

      implicit none

      ! I/O
      integer(4), intent(in) :: nexc
      real(8), intent(in) :: beval(:)
      complex(8), intent(in) :: bevec(:,:)
      logical, intent(in) :: fcoup

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
        call terminate
      end if  

      ! Write bin
      call getunit(unexc)
      open(unexc, file='EXCCOEFF.bin', action='write',form='unformatted', iostat=stat)

      if((stat/=0) .and. (rank==0)) then
        write(*,*) stat
        write(*,'("Error (storeexcitons): Error creating EXCCOEFF.bin")')
        write(*,*)
        call terminate
      end if

    !! Discuss with author
      ! Write
      !   Meta data
      write(unexc)&
        & input%xs%storeexcitons%minnumberexcitons,&
        & input%xs%storeexcitons%maxnumberexcitons,&
        & input%xs%broad,&
        & wl, wu, econv,&
        & nkptnr, nk_bse, kousize, koulims, hamsize,&
        & smap, smap_rel, fcoup 
      !   EVAL and EVEC
      do ievec = input%xs%storeexcitons%minnumberexcitons,&
        & input%xs%storeexcitons%maxnumberexcitons
        if(fcoup) then 
          write(unexc) beval(ievec), bevec(1:2*hamsize,ievec)
        else
          write(unexc) beval(ievec), bevec(1:hamsize,ievec)
        end if
      end do

      close(unexc)

    end subroutine storeexcitons
    
end module m_storeexcitons
