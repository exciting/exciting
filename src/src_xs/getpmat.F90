! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_getpmat
    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    implicit none

  contains

    !BOP
    ! !ROUTINE: getpmat
    ! !INTERFACE:
    subroutine getpmat(ik, vklt, i1, f1, i2, f2, tarec, filnam, pm)
    ! !USES:
      use modmain
      use modinput
      use modxs
      use modmpi
      use precision, only: i32, long_int, dp
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer :: ik           ! k point index
    ! integer :: i1,f1        ! Range of bands (final states)
    ! integer :: i2,f2        ! Range of bands (initial states)
    ! logical :: tarec        ! Use absolute ik as record index or
    !                         ! relative ik w.r.t k set of currrent MPI rank
    ! character(*) :: filname ! File name of direct access momentum matrix file
    ! Out:
    ! complex(8) :: pm(3, i1:f1, i2:f2) ! The requested slice of the saved momentum
    !                                   ! matrix elements
    ! !DESCRIPTION:
    !   Reads selectable portions of the momentum matrix form file.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      integer(i32), intent(in) :: ik, i1, f1, i2, f2
      real(dp), intent(in) :: vklt(:, :)
      logical(i32), intent(in) :: tarec
      character(*), intent(in) :: filnam
      complex(dp), intent(out) :: pm(:, :, :)

      ! Local variables
      character(*), parameter :: thisnam = 'getpmat'
      integer(i32) :: un, ikr, nstsv_, ierror
      integer(long_int) :: reclen
      real(dp) :: vkl_ (3)
      logical(i32) :: existent
      complex(dp), allocatable :: pmt(:, :, :)

      ! Functions
      real(dp), external :: r3dist

      ! Check if file exists
      inquire(file=trim(filnam), exist=existent)
      if( .not. existent) then
        write(unitout, '(a)') 'Error(' // thisnam // '):&
          & file does not exist: ' // trim(filnam)
        call terminate
      end if

      ! Record position for k-point
      ikr = ik
      if( .not. tarec) then 
        ! Get relative (to set of current MPI rank) k index
        call getridx(nkpt, ik, ikr)
      end if

      !!<-- Checking range limits
      ! Check band range
      ierror = 0
      if( (i1 .lt. 1) .or. (i1 .gt. nstsv)&
        & .or. (f1 .lt. 1) .or. (f1 .gt. nstsv)&
        & .or. (i2 .lt. 1) .or. (i2 .gt. nstsv)&
        & .or. (f2 .lt. 1) .or. (f2 .gt. nstsv)&
        & .or. (i1 .gt. f1) .or. (i2 .gt. f2)) then

        write(unitout,*)
        write(unitout, '("Error(", a, "): inconsistent limits for states:")') thisnam
        write(unitout, '(" limits(lo/hi) : ", 2(2i6, 2x))') i1, f1, i2, f2
        write(unitout, '(" maximum value  : ", i6)') nstsv
        write(unitout,*)

        call flushifc(unitout)
        ierror = ierror + 1
      end if

      ! Check compatibility of ranges with passed pm array
      if(size(pm, 2) .ne. (f1-i1+1) .or. size(pm, 3) .ne. (f2-i2+1)) then

        write(unitout,*)
        write(unitout, '("Error(", a, "):&
         & output array does not match for states:")') thisnam
        write(unitout, '(" limits		     : ", 2(2i6, 2x))') i1, f1, i2, f2
        write(unitout, '(" requested number of states : ", 2i6)')&
          & f1 - i1 + 1, f2 - i2 + 1
        write(unitout, '(" array sizes		     : ", 2i6)')&
          & size(pm, 2), size(pm, 3)
        write(unitout,*)

        call flushifc(unitout)
        ierror = ierror + 1
      end if

      if(ierror .gt. 0) call terminate
      !!-->

      !------------------------!
      !     get parameters     !
      !------------------------!

      ! Read saved vkl and nstsv values from file 
      call inquire_large( reclen, vkl_, [nstsv_] )
      call open_direct_unformatted_large( un, trim( filnam ), "read", reclen, "old" )
      read(un, rec=1) vkl_, nstsv_
      close(un)

      ! Check if all states can be read from file
      ierror = 0
      if((f1 .gt. nstsv_) .or. (f2 .gt. nstsv_)) then
        write(unitout,*)
        write(unitout, '("Error(", a, "): requested states out of&
          & range for k-point ", i8)') thisnam, ik
        write(unitout, '(" limits	  : ", 2(2i6, 2x))') i1, f1, i2, f2
        write(unitout, '(" cutoff from file: ", i8)') nstsv_
        write(unitout, '(" filename	  : ", a )') trim(filnam)
        write(unitout,*)
        call flushifc(unitout)
        ierror = ierror + 1
      end if
      if(ierror .gt. 0) call terminate

      !------------------!
      !     get data     !
      !------------------!

      ! Allocate local arrays
      allocate(pmt(3, nstsv_, nstsv_))

      ! I/o record length
      call inquire_large( reclen, vkl_, [nstsv_], pmt )
      call open_direct_unformatted_large( un, trim( filnam ), "read", reclen, "old" )
      ! Read from file
      read(un, rec=ikr) vkl_, nstsv_, pmt
      close(un)

      ! Check k-point
      if(r3dist(vkl_, vklt(:, ik)) .gt. input%structure%epslat) then
        write(unitout,*)
        write(unitout, '(a)') 'Error(' // thisnam // '):&
          & differring parameters for matrix elements(current/file): '
        write(unitout, '(a, i6)') ' k-point index  :', ik
        write(unitout, '(a, i6)') ' record position:', ikr
        write(unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :',&
          & vklt(:, ik), ', ', vkl_
        write(unitout, '(" filename	  : ", a )') trim(filnam)
        write(unitout,*)
        call flushifc(unitout)
        call terminate
      end if

      ! Retrieve data within cutoffs
      pm(:, :, :) = pmt(:, i1:f1, i2:f2)

      deallocate(pmt)
    end subroutine getpmat
    !EOC
       
    subroutine getpmatxas (ik, vklt, i1, f1, i2, f2, tarec, filnam, pm)
         use modmain
         use modinput
         use modxs
         use modmpi
         use modxas, only : ncg
         use precision, only: i32, long_int, dp
         implicit None
    ! arguments
         integer(i32), intent(in) :: ik, i1, f1, i2, f2
         real(dp), intent (in) :: vklt (:, :)
         logical(i32), intent(in) :: tarec
         character(*), intent(in) :: filnam
         complex(dp), intent (out) :: pm (:, :, :)
    ! local variables
         character (*), parameter :: thisnam = 'getpmat'
         integer(long_int) :: recl
         integer(i32) :: un, ikr, nstsv_, err
         real(dp) :: vkl_ (3)
         logical(i32) :: existent
         complex(dp), allocatable :: pmt (:, :, :)
    ! functions
         real(dp), external :: r3dist
    ! check if file exists
         inquire (File=trim(filnam), Exist=existent)
         if ( .Not. existent) then
            write (unitout, '(a)') 'Error(' // thisnam // '): file does&
           & not exist: ' // trim (filnam)
            call terminate
         end If
    ! record position for k-point
         ikr = ik
         if ( .Not. tarec) call getridx(nkpt, ik, ikr)
         err = 0
    ! check band range
         if ((i1 .Lt. 1) .Or. (i1 .Gt. ncg) .Or. (f1 .Lt. 1) .Or. (f1 &
        & .Gt. ncg) .Or. (i2 .Lt. 1) .Or. (i2 .Gt. nstsv) .Or. (f2 &
        & .Lt. 1) .Or. (f2 .Gt. nstsv) .Or. (i1 .Gt. f1) .Or. (i2 .Gt. &
        & f2)) then
            write (unitout,*)
            write (unitout, '("Error(", a, "): inconsistent limits for &
           &states:")') thisnam
            write (unitout, '(" limits (lo/hi) : ", 2(2i6, 2x))') i1, &
           & f1, i2, f2
            write (unitout, '(" maximum value  : ", 2(i6, x))') ncg,nstsv
            write (unitout,*)
            call flushifc (unitout)
            err = err + 1
         end if
         if ((size(pm, 2) .Ne. (f1-i1+1)) .Or. (size(pm, 3) .Ne. &
        & (f2-i2+1))) Then
            write (unitout,*)
            write (unitout, '("Error(", a, "): output array does not ma&
           &tch for states:")') thisnam
            write (unitout, '(" limits		     : ", 2(2i6, 2x))') i1, f1, &
           & i2, f2
            write (unitout, '(" requested number of states : ", 2i6)') &
           & f1 - i1 + 1, f2 - i2 + 1
            write (unitout, '(" array sizes		     : ", 2i6)') size (pm, &
           & 2), size (pm, 3)
            write (unitout,*)
            call flushifc (unitout)
            err = err + 1
         end if
         if (err .Gt. 0) Call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
         call inquire_large( Recl, vkl_, [nstsv_] )
         call open_direct_unformatted_large( un, trim( filnam ), "read", Recl, "old" )
         read (un, Rec=1) vkl_, nstsv_
         close (un)
         err = 0
    ! check if all states can be read from file
         if ((f1 .Gt. ncg) .Or. (f2 .Gt. nstsv_)) Then
            write (unitout,*)
            write (unitout, '("Error(", a, "): requested states out of &
           &range for k-point ", I8)') thisnam, ik
            write (unitout, '(" limits	  : ", 2(2I6, 2x))') i1, f1, i2, &
           & f2
            write (unitout, '(" cutoff from file: ", I8)') nstsv_
            write (unitout, '(" filename	  : ", a )') trim (filnam)
            write (unitout,*)
            call flushifc (unitout)
            err = err + 1
         end if
         if (err .Gt. 0) Call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! allocate local arrays
         allocate (pmt(3, ncg, nstsv_))
    ! I/O record length                  
         call inquire_large( Recl, vkl_, [nstsv_], pmt )
         call open_direct_unformatted_large( un, trim( filnam ), "read", Recl, "old" )
    ! read from file
         read (un, Rec=ikr) vkl_, nstsv_, pmt
         close (un)
    ! check k-point
         if (r3dist(vkl_, vklt(:, ik)) .Gt. input%structure%epslat) &
        & then
            write (unitout,*)
            write (unitout, '(a)') 'Error(' // thisnam // '): differrin&
           &g parameters for matrix elements (current/file): '
            write (unitout, '(a, i6)') ' k-point index  :', ik
            write (unitout, '(a, i6)') ' record position:', ikr
            write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt &
           & (:, ik), ', ', vkl_
            write (unitout, '(" filename	  : ", a )') trim (filnam)
            write (unitout,*)
            call flushifc (unitout)
            call terminate
         end If
    ! retrieve data within cutoffs
         pm (:, :, :) = pmt (:, i1:f1, i2:f2)
         deallocate (pmt)
      end subroutine getpmatxas
   
end module m_getpmat
