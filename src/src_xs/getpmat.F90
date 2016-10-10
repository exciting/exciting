! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_getpmat

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
      use m_getunit
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
      integer, intent(in) :: ik, i1, f1, i2, f2
      real(8), intent(in) :: vklt(:, :)
      logical, intent(in) :: tarec
      character(*), intent(in) :: filnam
      complex(8), intent(out) :: pm(:, :, :)

      ! Local variables
      character(*), parameter :: thisnam = 'getpmat'
      integer :: reclen, un, ikr, nstsv_, ierror
      real(8) :: vkl_ (3)
      logical :: existent
      complex(8), allocatable :: pmt(:, :, :)

      ! Functions
      real(8), external :: r3dist

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
      call getunit(un)
      inquire(iolength=reclen) vkl_, nstsv_
      open(unit=un, file=trim(filnam), form='unformatted',&
       & action='read', access='direct', recl=reclen)
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

      call getunit(un)
      ! I/o record length
      inquire(iolength=reclen) vkl_, nstsv_, pmt
      open(unit=un, file=trim(filnam), form='unformatted',&
        & action='read', access='direct', recl=reclen)
      ! Read from file
      read(un, rec=ikr) vkl_, nstsv_, pmt
      close(un)

      ! Check k-point
      if(r3dist(vkl_, vklt(1, ik)) .gt. input%structure%epslat) then
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
    
end module m_getpmat
