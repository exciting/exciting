! Copyright(C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_getapwcmt

  implicit none

  contains

    subroutine getapwcmt(iq, ik, isti, istf, lmax, apwlm, fname, vpl)
      use modmpi
      use modmain
      use modinput
      use modxs
      use m_getunit

      implicit none

      ! arguments
      integer, intent(in) :: iq, ik, isti, istf, lmax
      complex(8), intent(out) :: apwlm(:, :, :, :)
      character(*), intent(in), optional :: fname
      real(8), intent(in), optional :: vpl(3)

      ! local variables
      character(*), parameter :: thisnam = 'getapwcmt'
      character(256) :: filextt
      character(256) :: filename
      integer :: un, reclen, nerr, nstfv_, apwordmax_, lmaxapw_
      real(8) :: vql_(3), vkl_(3), vklt(3), vqlt(3)
      complex(8), allocatable :: apwlmt(:, :, :, :)
      real(8), external :: r3dist

      nerr = 0

      ! Check band range
      if(isti .lt. 1 .or. istf .gt. nstfv .or. istf .le. isti) then

        write(unitout,*)
        write(unitout, '("Error(getapwcmt): Inconsistent limits for bands:")')
        write(unitout, '("  Band limits	: ", 2i6)') isti, istf
        write(unitout, '("  Maximum value : ", i6)') nstfv
        write(unitout,*)

        call flushifc(unitout)
        nerr = nerr + 1
      end if

      if(size(apwlm, 1) .ne. (istf-isti+1)) then

        write(unitout,*)
        write(unitout, '("Error(getapwcmt): Output array does not match for bands:")')
        write(unitout, '("  Band limits		    : ", 2i6)') isti, istf
        write(unitout, '("  Requested number of bands : ", i6)') istf - isti + 1
        write(unitout, '("  Array size		    : ", i6)') size(apwlm, 1)
        write(unitout,*)

        call flushifc(unitout)
        nerr = nerr + 1
      end if

      ! Check lmax value
      if(lmax .gt. input%groundstate%lmaxapw .or. lmax .lt. 0) then

        write(unitout,*)
        write(unitout, '(a, i8)') 'Error(' // thisnam // '):&
          & lmax > input%groundstate%lmaxapw or lmax < 0:', lmax
        write(unitout,*)

        call flushifc(unitout)
        nerr = nerr + 1
      end if

      if(nerr .gt. 0) call terminate

      ! Generate file name from which to read the APWCMTs
      if(present(fname)) then 
        filename = trim(adjustl(fname))
      else
        ! Set file extension
        filextt = filext
        if(iq .eq. 0) then
          call genfilextread(task)
        end if
        filename = 'APWCMT'//trim(filext)
      end if

      !------------------------!
      !     get parameters     !
      !------------------------!

      inquire(iolength=reclen) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_

      call getunit(un)

      open(un, file=trim(adjustl(filename)), action='read',&
        & form='unformatted', status='old', access='direct', recl=reclen)
      read(un, rec=1) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_

      close(un)

      nerr = 0
      ! Check number of bands
      if(nstfv .gt. nstfv_) then

        write(unitout,*)
        write(unitout, '("Error(", a, "): Invalid nstfv for k-point ", i8)')&
          & thisnam, ik
        write(unitout, '("  q-point    : ", i8)') iq
        write(unitout, '("  current    : ", i8)') nstfv
        write(unitout, '("  file       : ", i8)') nstfv_
        write(unitout, '("  filename   : ", a )') trim(adjustl(filename))
        write(unitout,*)

        call flushifc(unitout)
        nerr = nerr + 1
      end if

      ! Check APW matching order
      if(apwordmax .ne. apwordmax_) then

        write(unitout,*)
        write(unitout, '("Error(", a, "): Invalid apwordmax for k-point ", i8)')&
          & thisnam, ik
        write(unitout, '("  q-point    : ", i8)') iq
        write(unitout, '("  current    : ", i8)') apwordmax
        write(unitout, '("  file       : ", i8)') apwordmax_
        write(unitout, '("  filename   : ", a )') trim(adjustl(filename))
        write(unitout,*)

        call flushifc(unitout)
        nerr = nerr + 1
      end if

      ! Check lmax
      if(input%groundstate%lmaxapw .gt. lmaxapw_) then

        write(unitout,*)
        write(unitout, '("error(", a, "): invalid lmaxapw for k-point ", i8)')&
          & thisnam, ik
        write(unitout, '(" q-point    : ", i8)') iq
        write(unitout, '(" current    : ", i8)') input%groundstate%lmaxapw
        write(unitout, '(" file	      : ", i8)') lmaxapw_
        write(unitout, '(" filename   : ", a )') trim(adjustl(filename))

        call flushifc(unitout)
        write(unitout,*)
        nerr = nerr + 1
      end if

      if(nerr .gt. 0) call terminate

      !------------------!
      !     get data     !
      !------------------!

      ! Assign to output array and apply cutoff
      allocate(apwlmt(nstfv_, apwordmax, (lmaxapw_+1)**2, natmtot))

      ! Read data from file
      inquire(iolength=reclen) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_, apwlmt

      call getunit(un)
      open(un, file=trim(adjustl(filename)), action='read',&
        & form='unformatted', status='old', access='direct', recl=reclen)
      read(un, rec=ik) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_, apwlmt
      close(un)

      ! Check q-point and k-point
      if(input%xs%bse%beyond) then 
        if(present(vpl)) then 
          vklt(:) = vpl
        else
          vklt(:) = vkl(:,ik)
        end if
        vqlt(:) = vqlmt(:, iq)
      else
        if(iq .eq. 0) then
          ! Gamma q-point
          vklt(:) = vkl0(:, ik)
          vqlt(:) = 0.d0
        else
          vklt(:) = vkl(:, ik)
          vqlt(:) = vql(:, iq)
        end if
      end if

      if( r3dist(vkl_, vklt) .gt. input%structure%epslat .or.&
        & (r3dist(vql_, vqlt) .gt. input%structure%epslat .and. .not. tscreen) .or.&
        & (r3dist(vql_, vqlt) .gt. input%structure%epslat .and. input%xs%bse%beyond ))&
        & then

        write(unitout,*)
        write(unitout, '(a)') 'Error(' // thisnam // '):&
          &  Differring parameters for APW MT coefficients (current/file):'
        write(unitout, '(a, i6)') '  q-point index :', iq
        write(unitout, '(a, i6)') '  k-point index :', ik
        write(unitout, '(a, 3f12.6, a, 3f12.6)') ' vql :', vqlt, ', ', vql_
        write(unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl :', vklt, ', ', vkl_
        write(unitout, '(a)') ' file		 : '//trim(adjustl(filename))
        write(unitout,*)

        call flushifc(unitout)
        call terminate
      end if

      ! Retrieve data within cutoffs
      apwlm(:, :, :, :) = apwlmt(isti:istf, :, 1:(lmax+1)**2, :)
      deallocate(apwlmt)

      ! restore file extension
      if(.not. present(fname)) then 
        filext = filextt
      end if
    end subroutine getapwcmt

end module m_getapwcmt
