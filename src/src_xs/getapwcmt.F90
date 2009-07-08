

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getapwcmt
  implicit none
contains

  ! APW functions


subroutine getapwcmt(iq, ik, isti, istf, lmax, apwlm)
    use modmain
use modinput
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq, ik, isti, istf, lmax
    complex(8), intent(out) :: apwlm(:, :, :, :)
    ! local variables
    character(*), parameter :: thisnam='getapwcmt'
    character(256) :: filextt
    integer :: un, recl, err, nstfv_, apwordmax_, lmaxapw_
    real(8) :: vql_(3), vkl_(3), vklt(3), vqlt(3)
    complex(8), allocatable :: apwlmt(:, :, :, :)
    real(8), external :: r3dist
    err=0
    ! check band range
    if ((isti.lt.1).or.(istf.gt.nstfv).or.(istf.le.isti)) then
       write(unitout, *)
       write(unitout, '("Error(getapwcmt): inconsistent limits for bands:")')
       write(unitout, '(" band limits	: ", 2i6)') isti, istf
       write(unitout, '(" maximum value : ", i6)') nstfv
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if (size(apwlm, 1).ne.(istf-isti+1)) then
       write(unitout, *)
       write(unitout, '("Error(getapwcmt): output array does not match for &
	    &bands:")')
       write(unitout, '(" band limits		    : ", 2i6)') isti, istf
       write(unitout, '(" requested number of bands : ", i6)') istf-isti+1
       write(unitout, '(" array size		    : ", i6)') size(apwlm, 1)
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    ! check lmax value
    if ((lmax.gt.input%groundstate%lmaxapw).or.(lmax.lt.0)) then
       write(unitout, *)
       write(unitout, '(a, i8)') 'Error('//thisnam//'): lmax > input%groundstate%lmaxapw or < 0:', &
	    lmax
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    if (err.gt.0) call terminate
    ! set file extension
    filextt=filext
    if (iq.eq.0) call genfilextread(task)
    !------------------------!
    !     get parameters     !
    !------------------------!
    inquire(iolength=recl) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_
    call getunit(un)
    open(un, file = 'APWCMT'//trim(filext), action = 'read', &
	 form = 'unformatted', status = 'old', access = 'direct', recl = recl)
    read(un, rec=1) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_
    close(un)
    err=0
    ! check number of bands
    if (nstfv.gt.nstfv_) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): invalid nstfv for k-point ", I8)') &
	    thisnam, ik
       write(unitout, '(" q-point    : ", I8)') iq
       write(unitout, '(" current    : ", I8)') nstfv
       write(unitout, '(" FILE	     : ", I8)') nstfv_
       write(unitout, '(" filename   : ", a )') 'APWCMT'//trim(filext)
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    ! check APW matching order
    if (apwordmax.ne.apwordmax_) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): invalid apwordmax for k-point ", I8)') &
	    thisnam, ik
       write(unitout, '(" q-point    : ", I8)') iq
       write(unitout, '(" current    : ", I8)') apwordmax
       write(unitout, '(" FILE	     : ", I8)') apwordmax_
       write(unitout, '(" filename   : ", a )') 'APWCMT'//trim(filext)
       write(unitout, *)
       call flushifc(unitout)
       err=err+1
    end if
    ! check lmax
    if (input%groundstate%lmaxapw.gt.lmaxapw_) then
       write(unitout, *)
       write(unitout, '("Error(", a, "): invalid lmaxapw for k-point ", I8)') &
	    thisnam, ik
       write(unitout, '(" q-point    : ", I8)') iq
       write(unitout, '(" current    : ", I8)') input%groundstate%lmaxapw
       write(unitout, '(" FILE	     : ", I8)') lmaxapw_
       write(unitout, '(" filename   : ", a )') 'APWCMT'//trim(filext)
       call flushifc(unitout)
       write(unitout, *)
       err=err+1
    end if
    if (err.gt.0) call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! assign to output array and apply cutoff
    allocate(apwlmt(nstfv_, apwordmax, (lmaxapw_+1)**2, natmtot))
    ! read data from file
    inquire(iolength=recl) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_, apwlmt
    call getunit(un)
    open(un, file = 'APWCMT'//trim(filext), action = 'read', form = 'unformatted', &
	 status = 'old', access = 'direct', recl = recl)
    read(un, rec=ik) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_, apwlmt
    close(un)
    ! check q-point and k-point
    if (iq.eq.0) then
       ! Gamma Q-point
       vklt(:)=vkl0(:, ik)
       vqlt(:)=0.d0
    else
       vklt(:)=vkl(:, ik)
       vqlt(:)=vql(:, iq)
    end if
    if ((r3dist(vkl_, vklt).gt.input%structure%epslat).or.((r3dist(vql_,&
    &vqlt).gt.input%structure%epslat).and.(.not.tscreen))) then
       write(unitout, *)
       write(unitout, '(a)') 'Error('//thisnam//'): differring parameters for &
	    &APW MT coefficients (current/file): '
       write(unitout, '(a, i6)') ' q-point index  :', iq
       write(unitout, '(a, i6)') ' k-point index  :', ik
       write(unitout, '(a, 3f12.6, a, 3f12.6)') ' vql		 :', vqlt, ', ', vql_
       write(unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt, ', ', vkl_
       write(unitout, '(a)')	' file		 : APWCMT'//trim(filext)
       write(unitout, *)
       call flushifc(unitout)
       call terminate
    end if
    ! retrieve data within cutoffs
    apwlm(:, :, :, :)=apwlmt(isti:istf, :, 1:(lmax+1)**2, :)
    deallocate(apwlmt)
    ! restore file extension
    filext=filextt
  end subroutine getapwcmt

end module m_getapwcmt
