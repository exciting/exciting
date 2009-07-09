


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: kernxc_bse2
! !INTERFACE:


subroutine kernxc_bse2
! !USES:
use modinput
  use modmain
  use modxs
  use m_getevalsvr
  use m_getpemat
  use m_xsgauntgen
  use m_findgntn0
  use m_genwgrid
  use m_xszoutpr3
  use m_genfilname
  use m_getunit
  use m_xszoutpr3
! !INPUT/OUTPUT PARAMETERS:
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
!
! !REVISION HISTORY:
!   Created February 2009 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='kernxc_bse2'
  real(8), parameter :: eps=1.d-5
  integer, parameter :: iqmt=1, noptc=3
  character(256) :: filnam
  integer :: iknr, jknr, iv, ic, jv, jc, si, sj, nv, nc, wsiz, n, nt, i, j, iop, iw, un, recl
  integer :: nvl, ncu
  real(8) :: t1, dei, dej
  complex(8) :: zt1
  complex(8), allocatable :: w(:), wmat(:, :), wmatq(:, :), wm(:, :, :, :)
  complex(8), allocatable :: resr(:, :), resq(:, :), oscr(:, :), oscq(:, :)
  complex(8), allocatable :: denr(:), denq(:), fxc(:, :, :)
  complex(8), allocatable :: xiout(:, :, :), pmout(:, :, :), me(:, :)
  real(8), allocatable :: ev(:), de(:), scisk(:, :)
  integer, allocatable :: widx(:, :, :)
  integer, external :: idxkkp

  write(*, *) 'initializing...'

  input%xs%emattype=2
  call init0
  call init1
  call init2

  write(*, *) 'preparing...'

  call readfermi
  call xssave0
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax), input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  call findgntn0(max(input%xs%lmaxapwwf, lolmax), max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  call findocclims(0, istocc0, istocc, istunocc0, istunocc, isto0, isto, istu0, istu)
  call ematbdcmbs(input%xs%emattype)
  allocate(w(nwdf))
  call genwgrid(nwdf, wdos, input%xs%tddft%acont, 0.d0, w_cmplx=w)

  write(*, *) 'done.'

  n=ngq(iqmt)
  nv=nst1
  nc=nst3
  wsiz=nv*nc*nkptnr
  ! same window for bands as bse-routine
  nvl=nv-nbfbse+1
  ncu=nafbse

  write(*, *) 'n, nv, nc, nvl, ncu, wsiz', n, nv, nc, nvl, ncu, wsiz

  allocate(wm(nv, nc, nv, nc), widx(nv, nc, nkptnr))
  allocate(de(wsiz), wmat(wsiz, wsiz), wmatq(wsiz, wsiz), me(-3:n, wsiz))
  allocate(resr(-3:n, wsiz), resq(-3:n, wsiz), oscr(-3:n, -3:n), oscq(-3:n, -3:n))
  allocate(ev(nstsv), denr(nwdf), denq(nwdf))
  allocate(fxc(-3:n, -3:n, nwdf))
  allocate(xiout(nv, nc, n), pmout(3, nv, nc), scisk(nst1, nst3))
  if (allocated(pmou)) deallocate(pmou)
  allocate(pmou(3, nv, nc))
  if (allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3, nc, nv))
  if (allocated(deou)) deallocate(deou)
  allocate(deou(nst1, nst3))
  if (allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst3, nst1))
  if (allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1, nst3))
  if (allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst3, nst1))


  ! set up indices
  si=0
  do iknr=1, nkptnr
  do iv=1, nv
  do ic=1, nc
    si=si+1
    widx(iv, ic, iknr)=si
  end do
  end do
  end do



  ! set up energies and their differences
  do iknr=1, nkptnr
    call getevalsvr('EVALSV_SCR.OUT', 1, nstsv, vkl(:, iknr), ev)
    do iv=1, nv
    do ic=1, nc
      si=widx(iv, ic, iknr)
      de(si)=ev(istocc+ic)-ev(iv)
    end do
    end do
  end do


  write(*, *) 'calculating matrix elements....'

  ! calculate matrix elements
  input%xs%emattype=1
  call ematbdcmbs(input%xs%emattype)
  call ematrad(iqmt)
  call ematqalloc
  do iknr=1, nkptnr
write(unitout, *) 'matrix elements iknr=', iknr
call flushifc(unitout)
    call ematqk1(iqmt, iknr)
    call getdevaldoccsv(iqmt, iknr, iknr, istl1, istu1, istl2, istu2, deou, &
	  docc12, scisk)
    call getpemat(iqmt, iknr, 'PMAT_SCR.OUT', '', m12=xiout, p12=pmout)
    do iv=1, nv
      do ic=1, nc
	si=widx(iv, ic, iknr)
	do iop=1, noptc
	  me(-iop, si)=pmout(iop, iv, ic)
	end do
	me(1:, si)=xiout(iv, ic, :)
      end do
    end do
  end do
  input%xs%emattype=2
  call ematbdcmbs(input%xs%emattype)


  write(*, *) 'done.'

  if ((input%xs%tddft%fxctypenumber.eq.7).or.(input%xs%tddft%fxctypenumber.eq.8)) then
     call getbsediag
     write(unitout, '("Info(", a, "): read diagonal of BSE kernel")') trim(thisnam)
     write(unitout, '(" mean value : ", 2g18.10)') bsed
  end if

  ! set up W-matrix
  wmat(:, :)=zzero
  wmatq(:, :)=zzero
  do iknr=1, nkptnr; do jknr=iknr, nkptnr
    call getbsemat('SCCLI.OUT', idxkkp(iknr, jknr, nkptnr), nv, nc, wm)
    do iv=nvl, nv;  do ic=1, ncu
      si=widx(iv, ic, iknr)
      dei=de(si)
      do jv=nvl, nv; do jc=1, ncu
	sj=widx(jv, jc, jknr)
	if (si.ne.sj) then
	  dej=de(sj)
	  zt1=-wm(iv, ic, jv, jc)	! correct
! * all those cases below lead to unreasonable results *
!zt1=wm(iv,ic,jv,jc)
!zt1=conjg(wm(iv,ic,jv,jc))
!zt1=-conjg(wm(iv,ic,jv,jc))
!zt1=wm(jv,jc,iv,ic)
!zt1=-wm(jv,jc,iv,ic)
!zt1=conjg(wm(jv,jc,iv,ic))
!zt1=-conjg(wm(jv,jc,iv,ic))
	  if (abs(dei-dej).lt.input%xs%tddft%fxcbsesplit) then
	    wmatq(si, sj)=zt1
	  else
	    wmat(si, sj)=zt1/(dei-dej)
	  end if
	end if
      end do; end do
    end do; end do
  end do; end do

  do si=1, wsiz
    do sj=si+1, wsiz
      wmat(sj, si)=-conjg(wmat(si, sj))
      wmatq(sj, si)=conjg(wmatq(si, sj))
    end do
  end do

  write(*, *) 'setting up residuals...'

  ! set up residuals
  resr=transpose(matmul(wmat, conjg(transpose(me))))
  resq=transpose(matmul(wmatq, conjg(transpose(me))))

  ! deallocate the wmat arrays
  deallocate(wmat, wmatq)

  write(*, *) 'setting up kernel...'

  t1=2.d0/(nkptnr*omega)
  fxc(:, :, :)=zzero
  nt=n+noptc+1
  do si=1, wsiz
    oscr(:, :)=zzero
    oscq(:, :)=zzero
    call ZGERU(nt, nt, zone, me(:, si), 1, resr(:, si), 1, oscr, nt)
    call ZGERU(nt, nt, zone, me(:, si), 1, resq(:, si), 1, oscq, nt)
    ! add Hermitian transpose
    forall(i=-3:n, j=-3:n)
      oscr(i, j)=oscr(i, j)+conjg(oscr(j, i))
    end forall
    denr(:)=1.d0/(w(:)+bsed-de(si)+zi*input%xs%broad)
    denq(:)=denr(:)**2
    do iw=1, nwdf
      fxc(:, :, iw)=fxc(:, :, iw)+t1*denr(iw)*oscr(:, :)+t1*denq(iw)*oscq(:, :)
    end do
  end do


  write(*, *) 'writing out kernel...'


  ! write out kernel
  inquire(iolength = recl) n, fxc( - 3: - 1, -3: - 1, 1), fxc( - 3: - 1, 1:, 1), fxc(1:, -3: - 1, 1), fxc(1:, 1:, 1)
  call genfilname(basename = 'FXC_BSE', asc = .false., bzsampl = 0, &
       acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresdf, iqmt = iqmt, filnam = filnam)
  call getunit(un)
  open(un, file = trim(filnam), form = 'unformatted', action = 'write', &
       status = 'replace', access = 'direct', recl = recl)
  do iw=1, nwdf
     write(un, rec = iw) n, fxc( - 3: - 1, -3: - 1, iw), fxc( - 3: - 1, 1:, iw), fxc(1:, -3: - 1, iw), fxc(1:, 1:, iw)
  end do
  close(un)

  deallocate(wm, widx)
  deallocate(de, me)
  deallocate(resr, resq, oscr, oscq)
  deallocate(ev, denr, denq)
  deallocate(xiout, pmout, scisk)

end subroutine kernxc_bse2
!EOC
