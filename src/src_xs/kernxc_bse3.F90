



! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: kernxc_bse3
! !INTERFACE:


subroutine kernxc_bse3
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
!   Created March 2009 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='kernxc_bse3'
  real(8), parameter :: eps=1.d-5
  integer, parameter :: iqmt=1, noptc=3
  character(256) :: filnam
  integer :: iknr, jknr, iv, ic, jv, jc, si, sj, nv, nc, wsiz, n, iop, iw, un, recl
  integer :: nvl, ncu
  complex(8), allocatable :: w(:), mat(:, :), wmat(:, :), wm(:, :, :, :)
  complex(8), allocatable :: l0mat(:, :), l0mata(:, :)
  complex(8), allocatable :: hmat(:, :), hmat2(:, :)
  complex(8), allocatable :: fxc(:, :, :)
  complex(8), allocatable :: xiout(:, :, :), pmout(:, :, :)
  complex(8), allocatable :: xiuot(:, :, :), pmuot(:, :, :)
  complex(8), allocatable :: me(:, :), mea(:, :)
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
  ! initialize states below and above the Fermi energy
  call initocc(nbfbse, nafbse)
  call xssave0
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax), input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  call findgntn0(max(input%xs%lmaxapwwf, lolmax), max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  call findocclims(0, istocc0, istocc, istunocc0, istunocc, isto0, isto, istu0, istu)
  call ematbdcmbs(input%xs%emattype)
  allocate(w(nwdf))
  call genwgrid(nwdf, input%xs%dosWindow%intv, input%xs%tddft%acont, 0.d0, w_cmplx=w)

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
  allocate(de(wsiz), wmat(wsiz, wsiz), mat(wsiz, wsiz))
  allocate(l0mat(wsiz, wsiz), me(-3:n, wsiz))
  allocate(l0mata(wsiz, wsiz), mea(-3:n, wsiz))
  allocate(hmat(-3:n, wsiz), hmat2(wsiz, -3:n))
  allocate(ev(nstsv))
  allocate(fxc(-3:n, -3:n, nwdf))
  allocate(scisk(nst1, nst3))
  allocate(xiout(nv, nc, n), pmout(3, nv, nc))
  allocate(xiuot(nc, nv, n), pmuot(3, nc, nv))
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
	  mea(-iop, si)=pmuot(iop, ic, iv)
	    end do
	    me(1:, si)=xiout(iv, ic, :)
	    mea(1:, si)=xiuot(ic, iv, :)
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

  ! set up W-matrix and L^0 matrix
  wmat(:, :)=zzero
  do iknr=1, nkptnr
  do jknr=iknr, nkptnr
    call getbsemat('SCCLI.OUT', idxkkp(iknr, jknr, nkptnr), nv, nc, wm)
    do iv=nvl, nv
    do ic=1, ncu
      si=widx(iv, ic, iknr)
      do jv=nvl, nv
      do jc=1, ncu
	sj=widx(jv, jc, jknr)
	    if (si.ne.sj) then
	      wmat(si, sj)=-wm(iv, ic, jv, jc)
	    end if
      end do
      end do
    end do
    end do
  end do
  end do

  do si=1, wsiz
    do sj=si+1, wsiz
      wmat(sj, si)=conjg(wmat(si, sj))
    end do
  end do

  ! loop over w-points
  do iw=1, nwdf

    write(*, *) 'w-loop, at: ', iw

    ! set up L0-matrix
    l0mat(:, :)=zzero
    l0mata(:, :)=zzero
    do si=1, wsiz
      l0mat(si, si)=1.d0/(w(iw)+bsed-de(si)+zi*input%xs%broad)
      l0mata(si, si)=-1.d0/(w(iw)+bsed+de(si)+zi*input%xs%broad)
    end do

!    ! do the 4 matrix multiplications
!    hmat=matmul(me,l0mat)
!    hmat=matmul(hmat,wmat)
!    hmat=matmul(hmat,l0mat)
!    fxc(:,:,iw)=2.d0*matmul(hmat,conjg(transpose(me)))

    ! resonant contribution
    do si=1, wsiz
      hmat(:, si)=me(:, si)*l0mat(si, si)
    end do
    hmat=matmul(hmat, wmat)
    do si=1, wsiz
      hmat2(si, :)=l0mat(si, si)*conjg(me(:, si))
    end do
    fxc(:, :, iw)=2.d0*matmul(hmat, hmat2)/(nkptnr*omega)

   if (input%xs%tddft%aresfxc) then
    ! anti-resonant contribution
    do si=1, wsiz
      hmat(:, si)=mea(:, si)*l0mata(si, si)
    end do
    hmat=matmul(hmat, -conjg(wmat))
    do si=1, wsiz
      hmat2(si, :)=l0mata(si, si)*conjg(mea(:, si))
    end do
    fxc(:, :, iw)=fxc(:, :, iw) + 2.d0*matmul(hmat, hmat2)/(nkptnr*omega)
   end if

  end do


   ! deallocate the wmat arrays
  deallocate(wmat, mat)

  write(*, *) 'writing out kernel...'


  ! write out kernel
  inquire(iolength=recl) n, fxc(-3:-1, -3:-1, 1), fxc(-3:-1, 1:, 1), fxc(1:, -3:-1, 1), fxc(1:, 1:, 1)
  call genfilname(basename = 'FXC_BSE', asc = .false., bzsampl = 0, &
       acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresfxc, tord=input%xs%tddft%tordfxc, iqmt = iqmt,&
    &filnam = filnam)
  call getunit(un)
  open(un, file = trim(filnam), form = 'unformatted', action = 'write', &
       status = 'replace', access = 'direct', recl = recl)
  do iw=1, nwdf
     write(un, rec=iw) n, fxc(-3:-1, -3:-1, iw), fxc(-3:-1, 1:, iw), fxc(1:, -3:-1, iw), fxc(1:, 1:, iw)
write(8888, '(i6, 6g18.10)') iw, fxc(-1, -1, iw), fxc(-2, -2, iw), fxc(-3, -3, iw)
  end do
  close(un)

  deallocate(wm, widx)
  deallocate(de, me)
  deallocate(ev, hmat, hmat2)
  deallocate(xiout, pmout, xiuot, pmuot, scisk, l0mat, l0mata)

end subroutine kernxc_bse3
!EOC
