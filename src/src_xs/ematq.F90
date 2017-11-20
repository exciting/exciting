! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematq(iq)
  use modinput
  use modmpi
  use mod_misc, only: filext
  use modxs, only: fnemat, fnemat_t, fnetim, unitout,&
    & ikmapikq, iqmtgamma, filext0,&
    & istocc0, istunocc0, isto0, istu0,&
    & istocc, istunocc, isto, istu,&
    & nst1, istl1, istu1, nst2, istl2, istu2,&
    & nst3, istl3, istu3, nst4, istl4, istu4,&
    & ngq, kpari, kparf, iqmt0, iqmt1,&
    & qvkloff, xiou, xiuo, bcbs
  use m_writegqpts
  use m_filedel
  use m_genfilname
  use m_ematqk
  use m_putemat

  implicit none

  ! arguments
  integer(4), intent (in) :: iq

  ! local variables
  character(*), parameter :: thisnam = 'ematq'
  integer(4) :: ik, n
  type(bcbs) :: bc

  ! filenames
  call genfilname(basename='EMAT', iqmt=iq, &
     & etype=input%xs%emattype, filnam=fnemat)
  call genfilname(basename='EMAT', iqmt=iq, &
     & etype=input%xs%emattype, procs=procs, rank=rank, &
     & filnam=fnemat_t)
  call genfilname(nodotpar=.true., basename='EMAT_TIMING', &
     & iqmt=iq, etype=input%xs%emattype, procs=procs, rank=rank, &
     & filnam=fnetim)

  ! file extension for gamma-point and q-point
  call genfilname(iqmt=iqmtgamma, fileext=filext0)
  iqmt0 = iqmtgamma
  call genfilname(iqmt=iq, setfilext=.true.)
  iqmt1 = iq

  ! calculate k+q and g+k+q related variables
  call init1offs(qvkloff(1:3, iq))

  ! write g+q-vectors
  if(rank .eq. 0) then
    call writegqpts(iq, filext)
    call writekmapkq(iq)
  end if

  ! find highest (partially) occupied and lowest (partially) unoccupied states
  call findocclims(iq, ikmapikq(:,iq), istocc0, istunocc0, isto0, isto, istu0, istu)
  istunocc = istunocc0
  istocc = istocc0

  call ematbdlims(1, nst1, istl1, istu1, nst2, istl2, istu2)

  ! generate radial integrals wrt. sph. bessel functions
  call ematrad(iq)

  ! delete timing information of previous runs
  call filedel(trim(fnetim))

  ! write information
  write(unitout, '(a, i6)') 'info(' // thisnam // '):&
    & Number of g+q vectors:', ngq(iq)

  ! Matrix size for response function
  ! Get number of G+q vectors for current q
  n = ngq(iq)

  call ematqalloc

  ! loop over k-points at fixed q, compute
  ! both ou and uo plane wave matrix elements and store them
  ! to file
  do ik = kpari, kparf

    ! The plane wave elements for ou and uo transitions are 
    ! calculated and stored in xiou and xiuo
    ! Set 12=ou 34=uo
    call ematbdcmbs(1)

    ! Get ou
    if(allocated(xiou)) deallocate(xiou)
    allocate(xiou(nst1, nst2, n))
    bc%n1 = nst1
    bc%n2 = nst2
    bc%il1 = istl1
    bc%il2 = istl2
    bc%iu1 = istu1
    bc%iu2 = istu2
    ikmapikq_ptr => ikmapikq
    call setptr01
    call ematqk(iq, ik, xiou, bc)

    ! Get uo
    if(allocated(xiuo)) deallocate(xiuo)
    allocate(xiuo(nst3, nst4, n))
    bc%n1 = nst3
    bc%n2 = nst4
    bc%il1 = istl3
    bc%il2 = istl4
    bc%iu1 = istu3
    bc%iu2 = istu4
    ikmapikq_ptr => ikmapikq
    call setptr01
    call ematqk(iq, ik, xiuo, bc)

    ! Store to file
    call putemat(iq, ik, .True.,&
      & trim(fnemat), istl1, istu1, istl2, istu2, xiou,&
      & istl3, istu3, istl4, istu4, xiuo)

  end do

  call ematqdealloc

  ! sharedfs defaults to true
  if(.not. input%sharedfs) call cpFileToNodes(fnemat)

  call barrier

end subroutine ematq
