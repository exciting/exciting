



! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine scrcoulint
  use modmain
use modinput
  use modmpi
  use modxs
  use summations
  use m_xsgauntgen
  use m_findgntn0
  use m_writevars
  use m_genfilname
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam='scrcoulint'
  character(256) :: fnsccli, fnscreeninv
  real(8), parameter :: epsortho=1.d-12
  integer :: ikkp,iknr,jknr,iqr,iq,iqrnr,jsym,jsymi,igq1,igq2,n,recl,un
  integer :: nsc,iv(3),ivgsym(3),j1,j2,nkkp
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24
  logical :: tq0,tphf
  real(8) :: vqr(3),vq(3),t1
  integer :: sc(maxsymcrys),ivgsc(3,maxsymcrys)
  integer, allocatable :: igqmap(:)
  complex(8) :: zt1
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:),scclid(:,:)
  complex(8), allocatable :: scieffg(:,:,:),tm(:,:),tmi(:,:),bsedt(:,:)
  complex(8), allocatable :: phf(:,:),emat12(:,:),emat34(:,:)
  ! external functions
  integer, external :: idxkkp
  logical, external :: tqgamma

  !---------------!
  !   main part   !
  !---------------!

  input%xs%emattype=2
  call init0
  call init1
  call init2
  ! read Fermi energy from file
  call readfermi
  ! save variables for the Gamma q-point
  call xssave0
  ! generate Gaunt coefficients
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax), input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  ! find indices for non-zero Gaunt coefficients
  call findgntn0(max(input%xs%lmaxapwwf, lolmax), max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
  write(unitout,'(a,3i8)') 'Info('//thisnam//'): Gaunt coefficients generated &
       &within lmax values:', input%groundstate%lmaxapw, input%xs%lmaxemat, input%groundstate%lmaxapw
  write(unitout, '(a, i6)') 'Info('//thisnam//'): number of q-input%xs%dosWindow%points: ', nqpt
  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  ! only for systems with a gap in energy
  if (.not.ksgap) then
     write(*,*)
     write(*,'("Error(",a,"): screened Coulomb interaction works only for &
          &systems with KS-gap.")') trim(thisnam)
     write(*,*)
     call terminate
  end if
  ! check number of empty states
  if (input%xs%screening%nempty.lt.input%groundstate%nempty) then
     write(*,*)
     write(*,'("Error(",a,"): too few empty states in screening eigenvector &
          &file - the screening should include many empty states &
	  &(BSE/screening)", 2i8)') trim(thisnam), input%groundstate%nempty, input%xs%screening%nempty
     write(*,*)
     call terminate
  end if
  call ematbdcmbs(input%xs%emattype)
  nst12=nst1*nst2
  nst34=nst3*nst4
  nst13=nst1*nst3
  nst24=nst2*nst4
  call genfilname(dotext='_SCI.OUT',setfilext=.true.)
  if (rank.eq.0) then
     call writekpts
     call writeqpts
  end if

  ! local arrays
  allocate(phf(ngqmax,ngqmax))
  allocate(sccli(nst1,nst3,nst2,nst4),scclid(nst1,nst3))
  allocate(scieffg(ngqmax,ngqmax,nqptr))
  sccli(:,:,:,:)=zzero
  scieffg(:,:,:)=zzero

  ! set file extension
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)

  !-----------------------------------!
  !     loop over reduced q-points    !
  !-----------------------------------!
  call getunit(un)
  call genparidxran('q',nqptr)

  do iqr=qpari,qparf
    call genfilname(basename='W_SCREEN',iq=iqr,filnam=fnscreeninv)
	  open(un,file=trim(fnscreeninv),form='formatted',action='write', &
			status='replace')
    call chkpt(3,(/task,1,iqr/),'task,sub,reduced q-point; generate effective &
          &screened Coulomb potential')
     ! locate reduced q-point in non-reduced set
     iqrnr=iqmap(ivqr(1,iqr),ivqr(2,iqr),ivqr(3,iqr))
     n=ngq(iqrnr)

     ! calculate effective screened Coulomb interaction
     call genscclieff(iqr,ngqmax,n,scieffg(1,1,iqr))
     do igq1=1,n
       do igq2=1,n
         write(un,'(2i8,3g18.10)') igq1,igq2,scieffg(igq1,igq2,iqr),abs(scieffg(igq1,igq2,iqr))
       end do
     end do
     call writevars(un,iqr,0)
     close(un)

     ! generate radial integrals for matrix elements of plane wave
     call putematrad(iqr,iqrnr)
  end do
  ! communicate array-parts wrt. q-points
  call zalltoallv(scieffg,ngqmax**2,nqptr)

  call barrier

  ! information on size of output file
  nkkp=(nkptnr*(nkptnr+1))/2
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4,sccli
  write(unitout,*)
  write(unitout,'(a,f12.3)') 'file size for screened Coulomb &
       &interaction (GB):',recl*nkkp/1024.d0**3
  write(unitout,*)

  !-------------------------------!
  !     loop over (k,kp) pairs    !
  !-------------------------------!
  call genfilname(basename='SCCLI',asc=.true.,filnam=fnsccli)
  call getunit(un)
	if (rank.eq.0) open(un,file=trim(fnsccli),form='formatted',action='write', &
		status='replace')
  call genparidxran('p',nkkp)
  allocate(bsedt(3,0:procs-1))
  bsedt(1,:)=1.d8
  bsedt(2,:)=-1.d8
  bsedt(3,:)=zzero
  ! loop over combinations of k-points
  do ikkp=ppari,pparf
     call chkpt(3,(/task,2,ikkp/),'task,sub,(k,kp)-pair; direct term of &
	  &associated(input%xs%BSE) - Hamiltonian')
     call kkpmap(ikkp,nkptnr,iknr,jknr)
     ! k-point difference
     iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
     iv(:)=modulo(iv(:),ngridq(:))
     ! q-point (reduced)
     iqr=iqmapr(iv(1),iv(2),iv(3))
     vqr(:)=vqlr(:,iqr)
     ! q-point (non-reduced)
     iq=iqmap(iv(1),iv(2),iv(3))
     vq(:)=vql(:,iq)
     ! local field effects size
     tq0=tqgamma(iq)
     n=ngq(iq)

     allocate(igqmap(n),emat12(nst12,n),emat34(nst34,n))
     allocate(tm(n,n),tmi(n,n))
     allocate(scclit(nst12,nst34))

     ! find symmetry operations that reduce the q-point to the irreducible
     ! part of the Brillouin zone
     call findsymeqiv(input%xs%BSE%fbzq, vq, vqr, nsc, sc, ivgsc)

     ! find the map that rotates the G-vectors
     call findgqmap(iq,iqr,nsc,sc,ivgsc,n,jsym,jsymi,ivgsym,igqmap)
     ! generate phase factor for dielectric matrix due to non-primitive
     ! translations
     call genphasedm(iq,jsym,ngqmax,n,phf,tphf)

     ! get radial integrals
     call getematrad(iqr,iq)
     ! rotate radial integrals
     call rotematrad(n,igqmap)
     ! rotate inverse of screening, Coulomb potential and radial integrals
     tmi(:,:)=phf(:n,:n)*scieffg(igqmap,igqmap,iqr)

     ! calculate matrix elements of the plane wave
     input%xs%emattype=2
     call ematbdcmbs(input%xs%emattype)
     call ematqalloc
     call ematqk1(iq,iknr)
     input%xs%emattype=2
     call ematbdcmbs(input%xs%emattype)
     call chkpt(3,(/task,2,ikkp/),'task,sub,(k,kp)-pair; direct term of &
	  &associated(input%xs%BSE) - Hamiltonian')

     ! select screening level
     tm(:,:)=zzero
     select case (trim(input%xs%screening%screentype))
     case('longrange')
        ! constant screening (q=0 average tensor)
        forall(igq1=1:n)
           tm(igq1,igq1)=fourpi*dot_product(vgqc(:,igq1,iq),matmul(dielten, &
                vgqc(:,igq1,iq)))
        end forall
     case('diag')
        ! only diagonal of screening
        forall(igq1=1:n)
           tm(igq1,igq1)=tmi(igq1,igq1)
        end forall
     case('full')
        ! full screening
        tm(:,:)=tmi(:,:)
     end select

     ! combine indices for matrix elements of plane wave
     j1=0
     do ist2=1,nst2
        do ist1=1,nst1
           j1=j1+1
           emat12(j1,:)=xiou(ist1,ist2,:)
        end do
     end do
     j2=0
     do ist4=1,nst4
        do ist3=1,nst3
           j2=j2+1
           emat34(j2,:)=xiuo(ist3,ist4,:)
        end do
     end do

     ! matrix elements of direct term (as in BSE-code of Peter and
     ! in the SELF-documentation of Andrea Marini)
     scclit=matmul(conjg(emat12),matmul(tm,transpose(emat34)))/omega/nkptnr

     ! map back to individual band indices
     j2=0
     do ist4=1,nst4
        do ist3=1,nst3
           j2=j2+1
           j1=0
           do ist2=1,nst2
              do ist1=1,nst1
                 j1=j1+1
                 sccli(ist1,ist3,ist2,ist4)=scclit(j1,j2)
              end do
           end do
        end do
     end do
     if ((rank.eq.0).and.(ikkp.le.3)) then
        ! write to ASCII file
        do ist1=1,nst1
           do ist3=1,nst3
              do ist2=1,nst2
                 do ist4=1,nst4
		    zt1=sccli(ist1,ist3,ist2,ist4)
                    write(un,'(i5,3x,3i4,2x,3i4,2x,4e18.10)') ikkp,iknr,ist1,&
                         ist3,jknr,ist2,ist4,zt1,abs(zt1)**2, &
			 atan2(aimag(zt1),dble(zt1))/pi
                 end do
              end do
           end do
        end do
     end if
     ! analyze BSE diagonal
     if (iknr.eq.jknr) then
        do ist1=1,nst1
           do ist3=1,nst3
              zt1=sccli(ist1,ist3,ist1,ist3)
	      scclid(ist1,ist3)=zt1
              t1=dble(zt1)
              bsedt(1,rank)=min(dble(bsedt(1,rank)),t1)
              bsedt(2,rank)=max(dble(bsedt(2,rank)),t1)
              bsedt(3,rank)=bsedt(3,rank)+zt1/(nst1*nst3)
           end do
        end do
     end if

     ! parallel write
     call putbsemat('SCCLI.OUT',sccli,ikkp,iknr,jknr,iq,iqr,nst1,nst3,nst2,nst4)

     deallocate(igqmap,emat12,emat34)
     deallocate(tm,tmi)
     deallocate(scclit)

     ! end loop over (k,kp)-pairs
  end do
  if (rank.eq.0) write(un,'("# ikkp, iknr,ist1,ist3, jknr,ist2,ist4,   &
  & Re(W),            Im(W),             |W|^2,           ang/pi")')
	if (rank.eq.0) close(un)

  call barrier

  ! communicate array-parts wrt. q-points
  call zalltoallv(bsedt,3,procs)
  ! BSE kernel diagonal parameters
  bsedl=minval(dble(bsedt(1,:)))
  bsedu=maxval(dble(bsedt(2,:)))
  bsedd=bsedu-bsedl
  bsed=sum(bsedt(3,:))/nkptnr
  deallocate(bsedt,scclid)
  ! write BSE kernel diagonal parameters
  if (rank.eq.0) call putbsediag('BSEDIAG.OUT')

  call findgntn0_clear

  write(unitout,'("Info(scrcoulint): Screened Coulomb interaction finished")')

end subroutine scrcoulint
