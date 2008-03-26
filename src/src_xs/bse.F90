
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bse
  use modmain
  use modxs
  implicit none
  ! local variables
  character(*), parameter :: thisnam='bse'
  integer, parameter :: iqmt=1
  character(256) :: fname
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iqrnr,isym,jsym,jsymi,igq1,igq2,n,iflg,recl
  integer :: ngridkt(3),iv(3),ivgsym(3),un,unsc,unex,j1,j2,s1,s2,hamsiz
  integer :: ist1,ist2,ist3,ist4,nst12,nst34,nst13,nst24,ikkp,bsetype
  logical :: nosymt,reducekt,tq0,nsc,tphf
  real(8) :: vklofft(3),vqr(3),vq(3),v2(3),s(3,3),si(3,3),t3,abstol
  real(8), allocatable :: potcl(:),rwork(:)
  integer :: igqmap(maxsymcrys),sc(maxsymcrys),ivgsc(3,maxsymcrys)
  integer, allocatable :: iwork(:),ifail(:)
  complex(8), allocatable :: scclit(:,:),sccli(:,:,:,:),ham(:,:),work(:)
  complex(8), allocatable :: scrni(:,:,:),tm(:,:),tmi(:,:)
  complex(8), allocatable :: phf(:,:),emat12(:,:),emat34(:,:)
  logical, allocatable :: done(:)
  integer :: ikkp_,iknr_,jknr_,iq_,iqr_,nst1_,nst2_,nst3_,nst4_
  ! external functions
  integer, external :: iplocnr
  logical, external :: tqgamma
  real(8), external :: dlamch

  ! type of contributions to BSE-Hamlitonian
  ! H = H_diag + 2H_x + H_c
  ! H_diag .......... diagonal term containing IP-energy-differences
  ! H_x ............. exchange term
  ! H_c ............. correlation term
  ! value of bsetype corresponds to
  ! 0........... H = H_diag                     IP-spectrum
  ! 1........... H = H_diag + 2H_x              RPA-spectrum
  ! 2........... H = H_diag + 2H_x + H_c        correlated, spin-singlet
  ! 3........... H = H_diag + H_c               correlated, spin-triplet
  bsetype=0

  !----------------!
  !   initialize   !
  !----------------!
  ! save global variables
  nosymt=nosym
  reducekt=reducek
  ngridkt(:)=ngridk(:)
  vklofft(:)=vkloff(:)
  ! map variables for screened Coulomb interaction
  call initbse
  nosym=nosymscr
  ! no symmetries implemented for screened Coulomb interaction
  reducek=.false.
  ! q-point set of screening corresponds to (k,kp)-pairs
  ngridk(:)=ngridq(:)
  vkloff(:)=vkloffbse(:)
  if (nemptyscr.eq.-1) nemptyscr=nempty
  emattype=1
  call init0
  call init1
  call init2xs
  ! read Fermi energy from file
  call readfermi
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  call ematbdcmbs(emattype)

  ! size of BSE-Hamiltonian
  hamsiz=nkptnr*nstcon*nstval

  write(*,'(a,4i6)') 'nst1,2,3,4',nst1,nst2,nst3,nst4
  allocate(sccli(nst1,nst2,nst1,nst2))
  allocate(excli(nst1,nst2,nst1,nst2))
  ! allocate BSE-Hamiltonian (large matrix, up to several GB)
  allocate(ham(hamsiz,hamsiz))
  ham(:,:)=zzero

  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,nst4, &
       sccli(:,:,:,:)
  call genfilname(basename='SCCLI',dotext='.OUT',filnam=fnamesc)
  call genfilname(basename='EXCLI',dotext='.OUT',filnam=fnameex)
  open(unsc,file=trim(fnamesc),form='unformatted',action='read', &
       status='old',access='direct',recl=recl)
  open(unex,file=trim(fnameex),form='unformatted',action='read', &
       status='old',access='direct',recl=recl)

  ! read in energies
  do iknr=1,nkptnr
     call getevalsv(vkl(1,iknr),evalsv(1,iknr))
  end do
  write(*,*) 'reading energies done'


  ! set up BSE-Hamiltonian
  do iknr=1,nkptnr
     do jknr=iknr,nkptnr
        ikkp=ikkp+1
        iv(:)=ivknr(:,jknr)-ivknr(:,iknr)
        iv(:)=modulo(iv(:),ngridk(:))
        ! q-point (reduced)
        iqr=iqmapr(iv(1),iv(2),iv(3))
        ! q-point (non-reduced)
        iq=iqmap(iv(1),iv(2),iv(3))

        if (bsetype.ne.0) then
           ! read screened Coulomb interaction
           read(unsc,rec=ikkp) ikkp_,iknr_,jknr_,iq_,iqr_,nst1_,nst2_,nst3_, &
                nst4_,sccli
           if ((ikkp.ne.ikkp_).or.(iknr.ne.iknr_).or.(jknr.ne.jknr_).or. &
                (iq.ne.iq_).or.(iqr.ne.iqr_).or.(nst1.ne.nst1_).or. &
                (nst2.ne.nst2_).or.(nst3.ne.nst3_).or.(nst4.ne.nst4_)) then
              write(*,*)
              write(*,'("Error(kernxc_bse): wrong indices for screened Coulomb&
                   & interaction")')
              write(*,'(" indices (ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,&
                   &nst4)")')
              write(*,'(" current:",i6,3x,2i4,2x,2i4,2x,4i4)') ikkp,iknr,jknr,&
                   iq,iqr,nst1,nst2,nst3,nst4
              write(*,'(" file   :",i6,3x,2i4,2x,2i4,2x,4i4)') ikkp_,iknr_,&
                   jknr_,iq_,iqr_,nst1_,nst2_,nst3_,nst4_
              write(*,*)
              call terminate
           end if

           ! read exchange Coulomb interaction
           read(unex,rec=ikkp) ikkp_,iknr_,jknr_,iq_,iqr_,nst1_,nst2_,nst3_, &
                nst4_,excli
           if ((ikkp.ne.ikkp_).or.(iknr.ne.iknr_).or.(jknr.ne.jknr_).or. &
                (iq.ne.iq_).or.(iqr.ne.iqr_).or.(nst1.ne.nst1_).or. &
                (nst2.ne.nst2_).or.(nst3.ne.nst3_).or.(nst4.ne.nst4_)) then
              write(*,*)
              write(*,'("Error(kernxc_bse): wrong indices for screened Coulomb&
                   & interaction")')
              write(*,'(" indices (ikkp,iknr,jknr,iq,iqr,nst1,nst2,nst3,&
                   &nst4)")')
              write(*,'(" current:",i6,3x,2i4,2x,2i4,2x,4i4)') ikkp,iknr,jknr,&
                   iq,iqr,nst1,nst2,nst3,nst4
              write(*,'(" file   :",i6,3x,2i4,2x,2i4,2x,4i4)') ikkp_,iknr_,&
                   jknr_,iq_,iqr_,nst1_,nst2_,nst3_,nst4_
              write(*,*)
              call terminate
           end if
        end if

        do ist1=1,nst1
           do ist2=1,nst2
              do ist3=1,nst1
                 do ist4=1,nst2
                    s1=hamidx(ist1,ist2,iknr)
                    s2=hamidx(ist3,ist4,jknr)
                    ! add diagonal term
                    if (ist1.eq.ist3.and.ist2.eq.ist4.and.iknr.eq.jknr) then
                       ham(s1,s2)=ham(s1,s2)+ &
                            evalsv(ist2+istocc,iknr)-evalsv(ist1,iknr)+scissor
                    end if
                    ! add exchange term
                    if ((bsetype.eq.1).or.(bsetype.eq.2)) then
                       ham(s1,s2)=ham(s1,s2)+ &
                            2.d0*excli(ist1,ist2,ist3,ist4)
                    end if
                    ! add correlation term
                    if ((bsetype.eq.2).or.(bsetype.eq.3)) then
                       ham(s1,s2)=ham(s1,s2)+ &
                            sccli(ist1,ist2,ist3,ist4)
                    end if
                 end do
              end do
           end do
        end do


        ! end loop over (k,kp)-pairs
     end do
  end do

  deallocate(excli,sccli)
  lwork=2*hamsiz
  allocate(work(lwork),rwork(7*hamsiz),iwork(5*hamsiz),ifail(hamsiz))
  
  
  ! assign lower triangle of Hamiltonian
  do s1=1,hamsiz
     do s2=1,s1-1
        ham(s1,s2)=conjg(ham(s2,s1))
     end do
  end do


  ! *** solve eigenvalue problem
  abstol=2.d0*dlamch('S')
  ! allocate help arrays

  

!!$ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
!!$                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK,
!!$                   IWORK, IFAIL, INFO )
!!$*     .. Scalar Arguments ..
!!$      CHARACTER          JOBZ, RANGE, UPLO
!!$      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
!!$      DOUBLE PRECISION   ABSTOL, VL, VU
!!$*     ..
!!$*     .. Array Arguments ..
!!$      INTEGER            IFAIL( * ), IWORK( * )
!!$      DOUBLE PRECISION   RWORK( * ), W( * )
!!$      COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )

  ! deallocate Hamiltonian array
  deallocate(ham,work,rwork,iwork,ifail)

  ! read momentum matrix elements


  ! calculate spectrum



!///////////////////////////////////////////////////////////////////////////////
contains

  integer function hamidx(iv,ic,ik)
    use modxs
    implicit none
    hamidx=iv + nstval*(ic-1) + nstval*nstcon*(ik-1)
  end function hamidx



end subroutine bse

