
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bse
  use modmain
  use modxs
  use m_genwgrid
  use m_getpmat
  use m_genfilname
  use m_getunit
  use m_writeeps
  implicit none
  ! local variables
  character(*), parameter :: thisnam='bse'
  integer, parameter :: iqmt=1
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iw,ngridkt(3),iv2(3),s1,s2,hamsiz,nexc, ne
  integer :: ist1,ist2,ist3,ist4,ikkp,oct,iv,ic,nvdif,ncdif
  logical :: nosymt,reducekt
  real(8) :: vklofft(3),de,egap,cpu0,cpu1
  ! allocatable arrays
  real(8), allocatable :: beval(:),w(:),oszsa(:)
  integer, allocatable :: sor(:)
  complex(8), allocatable :: excli(:,:,:,:),sccli(:,:,:,:),ham(:,:)
  complex(8), allocatable :: bevec(:,:),pm(:,:,:),pmat(:),oszs(:),spectr(:)
  ! external functions
  integer, external :: l2int

  !TODO: symmetrize head of DM for spectrum

  ! type of contributions to BSE-Hamlitonian
  ! H = H_diag + 2H_x + H_c
  ! H_diag .......... diagonal term containing IP-energy-differences
  ! H_x ............. exchange term
  ! H_c ............. correlation term
  ! value of bsetype corresponds to
  ! IP ................ H = H_diag                     IP-spectrum
  ! RPA ............... H = H_diag + 2H_x              RPA-spectrum
  ! BSE (singlet) ..... H = H_diag + 2H_x + H_c        correlated, spin-singlet
  ! BSE (triplet) ..... H = H_diag + H_c               correlated, spin-triplet

  ! reset file extension to default
  call genfilname(setfilext=.true.)
  egap=1.d8

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
  call initoccbse(nbfbse,nafbse)
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  call ematbdcmbs(emattype)
  ! file names
  fnpmat='PMAT_XS.OUT'

  nvdif=nstval-nbfbse
  ncdif=nstcon-nafbse


  write(*,'("number of states below Fermi energy:",i6)') nbfbse
  write(*,'("number of states above Fermi energy:",i6)') nafbse
  write(*,*) 'nvdif',nvdif
  write(*,*) 'ncdif',ncdif
  write(*,'(a,4i6)') 'nst1,2,3,4',nst1,nst2,nst3,nst4

  if ((nvdif.lt.0).or.(ncdif.lt.0)) stop 'bse: bad nbfbse,nafbse'

  ! size of BSE-Hamiltonian
  hamsiz=nbfbse*nafbse*nkptnr

  allocate(sccli(nst1,nst2,nst1,nst2))
  allocate(excli(nst1,nst2,nst1,nst2))
  ! allocate BSE-Hamiltonian (large matrix, up to several GB)
  allocate(ham(hamsiz,hamsiz))
  ham(:,:)=zzero

  ! read in energies
  do iknr=1,nkptnr
     call getevalsv(vkl(1,iknr),evalsv(1,iknr))
  end do
  write(*,*) 'reading energies done'

  ! set up BSE-Hamiltonian
  ikkp=0
  do iknr=1,nkptnr
     do jknr=iknr,nkptnr
        ikkp=ikkp+1
        iv2(:)=ivknr(:,jknr)-ivknr(:,iknr)
        iv2(:)=modulo(iv2(:),ngridk(:))
        ! q-point (reduced)
        iqr=iqmapr(iv2(1),iv2(2),iv2(3))
        ! q-point (non-reduced)
        iq=iqmap(iv2(1),iv2(2),iv2(3))
        select case(trim(bsetype))
        case('singlet','triplet')
           ! read screened Coulomb interaction
           call getbsemat('SCCLI.OUT',ikkp,nst1,nst3,sccli)
        end select
        ! read exchange Coulomb interaction
        select case(trim(bsetype))
        case('rpa','singlet')
           call getbsemat('EXCLI.OUT',ikkp,nst1,nst3,excli)
        end select
        ! set up matrix
        do ist1=1+nvdif,nst1
           do ist2=1,nst2-ncdif
              do ist3=1+nvdif,nst1
                 do ist4=1,nst2-ncdif
                    s1=hamidx(ist1-nvdif,ist2,iknr,nbfbse,nafbse)
                    s2=hamidx(ist3-nvdif,ist4,jknr,nbfbse,nafbse)
                    ! add diagonal term
                    if (s1.eq.s2) then
                       de=evalsv(ist2+istocc,iknr)-evalsv(ist1,iknr)+scissor
                       ham(s1,s2)=ham(s1,s2)+de
                       egap=min(egap,de)
                    end if
                    ! add exchange term
                    select case(trim(bsetype))
                    case('rpa','singlet')
                       ham(s1,s2)=ham(s1,s2)+ &
                            2.d0*excli(ist1,ist2,ist3,ist4)
                    end select
                    ! add correlation term
                    select case(trim(bsetype))
                    case('singlet','triplet')
                       ham(s1,s2)=ham(s1,s2)-               &
                            sccli(ist1,ist2,ist3,ist4)
                    end select
                 end do
              end do
           end do
        end do
        ! end loop over (k,kp)-pairs
     end do
  end do
  deallocate(excli,sccli)

  write(*,*) 'call to zheevx..............'
  write(*,*) 'Hamiltonian size is: ', hamsiz

  ! diagonalize Hamlitonian
  allocate(beval(hamsiz),bevec(hamsiz,hamsiz))
  ne=hamsiz
  call cpu_time(cpu0)
  call bsesoldiag(hamsiz,ne,ham,beval,bevec)
  call cpu_time(cpu1)
  deallocate(ham)

  write(*,*) 'zheevx finished..............'
  write(*,'("Number of requested solutions:",i8)') ne
  write(*,'("Time (in seconds) in zheevx:",f12.3)') cpu1-cpu0

  ! number of excitons to consider
  nexc=hamsiz
  allocate(oszs(nexc),oszsa(nexc),sor(nexc),pmat(hamsiz))
  allocate(w(nwdos),spectr(nwdos))
  call genwgrid(nwdf,wdos,acont,0.d0,w_real=w)

  do oct=1,noptcomp
     oszs(:)=zzero
     call genfilname(basename='EPSILON_BSE',tq0=.true.,oc1=oct,oc2=oct, &
          bsetype=bsetype,scrtype=screentype,filnam=fneps)

     ! read momentum matrix elements
     allocate(pm(3,nstsv,nstsv))
     do iknr=1,nkptnr
        call getpmat(iknr,vkl,.true.,trim(fnpmat),pm)
        do ist1=1+nvdif,nstsv-nstcon
           do ist2=nstval+1,nstsv-ncdif
              s1=hamidx(ist1-nvdif,ist2-nstval,iknr,nbfbse,nafbse)
              pmat(s1)=pm(oct,ist1,ist2)
           end do
        end do
     end do
     deallocate(pm)

     ! calculate oscillators for spectrum  
     do s1=1,nexc
        do iknr=1,nkptnr
           do iv=1,nbfbse
              do ic=1,nafbse
                 s2=hamidx(iv,ic,iknr,nbfbse,nafbse)
                 oszs(s1)=oszs(s1)+bevec(s2,s1)*pmat(s2)/ &
                      (evalsv(ic+istocc,iknr)-evalsv(iv+nvdif,iknr))
              end do
           end do
        end do
     end do

     spectr(:)=zzero
     do iw=1,nwdos
        do s1=1,nexc
           ! Lorentzian lineshape
           spectr(iw)=spectr(iw) + abs(oszs(s1))**2 * ( &
                1.d0/(w(iw)-beval(s1)+zi*brdtd) + &
                1.d0/(-w(iw)-beval(s1)-zi*brdtd) )
        end do
     end do
     spectr(:)=l2int(oct.eq.oct)*1.d0-spectr(:)*8.d0*pi/omega/nkptnr

     ! write BSE spectrum
     call writeeps(0,w,spectr,fneps)

     ! oscillator strengths
     do s2=1,hamsiz
        write(983,'(i8,5g18.10)') s2,beval(s2)*escale,(beval(s2)-egap)*escale, &
             abs(oszs(s2))
     end do

     ! oscillator strengths sorted
     oszsa=abs(oszs)
     call sortidx(hamsiz,oszsa,sor)
     sor=sor(hamsiz:1:-1)
     do s1=1,hamsiz
        s2=sor(s1)
        write(984,'(i8,4g18.10)') s1,beval(sor(s2))*escale, &
             (beval(sor(s2))-egap)*escale,abs(oszs(s2))
     end do


     ! end loop over optical components
  end do



  !/////////////////////////////////////////////////////////////////////////////
contains

  integer function hamidx(iv,ic,ik,nv,nc)
    use modxs
    implicit none
    integer, intent(in) :: iv,ic,ik,nv,nc
    hamidx=iv + nv*(ic-1) + nv*nc*(ik-1)
  end function hamidx



end subroutine bse

!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////


subroutine bsesoldiag(hamsiz,nev,ham,eval,evec)
  implicit none
  ! arguments
  integer, intent(in) :: hamsiz,nev
  complex(8), intent(in) :: ham(hamsiz,hamsiz)
  real(8), intent(out) :: eval(nev)
  complex(8), intent(out) :: evec(hamsiz,nev)
  ! local variables
  real(8) :: vl,vu,abstol
  integer :: il,iu,neval,lwork,info
  ! allocatable arrays
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:),ifail(:)
  ! external functions
  real(8), external :: dlamch
  ! smallest and largest eigenvalue indices
  il=1
  iu=nev
  ! tolerance parameter
  abstol=2.d0*dlamch('S')
  ! workspace size (*** improve later ***)
  lwork=(32+1)*hamsiz
  allocate(work(lwork),rwork(7*hamsiz),iwork(5*hamsiz),ifail(hamsiz))
  ! LAPACK 3.0 call
  call zheevx('V','I','U',hamsiz,ham,hamsiz,vl,vu,il,iu,abstol,neval,eval, &
       evec,hamsiz,work,lwork,rwork,iwork,ifail,info)
  if (info.ne.0) then
     write(*,*)
     write(*,'("Error(bse): zheevx returned non-zero info:",i6)') info
     write(*,*)
     call terminate
  end if
  ! deallocate Hamiltonian array
  deallocate(work,rwork,iwork,ifail)
end subroutine bsesoldiag
