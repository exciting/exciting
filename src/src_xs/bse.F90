
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bse
! !INTERFACE:
subroutine bse
! !USES:
  use modmain
  use modxs
  use m_genwgrid
  use m_getpmat
  use m_genfilname
  use m_getunit
  use m_writeeps
! !DESCRIPTION:
!   Solves the Bethe-Salpeter equation (BSE). The BSE is treated as equivalent
!   effective eigenvalue problem (thanks to the spectral theorem that can
!   be applied to the original BSE in the case of a statically screened Coulomb
!   interaction). The effective BSE-Hamiltonian consists of three parts
!   originating from different sources. It reads
!   $$ H^{\rm eff} = H^{\rm diag} + 2H^{\rm x} + H^{\rm c}, $$
!   where $H^{\rm diag}$ is the diagonal part stemming from the independent
!   particle transitions, $H^{\rm x}$ denotes the exchange-term caused by the
!   unscreened (bare) Coulomb interaction, whereas $H^{\rm c}$ accounts for the
!   particle-hole correlations and is originating from the screened Coulomb
!   interaction.
!   For the purpose of describing independent particle transitions with the
!   BSE only the diagonal term is referred to:
!   $$ H^{\rm eff} = H^{\rm diag}. $$
!   By neglecting the correlation part in the effective Hamiltonian we arrive
!   at the {\it random phase approximation} (RPA)
!   $$ H^{\rm eff} = H^{\rm diag} + 2H^{\rm x}. $$
!   Investigations on the spin-structure of the BSE-Hamiltonian show that there
!   are tow channels, namely the {\it singlet}-channel as solution to the
!   Hamiltonian
!   $$  H^{\rm eff} = H^{\rm diag} + 2H^{\rm x} + H^{\rm c} $$
!   and a {\it triplet} channel with the exchange-part being absent.
!   $$ H^{\rm eff} = H^{\rm diag} + H^{\rm c}. $$
!   The equation of the eigenvalue problem is given by
!   $$ \sum_{v'c'{\bf k'}} H^{\rm eff}_{vc{\bf k},v'c'{\bf k'}}
!       A^{\lambda}_{v'c'{\bf k'}}
!       =  \varepsilon_{\lambda} A^{\lambda}_{vc{\bf k}}. $$
!   For the diagonalization of the Hamiltonian, a LAPACK-routine ({\tt zheevx})
!   is invoked to obtain the eigenvalues $\varepsilon_{\lambda}$ and
!   eigenvectors $A^{\lambda}_{vc{\bf k}}$ (alternatively, a time-evolution
!   method shall be implemented to obtain the macroscopic dielectric function
!   directly).
!   Consequently, the transition amplitudes $t_{\lambda}$ are calculated
!   according to
!   $$ t^{i}_{\lambda} = \left|\sum_{vc{\bf k}} A^{\lambda}_{vc{\bf k}} 
!      \frac{ p^{i}_{vc{\bf k}} }{ \varepsilon_{c{\bf k}}- 
!                                  \varepsilon_{v{\bf k}} } \right|^2. $$
!   Here, the index $i$ labels the polarization and the matrix elements
!   $p^{i}_{vc{\bf k}}$ are the ones for the $i$-th component of the momentum
!   operator in Cartesian coordinates.
!   The macroscopic dielectric function (MDF) is obtained by the realation
!   $$ {\rm Im}\; \epsilon^{i}_{\rm M}(\omega) = \frac{8\pi^2}{V} 
!                     \sum_{\lambda} t^{i}_{\lambda}
!                     \delta(\omega-\varepsilon_{\lambda}+\Delta),$$
!   where $\epsilon^{i}_{\rm M}$ is the MDF for the $i$-th polarization, $V$
!   denotes the crystal volume and $\Delta$ is a constant shift of the
!   conduction bands (scissors shift). The delta-function in the latter
!   expression is convoluted with a (symmetrized) Lorentzian 
!   $$ \pi\delta(\omega-\omega_0) = \lim_{\eta\rightarrow 0} \left[
!                         \frac{\eta}{(\omega-\omega_0)^2+\eta^2} +
!                         \frac{\eta}{(-\omega-\omega_0)^2-\eta^2} \right] =
!     \pi\delta(\omega-\omega_0) +  \pi\delta(\omega+\omega_0)       $$
!   which is true for $\omega \ge 0$ if $\omega_0>0$. In doing so, the analytic
!   property ${\rm Im}\epsilon_{\rm M}(0)=0$ is fulfilled.
!   The broadening $\eta$ in the latter expression is adjusted by the
!   {\tt broad} parameter. (All parts of the documentation written by
!   S. Sagmeister are part of the author's PhD-thesis.)
!
! !REVISION HISTORY:
!   Created June 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  integer, parameter :: iqmt=0
  real(8), parameter :: epsortho=1.d-12
  integer :: iknr,jknr,iqr,iq,iw,iv2(3),s1,s2,hamsiz,nexc,ne
  integer :: ist1,ist2,ist3,ist4,ikkp,oct,iv,ic,nvdif,ncdif
  real(8) :: de,egap,ts0,ts1
  ! allocatable arrays
  integer, allocatable :: sor(:)
  real(8), allocatable :: beval(:),w(:),oszsa(:)
  complex(8), allocatable :: excli(:,:,:,:),sccli(:,:,:,:),ham(:,:)
  complex(8), allocatable :: bevec(:,:),pm(:,:,:),pmat(:),oszs(:),spectr(:)
  ! external functions
  integer, external :: l2int
  ! *** TODO: symmetrize head of DM for spectrum
  !---------------------------!
  !     exciton variables     !   USE this ****************************
  !---------------------------!
  !!if (allocated(excite)) deallocate(excite)
  !!allocate(excite(nexcitmax,3))
  !!excite(:,:)=0.d0
  !!if (allocated(excito)) deallocate(excito)
  !!allocate(excito(nexcitmax,3))
  !!excito(:,:)=0.d0
  call init0
  call init1
  call init2
  call xssave0
  ! read Fermi energy from file
  call readfermi
  ! initialize states below and above the Fermi energy
  call initocc(nbfbse,nafbse)
  ! use eigenvector files from screening-calculation
  call genfilname(dotext='_SCR.OUT',setfilext=.true.)
  call findocclims(iqmt,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  nvdif=nstocc0-nbfbse
  ncdif=nstunocc0-nafbse

write(*,*) 'nvdif,ncdif',nvdif,ncdif

  ! ****************************************************
  emattype=2
  call ematbdcmbs(emattype)
  write(unitout,*)
  write(unitout,'("Info(bse): information on number of states:")')
  write(unitout,'(" number of states below Fermi energy in Hamiltonian:",i6)') &
       nbfbse
  write(unitout,'(" number of states above Fermi energy in Hamiltonian:",i6)') &
       nafbse
  write(unitout,'(" ranges of states according to BSE matrix:")')
  write(unitout,'("  range of first index and number  :",2i6,3x,i6)') &
       istl1,istu1,nst1
  write(unitout,'("  range of second index and number :",2i6,3x,i6)') &
       istl2,istu2,nst2
  write(unitout,'("  range of third index and number  :",2i6,3x,i6)') &
       istl3,istu3,nst3
  write(unitout,'("  range of fourth index and number :",2i6,3x,i6)') &
       istl4,istu4,nst4
  if ((nvdif.lt.0).or.(ncdif.lt.0)) then
     write(unitout,*)
     write(unitout,'("Error(bse): inconsistency in ranges of states - check &
          & routine for nvdif,ncdif")')
     write(unitout,*)
     call terminate
  end if
  ! size of BSE-Hamiltonian
  hamsiz=nbfbse*nafbse*nkptnr
  ! allocate arrays for Coulomb interactons
  allocate(sccli(nst1,nst3,nst2,nst4))
  allocate(excli(nst1,nst3,nst2,nst4))
  ! allocate BSE-Hamiltonian (large matrix, up to several GB)
  allocate(ham(hamsiz,hamsiz))
  ham(:,:)=zzero
  ! read in energies
  do iknr=1,nkptnr
     call getevalsv(vkl(1,iknr),evalsv(1,iknr))
  end do
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
        egap=1.d8
        ! set up matrix
        do ist1=1+nvdif,nst1
           do ist3=1,nst3-ncdif
              do ist2=1+nvdif,nst2
                 do ist4=1,nst4-ncdif
                    s1=hamidx(ist1-nvdif,ist3,iknr,nbfbse,nafbse)
                    s2=hamidx(ist2-nvdif,ist4,jknr,nbfbse,nafbse)
                    ! add diagonal term
                    if (s1.eq.s2) then
                       de=evalsv(ist3+istocc,iknr)-evalsv(ist1,iknr)+scissor
                       ham(s1,s2)=ham(s1,s2)+de
                       egap=min(egap,de)
                    end if
                    ! add exchange term
                    select case(trim(bsetype))
                    case('rpa','singlet')
                       ham(s1,s2)=ham(s1,s2)+2.d0*excli(ist1,ist3,ist2,ist4)
                    end select
                    ! add correlation term
                    select case(trim(bsetype))
                    case('singlet','triplet')
!!!                       ham(s1,s2)=ham(s1,s2)-sccli(ist2,ist4,ist1,ist3)
if (s1.eq.s2)                       ham(s1,s2)=ham(s1,s2)-sccli(ist1,ist3,ist2,ist4)
                    end select
                 end do
              end do
           end do
        end do
        ! end loop over (k,kp)-pairs
     end do
  end do
  deallocate(excli,sccli)
  write(unitout,*)
  write(unitout,'("Info(bse): invoking Lapack routine ZHEEVX")')
  write(unitout,'(" size of BSE-Hamiltonian       : ",i8)') hamsiz
  write(unitout,'(" number of requested solutions : ",i8)') nexcitmax
  ! allocate eigenvector and eigenvalue arrays
  allocate(beval(hamsiz),bevec(hamsiz,hamsiz))
  ! set number of excitons
  ne=hamsiz
  call timesec(ts0)
  ! diagonalize Hamiltonian
  call bsesoldiag(hamsiz,ne,ham,beval,bevec)
  call timesec(ts1)
  ! deallocate BSE-Hamiltonian
  deallocate(ham)
  write(unitout,'(" timing (in seconds)           :",f12.3)') ts1-ts0
  ! number of excitons to consider
  nexc=hamsiz
  allocate(oszs(nexc),oszsa(nexc),sor(nexc),pmat(hamsiz))
  allocate(w(nwdos),spectr(nwdos))
  call genwgrid(nwdf,wdos,acont,0.d0,w_real=w)
  do oct=1,noptcomp
     oszs(:)=zzero
     call genfilname(basename='EPSILON',tq0=.true.,oc1=oct,oc2=oct, &
          bsetype=bsetype,scrtype=screentype,filnam=fneps)
     ! read momentum matrix elements
     allocate(pm(3,nstsv,nstsv))
     do iknr=1,nkptnr
        call getpmat(iknr,vkl,1,nstsv,1,nstsv,.true.,'PMAT_XS.OUT',pm)
        do ist1=1+nvdif,nstsv-nstunocc0
           do ist2=nstocc0+1,nstsv-ncdif
              s1=hamidx(ist1-nvdif,ist2-nstocc0,iknr,nbfbse,nafbse)
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
                1.d0/(w(iw)-beval(s1)+zi*broad) + &
                1.d0/(-w(iw)-beval(s1)-zi*broad) )
        end do
     end do
     spectr(:)=l2int(oct.eq.oct)*1.d0-spectr(:)*8.d0*pi/omega/nkptnr
     ! write BSE spectrum
     call writeeps(iqmt,oct,oct,w,spectr,fneps)
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
contains

  integer function hamidx(i1,i2,ik,n1,n2)
    implicit none
    integer, intent(in) :: i1,i2,ik,n1,n2
    hamidx=i2+n2*(i1-1)+n1*n2*(ik-1)
  end function hamidx

end subroutine bse
!EOC
