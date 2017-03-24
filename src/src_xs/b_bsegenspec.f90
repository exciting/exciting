!BOP
! !ROUTINE: b_bsegenspec
! !INTERFACE:
subroutine b_bsegenspec()
! !USES:
  use modmpi
  use modinput
  use m_readoscillator
  use modbse, only: nk_bse
  use m_makespectrum
  use mod_constants, only: zzero, h2ev
  use mod_kpoint, only: nkptnr
  use mod_lattice, only: omega
  use m_genwgrid
! !DESCRIPTION:
!   Collects BSE results from files EXCITON*.OUT and computes
!   spectra anew.
!
! !REVISION HISTORY:
!   Created March, 2017, BA
!EOP
!BOC

  implicit none

  integer(4) :: iqmt, io1, no, nexc, nk, nw
  logical :: fcoup, fti, foff
  real(8), allocatable :: evals(:), bindevals(:), evalsim(:), w(:)
  complex(8), allocatable :: oscir(:), oscia(:)
  complex(8), allocatable :: oscirmat(:,:), osciamat(:,:)
  complex(8), allocatable, dimension(:,:,:) :: symspectr

  character(*), parameter :: thisname = "b_bsegenspec"

  if(mpiglobal%rank == 0) then 

    ! General init
    call init0
    ! k-grid init
    call init1
    ! q-point and qmt-point setup
    !   Init 2 sets up (task 446):
    !   * A list of momentum transfer vectors form the q-point list 
    !     (modxs::vqlmt and mod_qpoint::vql)
    !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
    !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
    !   * G+qmt quantities (modxs)
    !   * The square root of the Coulomb potential for the G+qmt points
    !   * Reads STATE.OUT
    !   * Generates radial functions (mod_APW_LO)
    call init2

    ! Get Q point index
    iqmt = input%xs%bse%iqmt
    if(iqmt == -1) then 
      write(*,'("Warning(",a,"): iqmt = -1, setting it to 1")') trim(thisname)
    end if

    ! Get BSE structure flags
    fcoup = input%xs%bse%coupling
    fti = input%xs%bse%ti
    ! Use offdiagonal elements
    foff = input%xs%dfoffdiag

    ! Read in exciton energies and oscillator strengths

    if(iqmt == 1) then 
      no = 3
    else
      no = 1
    end if

    !! Read in directional components

    do io1 = 1, no

      if(fcoup .and. .not. fti) then 
        call readoscillator(iqmt, io1, evals, bindevals, oscir, oscia, evalsim)
      else
        call readoscillator(iqmt, io1, evals, bindevals, oscir)
      end if

      if(.not. allocated(oscirmat)) then 
        allocate(oscirmat(size(oscir), 3))
        oscirmat = zzero
        if(fcoup .and. .not. fti) then 
          allocate(osciamat(size(oscir), 3))
          osciamat = zzero
        end if
      end if

      oscirmat(:,io1) = oscir
      if(fcoup .and. .not. fti) then 
        osciamat(:,io1) = oscia
      end if
    
    end do
    if(allocated(oscir)) deallocate(oscir)
    if(allocated(oscia)) deallocate(oscia)

    !! Make the spectrum

    nexc = size(evals)
    nk = nkptnr
    
    nk_bse = nkptnr

    write(*,*) "nexc=", nexc
    write(*,*) "nk", nk
    write(*,*) "omega", omega

    if(input%xs%tevout) then 
      evals = evals/h2ev
    end if

    ! Calculate lattice symmetrized spectrum.
    if(fcoup) then 
      if(fti) then
        call makespectrum_ti(iqmt, nexc, nk, evals, oscirmat, symspectr)
      else
        call makespectrum_full(iqmt, nexc, nk, evals, oscirmat, osciamat, symspectr)
      end if
    else
      call makespectrum_tda(iqmt, nexc, nk, evals, oscirmat, symspectr)
    end if

    ! Generate an evenly spaced frequency grid 
    nw = input%xs%energywindow%points
    allocate(w(nw))
    call genwgrid(nw, input%xs%energywindow%intv,&
      & input%xs%tddft%acont, 0.d0, w_real=w)
    ! Generate and write derived optical quantities
    call writederived(iqmt, symspectr, nw, w)
    deallocate(w)

    ! Some cleaning up
    deallocate(symspectr)
    if(allocated(oscirmat)) deallocate(oscirmat)
    if(allocated(osciamat)) deallocate(osciamat)
    if(allocated(evals)) deallocate(evals)
    if(allocated(evalsim)) deallocate(evalsim)

    call barrier

  else

    write(*,'("Info(",a,"): Rank ", i4," is waiting...")') trim(thisname), mpiglobal%rank

    call barrier

  end if

end subroutine b_bsegenspec
!EOC
