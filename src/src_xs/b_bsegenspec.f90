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

  integer(4) :: io1, no, nexc, nk, nw
  integer(4) :: iqmt, iqmti, iqmtf, nqmt, iq1, iq2
  logical :: fcoup, foff
  real(8), allocatable :: evals(:), bindevals(:), w(:)
  complex(8), allocatable :: oscir(:)
  complex(8), allocatable :: oscirmat(:,:)
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

    ! Q-point list entries
    !   Use all
    iqmti = 1
    iqmtf = size(input%xs%qpointset%qpoint, 2)
    !   or range
    if(input%xs%bse%iqmtrange(1) /= -1) then 
      iqmti=input%xs%bse%iqmtrange(1)
      iqmtf=input%xs%bse%iqmtrange(2)
    end if
    nqmt = iqmtf-iqmti+1
    iq1 = 1
    iq2 = nqmt

    ! Get BSE structure flags
    fcoup = input%xs%bse%coupling
    ! Use offdiagonal elements
    foff = input%xs%dfoffdiag

    do iqmt = iqmti+iq1-1, iqmti+iq2-1

      ! Read in exciton energies and oscillator strengths

      if(iqmt == 1) then 
        no = 3
      else
        no = 1
      end if

      !! Read in directional components

      do io1 = 1, no

        call readoscillator(iqmt, io1, evals, bindevals, oscir)

        if(.not. allocated(oscirmat)) then 
          allocate(oscirmat(size(oscir), 3))
          oscirmat = zzero
        end if

        oscirmat(:,io1) = oscir
      
      end do
      if(allocated(oscir)) deallocate(oscir)

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
      call makespectrum(iqmt, nexc, nk, evals, oscirmat, symspectr)

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
      if(allocated(evals)) deallocate(evals)

    ! iqmt
    end do

    call barrier(callername=trim(thisname))

  else

    write(*,'("Info(",a,"): Rank ", i4," is waiting...")') trim(thisname), mpiglobal%rank

    call barrier(callername=trim(thisname))

  end if

end subroutine b_bsegenspec
!EOC
