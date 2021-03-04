! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! REVISION HISTORY:
! Created July 2019 (Ronaldo Rodrigues Pela)
! Reference: https://arxiv.org/abs/2102.02630

!> This module contains the global variables for the RT-TDDFT implementation
module rttddft_GlobalVariables
  use precision, only: dp

  implicit none

  private
  !> List of the many global variables can be used
  public :: fileavec, filejind, filetime, filepvec, filenexc, fileetot, &
    & filepmat, fileinfortddft, formatTime, fileevalevec, &
    & nsteps, time, tstep, tend, &
    & maxstepsPredictorCorrector, predictorCorrector, &
    & printTimesGeneral, printTimesDetailed, calculateTotalEnergy, calculateNexc, &
    & method, tolPredCorr, &
    & apwalm, evecfv_gnd, evecfv_time, evecfv_save, evecsv, &
    & overlap, ham_time, ham_past, ham_predcorr, &
    & aext, aind, atot, pvec, jpara, jparanext, jparaold, jparaspurious, &
    & jdia, jind, &
    & nkicks, dirkick, t0kick, wkick, amplkick, &
    & ntrapcos, dirtrapcos, ampltrapcos, omegatrapcos, phasetrapcos, &
    & t0trapcos, trtrapcos, wtrapcos, &
    & nsinsq, dirsinsq, amplsinsq, omegasinsq, phasesinsq, &
    & t0sinsq, tpulsesinsq, &
    & pmat, &
    & timesecRTTDDFT, TimingRTTDDFT

  !> This type stores the time (in seconds) spent in the procedures of RT-TDDFT
  type :: TimingRTTDDFT
    real(dp) :: t_wvf
    real(dp) :: t_dens, t_dens_rho, t_dens_symrf, t_dens_rfmtctof
    real(dp) :: t_dens_addrhocr, t_dens_charge, t_dens_rhonorm
    real(dp) :: t_uppot, t_poteff, t_genveffig, t_genmeffig
    real(dp) :: t_curr, t_obtaina
    real(dp) :: t_hmlint, t_ham, t_upham
    real(dp) :: t_predcorr
    real(dp) :: t_toten, t_nexc, t_screenshot
    real(dp) :: t_iteration
  end type TimingRTTDDFT

  !> Files
  integer                   :: fileavec, filejind, filetime, filepvec
  integer                   :: filenexc, fileetot, filepmat
  integer                   :: fileinfortddft
  character(50),parameter   :: formatTime = '(A30,F12.6)'
  character(50),parameter   :: fileevalevec = 'EVALEVEC_GS'

  !> Variables related to time \( t \)
  integer                   :: nsteps
  real(dp)                  :: time, tstep, tend

  !> Some general (useful) variables
  integer                   :: maxstepsPredictorCorrector
  logical                   :: predictorCorrector
  logical                   :: printTimesGeneral, printTimesDetailed
  logical                   :: calculateTotalEnergy, calculateNexc
  character(6)              :: method
  real(dp)                  :: tolPredCorr


  !> Matching coefficients (apwalm)
  complex(dp), allocatable  :: apwalm(:,:,:,:,:)

  !> Coefficients: KS WFs expanded in terms of basis functions
  complex(dp),allocatable   :: evecfv_gnd(:,:,:)  !> groundstate KS WFs
  complex(dp),allocatable   :: evecfv_time(:,:,:) !> current KS WFs
  complex(dp),allocatable   :: evecfv_save(:,:,:) !> auxiliary variable (predcorr)
  !> Second-variational coefficients
  complex(dp), allocatable  :: evecsv(:,:,:)

  !> Overlap of basis functions
  complex(dp),allocatable   :: overlap(:,:,:)
  !> Hamiltonian at current time
  complex(dp),allocatable   :: ham_time(:,:,:)
  !> Hamiltonian at previous time (t-Deltat)
  complex(dp),allocatable   :: ham_past(:,:,:)
  !> Auxiliary hamiltonian, to help with convergence of predictor-corrector
  complex(dp),allocatable   :: ham_predcorr(:,:,:)

  !> Fields
  real(dp)                  :: aext(3), aind(3), atot(3)
  real(dp)                  :: pvec(3)
  real(dp)                  :: jpara(3), jparanext(3), jparaold(3), jparaspurious(3)
  real(dp)                  :: jdia(3), jind(3)

  !> (Parameters of) the delta-Kicks
  integer                   :: nkicks
  character,allocatable     :: dirkick(:)
  real(dp), allocatable     :: t0kick(:), wkick(:), amplkick(:)
  !> (Parameters of) Cossine pulses modulated by a trapezoidal function
  integer                   :: ntrapcos
  character,allocatable     :: dirtrapcos(:)
  real(dp), allocatable     :: ampltrapcos(:), omegatrapcos(:), phasetrapcos(:)
  real(dp), allocatable     :: t0trapcos(:), trtrapcos(:), wtrapcos(:)
  !> (Parameters of) Cossine pulses modulated by sin squared
  integer                   :: nsinsq
  character,allocatable     :: dirsinsq(:)
  real(dp), allocatable     :: amplsinsq(:), omegasinsq(:), phasesinsq(:)
  real(dp), allocatable     :: t0sinsq(:), tpulsesinsq(:)

  !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
  complex(dp), allocatable  :: pmat(:,:,:,:)

contains

  !> Checks the clock (current execution time in seconds)
  !> and stores the difference between this time and timei (passed as parameter).
  !> This helps to evaluate how many seconds the execution of a subroutine is taking.
  !> In the end, we make assign the current time to timei
  !> In the code that calls timesecRTTDDFT we usually have e.g.:
  !> call timesec(timei)
  !> call subroutine
  !> call timesecRTTDDFT(timei, timef, timediff )
  !> (code to treat/store timediff)
  !> Now, timediff contains the time elapsed to execute the subroutine
  !> The code could follow as:
  !> call second_subroutine
  !> call timesecRTTDDFT(timei, timef, timediff )
  !> (code to treat/store timediff)
  !> call third_subroutine
  !> call timesecRTTDDFT(timei, timef, timediff )
  !> (code to treat/store timediff)
  subroutine timesecRTTDDFT(timei, timef, timediff )
    real(dp),intent(inout) :: timei
    real(dp),intent(out)   :: timef, timediff
    call timesec(timef)
    timediff = timef - timei
    timei = timef
  end subroutine timesecRTTDDFT

end module rttddft_GlobalVariables
