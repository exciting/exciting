! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! REVISION HISTORY:
! Created July 2019 (Ronaldo Rodrigues Pela)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module contains the global variables for the RT-TDDFT implementation
module rttddft_GlobalVariables
  use MD, only: MD_timing
  use precision, only: dp

  implicit none

  private
  ! List of the many global variables can be used publicly
  public :: &
    & apwalm, evecfv_gnd, evecfv_time, evecfv_save, evecsv, &
    & overlap, ham_time, ham_past, ham_predcorr, &
    & pvec, jpara, jparanext, jparaold, jparaspurious, jdia, jind, &
    & pmat, mathcalH, mathcalB, B_time, B_past, pmatmt

  ! Global variables of general purpose
  !> Matching coefficients of the (L)APWs
  complex(dp), allocatable  :: apwalm(:,:,:,:,:)

  !> Basis-expansion coefficients of the groundstate KS-WFs 
  complex(dp),allocatable   :: evecfv_gnd(:,:,:)  
  !> Basis-expansion coefficients of the KS-WFs at time \(t\)
  complex(dp),allocatable   :: evecfv_time(:,:,:) 
  !> Basis-expansion coefficients of the KS-WFs at time \(t\) - auxiliary 
  !> variable used in the predictor-corrector loop
  complex(dp),allocatable   :: evecfv_save(:,:,:) 
  !> Basis-expansion coefficients of the KS-WFs: second-variational coefficients
  complex(dp), allocatable  :: evecsv(:,:,:)

  !> Overlap matrix (of basis functions)
  complex(dp),allocatable   :: overlap(:,:,:)
  !> Hamiltonian matrix at current time \(t\)
  complex(dp),allocatable   :: ham_time(:,:,:)
  !> Hamiltonian matrix at previous time \(t - \Delta t \)
  complex(dp),allocatable   :: ham_past(:,:,:)
  !> Auxiliary hamiltonian matrix, for the predictor-corrector loop
  complex(dp),allocatable   :: ham_predcorr(:,:,:)

  !> Polarization field
  real(dp)                  :: pvec(3)
  !> Paramagnetic component of the current density
  real(dp)                  :: jpara(3)
  !> Auxiliary variable, used to evolve `[[rttddft_GlobalVariables:jpara]]`
  real(dp)                  :: jparanext(3)
  !> Auxiliary variable, used to evolve `[[rttddft_GlobalVariables:jpara]]`
  real(dp)                  :: jparaold(3)
  !> Spurious paramagnetic current density (obtained for \(t=0\) - this should
  !> ideally be zero for a dense `k-grid` mesh)
  real(dp)                  :: jparaspurious(3)
  !> Diamagnetic component of the current density
  real(dp)                  :: jdia(3)
  !> Total (induced) current density
  real(dp)                  :: jind(3)

  !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
  complex(dp), allocatable  :: pmat(:,:,:,:)
  
  !> Muffin-tin part of the Momentum matrix
  complex(dp), allocatable  :: pmatmt(:,:,:,:,:)

  !> `mathcalH` gives the impact of an ion displacement on the hamiltonian matrix
  !> \[ \left[ \left\langle 
  !> \frac{\partial \phi_{\mu'}^{\mathbf{k}}}{\partial \mathbf{R}_J}
  !> \Bigg|\hat{H}\Bigg|\phi_{\mu}^{\mathbf{k}}\right\rangle +
  !> \left\langle\phi_{\mu'}^{\mathbf{k}}\Bigg|\hat{H}\Bigg|\frac{\partial 
  !> \phi_{\mu}^{\mathbf{k}}}{\partial \mathbf{R}_J}\right\rangle \right] 
  !> \]
  complex(dp), allocatable  :: mathcalH(:,:,:,:,:)
  
  !> `mathcalB` measures how the ions displacements affect overlap elements
  !> \[ \mathcal{B}_{J\mu'\mu}^{\mathbf{k}} = \left \langle
  !> \phi_{\mu'}^{\mathbf{k}}\bigg| \frac{\partial}{\partial \mathbf{R}_J}
  !> \phi_{\mu}^{\mathbf{k}} \right\rangle \]
  complex(dp), allocatable  :: mathcalB(:,:,:,:,:)
  
  !> `B_time` quantifies the impact of the Ehrenfest molecular dynamics
  !> on the time evolution of the electonic wavefunctions.
  !> \[B_{\mu'\mu}^{\mathbf{k}} = \left \langle
  !> \phi_{\mu'}^{\mathbf{k}}\left|\frac{d}{d t}\right.
  !> \phi_{\mu}^{\mathbf{k}}\right\rangle =
  !> \sum_J \dot{\mathbf{R}}_J\cdot \mathcal{B}_{J\mu'\mu}^{\mathbf{k}} \]
  complex(dp), allocatable  :: B_time(:,:,:)
  !> Same as `B_time`, but at the previous time step: \(t-\Delta t\)
  complex(dp), allocatable  :: B_past(:,:,:)

end module rttddft_GlobalVariables
