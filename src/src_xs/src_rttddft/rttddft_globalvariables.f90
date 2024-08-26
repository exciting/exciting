! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! REVISION HISTORY:
! Created July 2019 (Ronaldo Rodrigues Pela)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module contains the global variables for the RT-TDDFT implementation
module rttddft_GlobalVariables
  use precision, only: dp

  implicit none

  private
  ! List of the many global variables can be used publicly
  public :: mathcalH, mathcalB, B_time, B_past

  ! Global variables of general purpose

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
