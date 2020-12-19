
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Module modmain
      !TODO(Alex) Issue #21. Restrict the use of globals from modmain
      use constants, only: pi, twopi, fourpi
      Use mod_atoms
      Use mod_lattice
      Use mod_muffin_tin
      Use mod_spin
      Use mod_Gvector
      Use mod_symmetry
      Use mod_kpoint
      Use mod_SHT
      Use mod_qpoint
      Use mod_Gkvector
      Use mod_potential_and_density
      Use mod_charge_and_moment
      Use mod_APW_LO
      Use mod_eigensystem
      Use mod_eigenvalue_occupancy
      Use mod_corestate
      Use mod_energy
      Use mod_force
      Use mod_plotting
      Use mod_DOS_optics_response
      Use mod_LDA_LU
      Use mod_RDMFT
      Use mod_misc
      Use mod_timing
      Use constants
      Use mod_phonon
      Use mod_OEP_HF
      Use mod_convergence
      Use mod_names
      Logical :: lwarning = .False.
End Module
!
