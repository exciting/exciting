!> Short-lived module to house refactored initialisation routines
!> for init0.f90
module tmp_mod_init0
    use precision, only: dp
    use modinput, only: input_type

    implicit none

    private
    public :: map_atoms_per_species_to_atomic_index


contains
   !> Map atoms per species to an atomic index over all atoms in the system
   subroutine map_atoms_per_species_to_atomic_index (nspecies, natoms, idxas, natmmax, natmtot)

      use constants, only: maxatoms, maxspecies

      !> Number of species
      integer, intent(in) :: nspecies
      !> Number of atoms for each species
      integer, intent(in) :: natoms (:)
      !> Maximum number of atoms over all the species
      integer, intent(out) :: natmmax
      !> Total number of atoms
      integer, intent(out):: natmtot
      !> Map atoms per species to an atomic index over all atoms in the system
      integer, intent(out) :: idxas (maxatoms, maxspecies)

      !> Local variables
      integer :: ias, is, ia

      natmmax = 0
      ias = 0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = ias + 1
            idxas (ia, is) = ias
         End Do
         natmmax = Max (natmmax, natoms(is))
      End Do
      natmtot = ias

   end subroutine

end module 
