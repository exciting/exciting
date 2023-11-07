module trial_energy_selection

   use constants, only: y00
   use wigner_seitz_trial_energy_generator, only: generate_wigner_seitz_trial_energies
   use precision, only: dp
   use mod_APW_LO, only: default_lorbn, default_apwn

   implicit none

   private
   public :: select_local_orbital_trial_energies, select_apw_trial_energies

contains


   !> If a principal quantum number other than the default value is specified in the species file, this subroutine will 
   !> automatically calculate the trial energy for the corresponding local orbital radial function using the Wigner-Seitz rules.
   !> On the other hand, if the principal quantum number is set to its default value, the subroutine will not modify the trial energy. 
   subroutine select_local_orbital_trial_energies(nlorb, lorbord, lorbl, lorbn, nspecies, idxas, nrmt, spr, veffmt, lorbe0)

      !> number of local orbitals
      integer, intent(in) :: nlorb(:)
      !> local orbital order
      integer, intent(in) :: lorbord(:, :)
      !> local-orbital angular momentum
      integer, intent(in) :: lorbl(:, :)
      !> local-orbital principal quantum number
      integer, intent(in) :: lorbn(:, :, :)
      !> number of species
      integer, intent(in) :: nspecies
      !> index to atoms and species
      integer, intent(in) :: idxas(:, :)
      !> species radial mesh
      real(dp), intent(in) :: spr(:, :)
      !> number of muffin-tin radial points for each species
      integer, intent(in) :: nrmt(:)
      !> muffin-tin effective potential
      real(dp), intent(in) :: veffmt(:, :)
      !> local orbital energy
      real(dp), intent(inout) :: lorbe0(:, :, :)
      
      ! local variables
      Integer :: is, ilo, io
      Real(dp) :: v(maxval(nrmt))

      !> tolerance for considering the slope of a wave function as zero
      real(dp), parameter :: wf_slope_tolerance = 1.e-4_dp 
      !> tolerance for the trial energy in the bisection 
      real(dp), parameter :: energy_tolerance = 1.e-6_dp


      do is = 1, nspecies
         do ilo = 1, nlorb(is)
            do io = 1, lorbord(ilo, is) ! Looping over all local orbitals.
               ! If the number of nodes of a local orbital is given, the trial energy is computed automatically.
               if (lorbn(io, ilo, is) /= default_lorbn) then
                  v = veffmt(1:nrmt(is), idxas(1, is))*y00
                  ! Computing trial energy with Wigner-Seitz algorithm.
                  lorbe0(io, ilo, is) = generate_wigner_seitz_trial_energies(lorbl(ilo, is), lorbn(io, ilo, is), &
                                                                             spr(:, is), nrmt(is), v, wf_slope_tolerance, energy_tolerance)
               end if
            end do
         end do
      end do

   end subroutine select_local_orbital_trial_energies

   
   !> If a principal quantum number other than the default value is specified in the species file, this subroutine will 
   !> automatically calculate the trial energy for the corresponding (L)APW radial function using the Wigner-Seitz rules.
   !> On the other hand, if the principal quantum number is set to its default value, the subroutine will not modify the trial energy.
   subroutine select_apw_trial_energies(maxapword, maxlapw, apwn, nspecies, idxas, nrmt, spr, veffmt, apwe0)

      !> maximum APW order
      integer, intent(in) :: maxapword
      !> maximum APW angular momentum
      integer, intent(in) :: maxlapw
      !> APW principal quantum number
      integer, intent(in) :: apwn(:, 0:, :)
      !> number of species
      integer, intent(in) :: nspecies
      !> index to atoms and species
      integer, intent(in) :: idxas(:, :)
      !> species radial mesh
      real(dp), intent(in) :: spr(:, :)
      !> number of muffin-tin radial points for each species
      integer, intent(in) :: nrmt(:)
      !> muffin-tin effective potential
      real(dp), intent(in) :: veffmt(:, :)

      !> APW energy
      real(dp), intent(inout) :: apwe0(:, 0:, :)

      ! local variables
      Integer :: is, l, io, ilx, nlx
      Real(dp) :: v(maxval(nrmt))

      !> tolerance for considering the slope of a wave function as zero        
      real(dp), parameter :: wf_slope_tolerance = 1.e-4_dp 
      !> tolerance for the trial energy in the bisection                                                          
      real(dp), parameter :: energy_tolerance = 1.e-6_dp


      do is = 1, nspecies
         Do l = 0, maxlapw
            do io = 1, maxapword
               ! If the number of nodes of a basis function is given, the trial energy is computed automatically.
               if (apwn(io, l, is) /= default_apwn) then
                  v = veffmt(1:nrmt(is), idxas(1, is))*y00
                  ! Computing trial energy with Wigner-Seitz algorithm.
                  apwe0(io, l, is) = generate_wigner_seitz_trial_energies(l, apwn(io, l, is), spr(:, is), nrmt(is), &
                                                                           v, wf_slope_tolerance, energy_tolerance)
               end if
            end do
         end do
      end do

   end subroutine select_apw_trial_energies

end module trial_energy_selection
