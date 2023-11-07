module wigner_seitz_trial_energy_generator
   use modmpi, only: terminate
   use precision, only: dp

   implicit none
   
   private
   public :: generate_wigner_seitz_trial_energies

contains

!> Generates the trial energy for a given local orbital from
!> its number of nodes by using the Wigner-Seitz rules.
   function generate_wigner_seitz_trial_energies(l, principal_n, spr, nr, vr, wf_tolerance, energy_tolerance) result(e_trial)

      !> angular momentum
      integer, intent(in) :: l
      ! prinicipal quantum number
      integer, intent(in) :: principal_n
      !> number of muffin-tin radial points for the species
      integer, intent(in) :: nr
      !> species radial mesh
      real(dp), intent(in) :: spr(nr)
      !> radial component of the muffin-tin effective potential
      real(dp), intent(in) :: vr(nr)
      !> tolerance for considering the slope of a wave function as zero
      real(dp), intent(in) :: wf_tolerance
      !> target accuracy for the trial energy in the bisection 
      real(dp), intent(in) :: energy_tolerance

      !> computed trial energy
      real(dp) :: e_trial

      ! local variables
      !> number of nodes
      integer :: nodes
      ! number of nodes
      integer :: nn
      ! major component of the radial wavefunction
      real(dp) :: p0(nr)
      ! radial derivative of p0
      real(dp) :: p1(nr)
      ! minor component of the radial wavefunction
      real(dp) :: q0(nr)
      ! radial derivative of q0
      real(dp) :: q1(nr)
      ! major component of the radial wavefunction that gets zero at muffin-tin boundary
      real(dp) :: p0_zero_at_mt(nr)
      ! energy for which the wavefunction is zero at muffin-tin boundary
      real(dp) :: energy_for_zero_wf_at_mt
      ! lower bound energy
      real(dp) :: e_lower_bound
      ! upper bound energy
      real(dp) :: e_upper_bound
      ! radial derivative of the wavefunction for lower bound energy at muffin-tin boundary
      real(dp) :: p1_lower_bound
      ! radial derivative of the wavefunction for upper bound energy at muffin-tin boundary
      real(dp) :: p1_upper_bound
      ! mean of lower and upper bound energy
      real(dp) :: e_mean
      ! radial derivative of the wavefunction for mean energy at muffin-tin boundary
      real(dp) :: p1_mean
      ! flag to pick the equation (Dirac or Schroedinger) used in rdirac
      Logical  :: dirac_eq
      ! flag to pick a quick-and-dirty algorithm for integrating the Dirac equation in rdirac
      Logical  :: sloppy


      ! Error message
      character(1024) :: message

      if (principal_n < 1) Then
         call terminate("Error(wigner_seitz_trial_energy_generator): principal quantum number < 1")
      end if

! Schroedinger equation is used in rdirac
      dirac_eq = .false.
      sloppy = .false.

! Compute number of nodes from principal quantum number
      nodes = principal_n - l - 1

! Compute energy for which the wave function becomes 0 at the muffin-tin boundary
      energy_for_zero_wf_at_mt = 0._dp
      Call rdirac(0, principal_n, l, l + 1, nr, spr, vr, energy_for_zero_wf_at_mt, p0_zero_at_mt, q0, dirac_eq, sloppy)

! Choose the lower bound energy.
      e_lower_bound = 0._dp
      if (nodes == 0) then
         ! If the number of nodes is already equal to zero, we can't further reduce the nodes.
         ! In this case the lower bound energy is computed by taking the difference of two times the
         ! energies with zero nodes and one node.
         Call rdirac(0, principal_n + 1, l, l + 1, nr, spr, vr, e_lower_bound, p0, q0, dirac_eq, sloppy)
         e_lower_bound = 2*energy_for_zero_wf_at_mt - e_lower_bound
      else
         ! Compute energy for which the wave function with one node less becomes 0 at the muffin-tin boundary
         Call rdirac(0, principal_n - 1, l, l + 1, nr, spr, vr, e_lower_bound, p0, q0, dirac_eq, sloppy)
      end if

      e_upper_bound = energy_for_zero_wf_at_mt

      ! If the slope of the wavefunction that vanishes at the muffin-tin boundary vanishes as well,
      ! we already found the trial energy.
      if (abs(p0_zero_at_mt(nr) - p0_zero_at_mt(nr - 1)) < wf_tolerance) then
         e_trial = e_upper_bound
         return
      end if

      ! Else, we search for the energy via bisection.
      ! Compute radial derivative of the wave function of upper bound energy. The upper bound energy is
      ! the one for which the wavefunction vanishes at the muffin-tin boundary.
      Call rschroddme(0, l, 0, e_upper_bound, nr, spr, vr, nn, p0, p1, q0, q1)
      p1_upper_bound = p1(nr)

      ! Compute the radial derivatice of the wavefunction of the lower_bound energy.
      Call rschroddme(0, l, 0, e_lower_bound, nr, spr, vr, nn, p0, p1, q0, q1)
      p1_lower_bound = p1(nr)

      if (e_upper_bound < e_lower_bound) then
         call terminate("Error(gentrialenergy): e_upper_bound < e_lower_bound")
      end if

      
      do while (e_upper_bound - e_lower_bound > energy_tolerance)
         e_mean = 0.5_dp*(e_upper_bound + e_lower_bound)
         ! Compute the radial derivative of the wavefunction with energy e_mean
         Call rschroddme(0, l, 0, e_mean, nr, spr, vr, nn, p0, p1, q0, q1)
         p1_mean = p1(nr)

         ! Perform bisection to find energy for which the radial derivative
         ! vanishes at the muffin-tin boundary.
         if (p1_mean*p1_upper_bound < 0) then
            p1_lower_bound = p1_mean
            e_lower_bound = e_mean
         else
            p1_upper_bound = p1_mean
            e_upper_bound = e_mean
         end if
      end do

      ! Compute trial energy by taking the mean of the energy for which the
      ! wavefunction vanishes at muffin-tin boundary and the energy for which the
      ! radial derivative vanishes at the muffin-tin boundary.
      e_mean = 0.5_dp * (e_upper_bound + e_lower_bound)
      e_trial = 0.5_dp*(energy_for_zero_wf_at_mt + e_mean)

      ! checks if the radial derivative of the wavefunction with e_mean diverges
      Call rschroddme(0, l, 0, e_mean, nr, spr, vr, nn, p0, p1, q0, q1)
      if ((p1(nr)) > (maxval(abs(p1(1:nr - 1))))) then
         write (*, '("Error(wigner_seitz_trial_energy_generator): ")')
         write (*, '("Wavefunction derivative diverges for l =", I2, " and nodes =", I2, ".")') l, nodes
         write (message, '("Remove corresponding local orbital or choose trial energy by hand.")')
         call terminate(message)
      end if

      ! checks if the wavefunction with e_trial diverges
      Call rschroddme(0, l, 0, e_trial, nr, spr, vr, nn, p0, p1, q0, q1)
      if ((p0(nr)) > (maxval(abs(p0(1:nr - 1))))) then
         write (*, '("Error(wigner_seitz_trial_energy_generator): ")')
         write (*, '("Wavefunction diverges for l =", I2, " and nodes =", I2, ".")') l, nodes
         write (message, '("Remove corresponding local orbital or choose trial energy by hand.")')
         call terminate(message)
      end if

   end function generate_wigner_seitz_trial_energies

end module wigner_seitz_trial_energy_generator
