module autormt
   Use modmpi, only: mpiglobal
   Use errors_warnings, only: terminate_if_false
   use precision, only: dp
   use predefined_rgkmax, only: get_predefined_rgkmax
#include "asserts.fpp"
   implicit none 
   private
   public :: get_initial_rmt_rgkmax, &
             scale_rmt_lattice, & 
             optimal_rmt

   contains
      
      !> For a given species, calculates initial muffin-tin radii using the formula
      !> \[ 
      !>    \text{R}_{i} = 1 + 0.25 |\text{Z}_{i}|^{1/3}
      !> \]
      !> where \(\text{Z}_{i}\) is the atomic number of the \(i\)-th species. 
      !> 
      !> REVISION HISTORY:
      !> Created March 2005 (JKD)
      !> Changed the formula, September 2006 (JKD)
      function get_initial_rmt_formula(spzn_is) result(rmt_is)
         !> Atomic number   
         real(dp), intent(in) :: spzn_is 
         !> Scaling factor to calculate initial muffin-tin radii
         real(dp), parameter :: atomic_nr_scale = 0.25_dp
         !> Initial muffin-tin radius
         real(dp) :: rmt_is

         rmt_is =  1.0_dp + atomic_nr_scale * idnint(Abs(spzn_is))**(1.0_dp/3.0_dp)
      
      end function get_initial_rmt_formula

      !> For a given species, gets the initial value for the calculation of the muffin-tin radii.
      !> These initial values are rgkmax values that, for the given species, give a certain predefined 
      !> error in the energy (precision).
      function get_initial_rmt_rgkmax(spzn_is) result(rmt_is)
         !> Atomic number of given species  
         real(dp), intent(in) :: spzn_is 
         !> Initial muffin-tin radius
         real(dp) :: rmt_is
         
         ! Label for atomic number
         character(20) :: z_label, z_real_label
         integer :: z_int
         
         rmt_is = get_predefined_rgkmax(spzn_is)
         
         z_int = idnint(abs(spzn_is))
         write(z_label,'(I0)') z_int
         write(z_real_label,'(F12.6)') abs(spzn_is)

         call terminate_if_false(mpiglobal, rmt_is > 0.0_dp, &
              "Error (autormt): No predefined initial muffin-tin radius exists for nuclear charge Z = " // &
              trim(adjustl(z_real_label)) // " (nearest integer: " // trim(z_label) // "). " // &
              "Predefined values are available only for integer atomic numbers 1–86. " // &
              "For fractional or out-of-range Z, please specify the muffin-tin radius manually." )
      end function get_initial_rmt_rgkmax

      !> Scales the given initial muffin-tin radii in `rmt` using `fixed_rmt`.
      !> The scaling factor is obtained from the given fixed muffin-tin radius `fixed_rmt` and its
      !> given initial muffin-tin value (initial muffin-tin values are obtained from
      !> [[get_initial_rmt_rgkmax(function)]] or [[get_initial_rmt_formula(function)]]).
      !> \[ 
      !>    \text{scaling_factor} = \text{fixed_rmt} / \text{initial_rmt_fixedrmt}
      !> \]
      subroutine scale_rmt_fixedrmt(rmt, fixed_rmt, idx_fixedrmt, nspecies) 
         !> Muffin-tin radii
         real(dp), intent(inout) :: rmt(:)
         !> Value of fixed muffin-tin radius
         real(dp), intent(in) :: fixed_rmt
         !> Index for the species with fixed rmt
         integer, intent(in) :: idx_fixedrmt
         !> Number of species
         integer, intent(in) :: nspecies

         ! Scaling factor 
         real(dp) :: scaling_factor
         ! Initial muffin-tin value of the species with fixed rmt  
         real(dp) :: initial_rmt_fixedrmt

         initial_rmt_fixedrmt = rmt(idx_fixedrmt)

         CALL_ASSERT((fixed_rmt > 0.0_dp),  message='Value for fixed rmt has to be greater than zero.')
         CALL_ASSERT((idx_fixedrmt > 0) .and. (idx_fixedrmt < nspecies),  message='The species index has to be larger than zero and lower or equal to the total number of species.')
         
         scaling_factor = fixed_rmt / initial_rmt_fixedrmt 
         rmt(1:nspecies) = rmt(1:nspecies) * scaling_factor
         rmt(idx_fixedrmt) = fixed_rmt
      end subroutine 

      !> Scales the given initial muffin-tin radii in `rmt`. 
      !> Given the atomic positions and the lattice vectors, a supercell is constructed by translating the unit cell 
      !> and by adding the atomic positions to all translations. 
      !> Next, it computes the distance of the atoms of the unit cell and their neighbours of the translated cells. 
      !> The given initial muffin-tin radii in `rmt` are then all scaled by the calculated minimal distance. 
      !> The `rmt` values can also be further scaled by a global factor, where the maximal value can 
      !> be `autormtscaling` = 1, such that that the closest muffin-tin radii touch.
      subroutine scale_rmt_lattice(rmt, lattice_vect, atomic_positions, autormtscaling, nspecies, natoms)
         !> Number of species
         integer, intent(in) :: nspecies
         !> Number of atoms for each species
         integer, intent(in) :: natoms(nspecies)
         !> Governs distance between muffin-tin spheres. If equal to one, closest muffin-tins touch.
         real(dp), intent(in) :: autormtscaling
         !> Lattice vectors 
         real(dp), intent(in):: lattice_vect(3,3)
         !> Atomic positions 
         real(dp), intent(in) :: atomic_positions(:, :, :) 
         !> Array to save muffin-tin radii
         real(dp), intent(inout) :: rmt(:)

         ! local variables
         integer :: is, js, ia, ja, i1, i2, i3
         real(dp) :: s, v1 (3), v2 (3), t1, t2, t3
            
         ! external function
         real(dp) :: r3dist

         ! determine scaling factor
         s = 1.d10
         do i1 = - 1, 1
            do i2 = - 1, 1
               do i3 = - 1, 1
                  v1 (:) = dble (i1) * lattice_vect(:, &
                  & 1) + dble (i2) * lattice_vect(:, 2) &
                  & + dble (i3) * lattice_vect(:, 3)
                  do is = 1, nspecies
                     do ia = 1, natoms (is)
                        v2 (:) = v1 (:) + atomic_positions(:, ia, is)
                        do js = 1, nspecies
                           t1 = 1.d0 / (rmt(is)+rmt(js))
                           do ja = 1, natoms (js)
                              if ((i1 .Ne. 0) .Or. (i2 .Ne. 0) .Or. (i3 & 
                              & .Ne. 0) .Or. (is .Ne. js) .Or. (ia .Ne. & 
                              & ja)) then                                 
                                 t2 = r3dist (v2, atomic_positions(:, ja, js))
                                 t3 = t1 * t2
                                 if (t3 .Lt. s) then
                                     s = t3 
                                 end if
                              end if
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         s = s * autormtscaling
         ! scale all radii
         do is = 1, nspecies
            ! limit number of decimal digits
            t1 = s * rmt (is) * 10000.0_dp
            t1 = dble (Int(t1)) / 10000.0_dp
            rmt (is) = t1
         end do
 
      end subroutine

      !> Automatically determines the optimal muffin-tin radii. For a given species, either calculates the initial value
      !> for the muffin-tin radius (`init_rad_version` = 0, see function [[get_initial_rmt_formula(function)]]) or uses a predefined rgkmax-value as the 
      !> initial muffin-tin radius (`init_rad_version` = 1, see function [[get_initial_rmt_rgkmax(function)]]).
      !> If `idx_fixedrmt` is present, the initial muffin-tin radii are scaled using the fixed muffin-tin radius
      !> and the initial muffin-tin radius of the species with index `idx_fixedrmt`.
      !> If `idx_fixedrmt` is not present, the muffin-tin radii are scaled in such a way that the closest muffin-tin radii 
      !> in the unit cell touch. The value which then can further govern the distance between the muffin-tins is stored in `autormtscaling`. 
      !> When `autormtscaling` = 1, the closest muffin-tins will continue to touch. If it is smaller the interstitial region enlarges. 
      subroutine optimal_rmt(rmt, spzn, lattice_vect, atomic_positions, autormtscaling, natoms, nspecies, init_rad_version, idx_fixedrmt)
         !> Number of species
         integer, intent(in):: nspecies
         !> Number of atoms for each species
         integer, intent(in) :: natoms(nspecies)
         !> Scaling factor
         real(dp), intent(in) :: autormtscaling
         !> Lattice vectors column wise
         real(dp), intent(in):: lattice_vect(3, 3)
         !> Atomic positions 
         real(dp), intent(in) :: atomic_positions(:, :, :) 
         !> Atomic number for each species  
         real(dp), intent(in) :: spzn(:) 
         !> Input array to save muffin-tin radii
         real(dp), intent(inout) :: rmt(:)
         !> Case for initial muffin-tin radii     
         integer, intent(in) :: init_rad_version
         !> Index for the species with fixed rmt
         integer, optional, intent(in) :: idx_fixedrmt

         integer :: is
         real(dp) :: fixed_rmt 
         
         CALL_ASSERT((init_rad_version == 0) .or. (init_rad_version == 1),  message='Error(autormt): Invalid case for initial muffin-tin radii')

         if (present(idx_fixedrmt)) then 
            fixed_rmt = rmt(idx_fixedrmt)
         end if 

         ! Get initial muffin-tin radii
         select case (init_rad_version) 
            case(0)
               rmt(1:nspecies) = [(get_initial_rmt_formula((spzn(is))), is = 1, nspecies)] 
            case(1)
               rmt(1:nspecies) = [(get_initial_rmt_rgkmax((spzn(is))), is = 1, nspecies)]
         end select

         ! Scale muffin-tin radii
         if (present(idx_fixedrmt)) then 
            call scale_rmt_fixedrmt(rmt, fixed_rmt, idx_fixedrmt, nspecies)
         else 
            call scale_rmt_lattice(rmt, lattice_vect, atomic_positions, autormtscaling, nspecies, natoms)
         end if 

      end subroutine optimal_rmt 


end module autormt

