module autormt
   Use modmpi, only: mpiglobal
   Use errors_warnings, only: terminate_if_false
   use precision, only: dp
   use predefined_rgkmax, only: get_predefined_rgkmax

   implicit none 
   private
   public :: get_inital_rmt_rgkmax, &
             scale_rmt, & 
             optimal_rmt

   contains
      
      !> For a given species, calculates initial muffin-tin radii using the formula
      !> \[ 
      !>    \text{R}_{i} = 1 + 0.25 |\text{Z}_{i}|^{1/3}
      !> \]
      !> where \(Z_i\) is the atomic number of the \(i\)-th species. 
      !> 
      !> REVISION HISTORY:
      !> Created March 2005 (JKD)
      !> Changed the formula, September 2006 (JKD)
      function get_inital_rmt_formula(spzn_is) result(rmt_is)
         !> Atomic number   
         real(dp), intent(in) :: spzn_is 
         !> Scaling factor to calculate initial muffin-tin radii
         real(dp), parameter :: atomic_nr_scale = 0.25_dp
         !> Initial muffin-tin radius
         real(dp) :: rmt_is

         rmt_is =  1.0_dp + atomic_nr_scale * idnint(Abs(spzn_is))**(1.0_dp/3.0_dp)
      
      end function get_inital_rmt_formula

      !> For a given species, gets the initial value for the calculation of the muffin-tin radii.
      !> These initial values are rgkmax values that, for the given species, give a certain predefined 
      !> error in the energy (precision).
      function get_inital_rmt_rgkmax(spzn_is) result(rmt_is)
         !> Atomic number of given species  
         real(dp), intent(in) :: spzn_is 
         !> Initial muffin-tin radius
         real(dp) :: rmt_is
         
         ! Label for atomic number
         character(10) :: z_label
         
         rmt_is = get_predefined_rgkmax(spzn_is)
         write(z_label,'(I3)') idnint(Abs(spzn_is))
         call terminate_if_false(mpiglobal, rmt_is /= -1.0_dp ,"(Error (autormt): initial muffin-tin& 
                                 & radius for given atomic number"// trim(z_label) // " does not exist.")
      end function get_inital_rmt_rgkmax

      !> Scales the given initial muffin-tin radii in \(\text{rmt}\). 
      !> Given the atomic positions and the lattice vectors, a supercell is constructed by translating the unit cell 
      !> and by adding the atomic positions to all translations. 
      !> Next, it computes the distance of the atoms of the unit cell and their neighbours of the translated cells. 
      !> The given initial muffin-tin radii in \(\text{rmt}\) are then all scaled by the calculated minimal distance. 
      !> The \(\text{rmt}\) values can also be further scaled by a global factor, where the maximal value can 
      !> be \(\text{scale_global_rmt} = 1 \), such that that the clostest muffin-tin radii touch.
      subroutine scale_rmt(rmt, lattice_vect, atomic_positions, scale_global_rmt, nspecies, natoms)
         !> Number of species
         integer, intent(in) :: nspecies
         !> Number of atoms for each species
         integer, intent(in) :: natoms(nspecies)
         !> Governs distance between muffin-tin spheres. If equal to one, closest muffin-tins touch.
         real(dp), intent(in) :: scale_global_rmt
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
         s = s * scale_global_rmt
         ! scale all radii
         do is = 1, nspecies
            ! limit number of decimal digits
            t1 = s * rmt (is) * 10000.0_dp
            t1 = dble (Int(t1)) / 10000.0_dp
            rmt (is) = t1
         end do
 
      end subroutine

      !> Automatically determines the optimal muffin-tin radii. For a given species, calculates the initial value
      !> for the muffin-tin radius and then scales them, in such a way that the 
      !> clostest muffin-tin radii in the unit cell touch.  
      !> The value which then can further govern the distance between the muffin-tins is stored in 
      !> scale_global_rmt. When \(\text{scale_global_rmt} = 1\), the closest muffin-tins will continue to touch. If it is smaller
      !> the interstitial region enlarges. 
      subroutine optimal_rmt(rmt, spzn, lattice_vect, atomic_positions, scale_global_rmt, natoms, nspecies, init_rad_version)
         !> Number of species
         integer, intent(in):: nspecies
         !> Number of atoms for each species
         integer, intent(in) :: natoms(nspecies)
         !> Scaling factor
         real(dp), intent(in) :: scale_global_rmt
         !> Lattice vectors column wise
         real(dp), intent(in):: lattice_vect(3,3)
         !> Atomic positions 
         real(dp), intent(in) :: atomic_positions(:, :, :) 
         !> Atomic number for each species  
         real(dp), intent(in) :: spzn(:) 
         !> Input array to save muffin-tin radii
         real(dp), intent(inout) :: rmt(:)
         !> Case for initial muffin-tin radii     
         integer, intent(in) :: init_rad_version

         integer :: is

         ! Get initial muffin-tin radii
         select case (init_rad_version) 
            case(0)
               rmt(1:nspecies) = [(get_inital_rmt_formula((spzn(is))), is = 1, nspecies)]   
            case(1)
               rmt(1:nspecies) = [(get_inital_rmt_rgkmax((spzn(is))), is = 1, nspecies)]
         end select
     
         ! Scale muffin-tin radii 
         call scale_rmt(rmt, lattice_vect, atomic_positions, scale_global_rmt, nspecies, natoms)
      
      end subroutine optimal_rmt 

end module autormt

