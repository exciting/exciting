
module init_autormt_test
    use precision, only: dp
    use modmpi, only: mpiinfo
    use unit_test_framework, only : unit_test_type
    use math_utils, only: all_close, round_down
    use autormt, only: get_inital_rmt_rgkmax, &
                            scale_rmt, & 
                            optimal_rmt

                    
    implicit none
    private
    public :: init_autormt_test_driver
  
contains 

    subroutine init_autormt_test_driver(mpiglobal, kill_on_failure)
        !> mpi environment
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program upon failure of an assertion
        logical, intent(in), optional :: kill_on_failure
  
        !> Test object
        type(unit_test_type) :: test_report
        
        !> Number of assertions
        integer, parameter :: n_assertions = 6
  
        call test_report%init(n_assertions, mpiglobal)
  
        ! Run unit tests
        call test_initrmt(test_report)
        call test_rmtscaling(test_report)

        if (present(kill_on_failure)) then
            call test_report%report('autormt', kill_on_failure)
         else
            call test_report%report('autormt')
         end if
   
        call test_report%finalise()
   
    end subroutine init_autormt_test_driver

    !> Tests whether correct initial rgkmax values are returned. 
    subroutine test_initrmt(test_report)
        !> Test object
        type(unit_test_type), intent(inout) :: test_report
          
        !> Number of species
        integer :: nspecies = 7
        !> Muffin-tin radii vector
        real(dp) :: rmt(7)
        !> Atomic Numbers: can be negative in exciting
        real(dp) :: spzn(7) = [2, -4, 8, -84, 44, 56, -80] 

        integer :: is

        !> Reference 
        real(dp) :: ref(7) = [8.226318_dp, 8.307929_dp, 10.239864_dp, 13.579543_dp, 12.551024_dp, 13.011835_dp, 13.261342_dp]
        
        do is = 1, nspecies
            rmt(is) = get_inital_rmt_rgkmax(spzn(is))
        end do

        call test_report%assert(all_close(rmt, ref), & 
                                'Tests that correct initial muffin-tin radii for given atomic number are returned.')
    end subroutine test_initrmt
    
    !> Tests whether correctly scaled muffin-tin radii are returned.
    subroutine test_rmtscaling(test_report)
        !> Test object
        type(unit_test_type), intent(inout) :: test_report
        
        !> Case for initial muffin-tin radii: empirical version using rgkmax values 
        integer :: init_rad_version = 1
        !> Number of species 
        integer :: nspecies = 2
        !> Atomic numbers for each species 
        real(dp) :: spzn(2) 
        !> Number of atoms per species
        integer:: natoms(2) = [1, 1]
        !> Muffin-tin radii arrays
        real(dp) :: rmt1(2), rmt2(2), rmt3(2), rmt4(2), rmt5(2)
        !> Scaling factor
        real(dp) :: scale_global_rmt = 1.0_dp
        !> Atomic positions
        real(dp) :: atomic_pos1(3, 1, 2), atomic_pos2(3, 1, 2), atomic_pos3(3, 1, 2)
        !> Lattice vector
        real(dp) :: lattice_vect(3, 3) = reshape([0.5_dp, 0.0_dp, 0.5_dp, &
                                                  0.5_dp, 0.5_dp, 0.0_dp, &
                                                  0.0_dp, 0.5_dp, 0.5_dp], [3,3])
        !> Bond length
        real(dp) :: bond_length
        !> References
        real(dp) :: reference1(2), reference2(2), reference3(2), reference4(2), reference5(2)

        ! Case 1: Test on equal initial rmt values
        spzn = [3, 3]
        atomic_pos1 = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                               0.25_dp, 0.25_dp, 0.25_dp], [3,1,2])
        bond_length = sqrt(3.0_dp)*0.25_dp
        reference1 = [round_down(0.5_dp * bond_length, 4), round_down(0.5_dp * bond_length, 4)]
        call optimal_rmt(rmt1, spzn, lattice_vect, atomic_pos1, scale_global_rmt, natoms, nspecies, init_rad_version)
        call test_report%assert(all_close(rmt1, reference1), & 
                                'Tests that correct muffin-tin radius is returned. &
                                Expected: sum of final rmt values equal the bond length and therefore each rmt value is half the bond length &
                                because the species are the same.') 
        
        ! Case 2: Test on different initial rmt values
        spzn = [1, 2] ! This gives initial rgkmax values: [5.835430, 8.226318]
        reference2 = [round_down(bond_length / (8.226318_dp / 5.835430_dp + 1), 4), round_down(bond_length / (5.835430_dp / 8.226318_dp + 1), 4)]
        call optimal_rmt(rmt2, spzn, lattice_vect, atomic_pos1, scale_global_rmt, natoms, nspecies, init_rad_version)
        call test_report%assert(all_close(rmt2, reference2), & 
                                'Tests that correct muffin-tin radius is returned. &
                                Expected: sum of final rmt values equal the bond length and with a ratio determined by the given initial rgkmax values.')

        ! Case 3: Test for equal atomic positions
        atomic_pos2 = reshape([0.5_dp, 0.5_dp, 0.5_dp, &
                                0.5_dp, 0.5_dp, 0.5_dp], [3,1,2]) 
        reference3 = [0.0_dp, 0.0_dp]
        call optimal_rmt(rmt3, spzn, lattice_vect, atomic_pos2, scale_global_rmt, natoms, nspecies, init_rad_version)
        call test_report%assert(all_close(rmt3, reference3), & 
                                'Tests that scaled muffin-tin radii are zero for equal atomic positions. &
                                Expected: all rmt values are zero.')
        
        ! Case 4: Test on equal initial rmt values
        spzn = [3, 3]
        atomic_pos3 = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                               0.5_dp, -0.5_dp, 0.5_dp], [3,1,2])
        reference4 = [0.25_dp, 0.25_dp]
        call optimal_rmt(rmt4, spzn, lattice_vect, atomic_pos3, scale_global_rmt, natoms, nspecies, init_rad_version)
        call test_report%assert(all_close(rmt4, reference4), & 
                                'Tests that correct muffin-tin radius is returned. &
                                Expected: rmt values are a quarter of the length of given cube.') 

        ! Case 5: Test on different initial rmt values
        rmt5 = [1.0_dp, 1.5_dp] 
        bond_length = 0.5 
        reference5 = [bond_length / (1.5_dp / 1.0_dp + 1), bond_length / (1.0_dp / 1.5_dp + 1)]
        call scale_rmt(rmt5, lattice_vect, atomic_pos3, scale_global_rmt, nspecies, natoms)
        call test_report%assert(all_close(rmt5, reference5), & 
                                'Tests that correct muffin-tin radius is returned. &
                                Expected: sum of final rmt values equal the bond length and they have a ratio determined by &
                                the given initial rmt values.')

    end subroutine test_rmtscaling
     
end module init_autormt_test


