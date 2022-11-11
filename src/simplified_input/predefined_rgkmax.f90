module predefined_rgkmax
    use precision, only: dp
    use modmpi, only: mpiglobal
    use errors_warnings, only: terminate_if_false

    implicit none 
    private
    public :: get_predefined_rgkmax

    contains
        
    !> Returns a predefined rgkmax value for a given atomic number.
    !> These rgkmax values have been obtained by studying the convergence of total energies for elemental solids as reported in
    !> Carbogno, Christian, et al. "Numerical quality control for DFT-based materials databases." npj Computational Materials 8.1 (2022): 1-8. 
    !> Applying these rgkmax values to elemental solids should yield a precision of about 0.1 meV/atom in total energy
    !> (Assuming all other numerical parameters allow this kind of precision).
    function get_predefined_rgkmax(spzn) result(rgkmax)
        !> Atomic number   
        real(dp), intent(in) :: spzn
        !> Returned rgkmax value 
        real(dp) :: rgkmax

        ! Label for atomic number
        character(10) :: z_label

        !> Number of total rgkmax values 
        integer, parameter :: n_initial_rgkmax = 86 
        !> Rgkmax values for atomic numbers going from 1 to 86 
        real(dp), parameter :: inital_rgkmax_params(n_initial_rgkmax) = &
            [5.83543_dp, & ! H, 1
            8.226318_dp, & ! He, 2
            8.450962_dp, & ! Li, 3
            8.307929_dp, & ! Be, 4
            8.965808_dp, & ! B, 5
            9.376204_dp, & ! C, 6
            9.553568_dp, & ! N, 7
            10.239864_dp, & ! O, 8
            10.790975_dp, & ! F, 9
            10.444355_dp, & ! Ne, 10
            10.636286_dp, & ! Na, 11
            10.579793_dp, & ! Mg, 12
            10.214125_dp, & ! Al, 13
            10.605334_dp, & ! Si, 14
            10.356352_dp, & ! P, 15
            9.932381_dp, & ! S, 16
            10.218153_dp, & ! Cl, 17
            10.466519_dp, & ! Ar, 18
            10.877475_dp, & ! K, 19
            10.774763_dp, & ! Ca, 20
            11.580691_dp, & ! Sc, 21
            11.800971_dp, & ! Ti, 22
            11.919804_dp, & ! V, 23
            12.261896_dp, & ! Cr, 24
            12.424606_dp, & ! Mn, 25
            12.571031_dp, & ! Fe, 26
            12.693836_dp, & ! Co, 27
            12.781331_dp, & ! Ni, 28
            12.619806_dp, & ! Cu, 29
            12.749802_dp, & ! Zn, 30
            12.68135_dp, & ! Ga, 31
            12.802838_dp, & ! Ge, 32
            12.78568_dp, & ! As, 33
            12.898916_dp, & ! Se, 34
            12.4_dp, & ! Br, 35
            10.596757_dp, & ! Kr, 36
            11.34606_dp, & ! Rb, 37
            10.857573_dp, & ! Sr, 38
            11.324413_dp, & ! Y, 39
            11.6642_dp, & ! Zr, 40
            11.859519_dp, & ! Nb, 41
            11.892673_dp, & ! Mo, 42
            12.30847_dp, & ! Tc, 43
            12.551024_dp, & ! Ru, 44
            12.740728_dp, & ! Rh, 45
            12.879424_dp, & ! Pd, 46
            13.02709_dp, & ! Ag, 47
            13.080576_dp, & ! Cd, 48
            13.230621_dp, & ! In, 49
            13.450665_dp, & ! Sn, 50
            13.495632_dp, & ! Sb, 51
            13.261039_dp, & ! Te, 52
            13.432654_dp, & ! I, 53
            11.329591_dp, & ! Xe, 54
            13.343047_dp, & ! Cs, 55
            13.011835_dp, & ! Ba, 56
            -1.0_dp, & ! La, 57
            -1.0_dp, & ! Ce, 58
            -1.0_dp, & ! Pr, 59
            -1.0_dp, & ! Nd, 60 
            -1.0_dp, & ! Pm, 61
            -1.0_dp, & ! Sm, 62
            -1.0_dp, & ! Eu, 63
            -1.0_dp, & ! Gd, 64
            -1.0_dp, & ! Tb, 65
            -1.0_dp, & ! Dy, 66
            -1.0_dp, & ! Ho, 67
            -1.0_dp, & ! Er, 68
            -1.0_dp, & ! Tm, 69
            -1.0_dp, & ! Yb, 70
            15.134859_dp, & ! Lu, 71
            14.955721_dp, & ! Hf, 72
            14.607311_dp, & ! Ta, 73
            13.930505_dp, & ! W, 74
            13.645267_dp, & ! Re, 75
            13.629439_dp, & ! Os, 76
            13.450805_dp, & ! Ir, 77
            13.069046_dp, & ! Pt, 78
            13.226699_dp, & ! Au, 79
            13.261342_dp, & ! Hg, 80
            13.365992_dp, & ! Tl, 81
            13.557571_dp, & ! Pb, 82
            13.565048_dp, & ! Bi, 83
            13.579543_dp, & ! Po, 84
            -1.0_dp, & ! At, 85
            12.273924_dp] ! Rn, 86
        
        write(z_label,'(I3)') idnint(Abs(spzn))
        call terminate_if_false(mpiglobal, ((idnint(Abs(spzn)) .gt. 0) .and. (idnint(Abs(spzn)) .le.  n_initial_rgkmax)),&
                                & "Error: Out of bounds (predefined_rgkmax): for given atomic number "// trim(z_label) // " no &
                                & predefined rgkmax value exists.")
        
        rgkmax = inital_rgkmax_params(idnint(Abs(spzn)))

    end function get_predefined_rgkmax

end module predefined_rgkmax