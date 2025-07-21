module mod_gen_lo

    use precision, only: dp, i32
    use modmpi, only: terminate
    use mod_atoms, only: spr, idxas, natmtot, natoms, nspecies
    use mod_APW_LO, only: lorbl, wfkappa, lorbe, lorbdm, lorbord, nlomax, maxlorbord, nlorb
    use mod_potential_and_density, only: veffmt
    use mod_muffin_tin, only: nrmtmax, nrmt, rmt
    use constants, only: y00
    use mgga_poteff, only: veffmt_gga
    use mod_convergence, only: iscl
    use modinput, only: input

    implicit none
 
    private
    public :: genlofr, get_normalized_radial_functions_and_matching_coefficients
 
contains

!> Generates the local-orbital radial functions and matching coefficients. This is done by integrating                                                            
!> the scalar relativistic Schr\"{o}dinger equation or Dirac equation                                                                   
!> (or its energy deriatives) at the current linearisation energies using the                                                           
!> spherical part of the effective potential. Dirac-type local orbitals are                                                             
!> useful in the context of spin-orbit coupling, but should only be used along                                                          
!> with second variation with local orbitals. For more details see:                                                                     
!> arXiv:2306.02965 [cond-mat.mtrl-sci]. 
!> For each local-orbital, a linear combination of {\tt lorbord} radial functions is                                                    
!> constructed such that its radial derivatives up to order ${\tt lorbord}-1$ are zero                                                  
!> at the muffin-tin radius. This function is normalized. The resulting matching coefficients
!> and normalized radial functions are stored in sepperate arrays 
subroutine get_normalized_radial_functions_and_matching_coefficients(matching_coefficients, p0, p1)

    external polynom

    !> local orbital matching coefficients
    real(dp), intent(out) :: matching_coefficients(natmtot, nlomax, maxlorbord)
    !> radial function multiplied by r
    real(dp), intent(out) :: p0(nrmtmax, natmtot, nlomax, maxlorbord)
    !> first radial derivative of p0
    real(dp), intent(out) :: p1(nrmtmax, natmtot, nlomax, maxlorbord)

    ! local variables
    integer(i32) :: l
    integer(i32) :: is, ia, ias
    integer(i32) :: ilo, io1, io2
    integer(i32) :: j, np
    integer(i32) :: ir, nr
    integer(i32) :: info
    integer(i32) :: nn
    real(dp) :: vr(nrmtmax)
    real(dp) :: norm_coeff
    real(dp) :: p0_temp(nrmtmax, maxlorbord), p1_temp(nrmtmax, maxlorbord)
    real(dp) :: q0_temp(nrmtmax, maxlorbord), q1_temp(nrmtmax, maxlorbord)
    real(dp) :: p0s(nrmtmax)
    real(dp) :: r_inv(nrmtmax)
    real(dp) :: fr(nrmtmax), gr(nrmtmax), cf(3, nrmtmax)
    real(dp) :: polynom
 
    ! allocatables
    real(dp), allocatable :: r_polynom(:), u_polynom(:)
    real(dp), allocatable :: a(:, :), b(:), c(:)
    real(dp), allocatable :: ipiv(:)

    ! number of points to which the polynomial is fitted 
    np = max(maxlorbord+1, 4)
    allocate(ipiv(np))
    allocate(r_polynom(np), u_polynom(np), c(np))
    allocate(a(np, np), b(np))

    do is = 1, nspecies
        nr = nrmt(is)
        r_inv(1:nr) = spr(1:nr, is)
        do ia = 1, natoms(is)
            ias = idxas(ia, is)
            if (associated(input%groundstate%mgga) .and. iscl > 1) then 
               vr(1:nr) = veffmt_gga(1, 1:nr, ias) * y00
            else 
               vr (1:nr) = veffmt(1, 1:nr, ias) * y00
            end if
            do ilo = 1, nlorb(is)
                l = lorbl(ilo, is)
                do io2 = 1, lorbord(ilo, is)

                    if (wfkappa(io2, ilo, is) /= 0) then
                        ! integrate the radial Dirac equation                                                                                                           
                        call rdiracdme(lorbdm(io2, ilo, is), wfkappa(io2, ilo, is), lorbe(io2, ilo, ias), nr, &
                                       spr(:, is), vr, nn, p0_temp(:, io2), p1_temp(:, io2), q0_temp(:, io2), &
                                       q1_temp(:, io2), .false.)
                    else
                        ! integrate the radial Schrodinger equation
                        call rschroddme(lorbdm(io2, ilo, is), l, 0, lorbe(io2, ilo, ias), nr, spr(:, is), vr, nn, &
                                        p0_temp(:, io2), p1_temp(:, io2), q0_temp(:, io2), q1_temp(:, io2))
                    end if

                    ! normalize radial functions
                    fr(1:nr) = p0_temp(1:nr, io2) ** 2
                    ! integrate [u(r)*r]^2 over r
                    call fderiv (-1, nr, spr(1:nr, is), fr, gr, cf)
                    norm_coeff = 1.0_dp / sqrt( abs( gr(nr) ) )
                    p0_temp(1:nr, io2) = norm_coeff * p0_temp(1:nr, io2)
                    p1_temp(1:nr, io2) = norm_coeff * p1_temp(1:nr, io2)

                    p0(1:nr, ias, ilo, io2) = p0_temp(1:nr, io2)
                    p1(1:nr, ias, ilo, io2) = p1_temp(1:nr, io2)

                    ! set up the matrix of radial derivatives
                    r_polynom(1:np) = spr(nr-np+1:nr, is)
                    u_polynom(1:np) = p0_temp(nr-np+1:nr, io2) * r_inv(nr-np+1:nr)

                    do io1 = 1, lorbord(ilo, is)
                        a(io1, io2) = polynom(io1-1, np, r_polynom, u_polynom, c, rmt(is))
                    end do ! io1
                    
                end do ! io2

                ! set up the target vector
                b(:) = 0.0_dp
                b(lorbord(ilo, is)) = 1.0_dp

                ! Solve system of linear equations
                call dgesv(lorbord(ilo, is), 1, a, np, ipiv, b, np, info)

                if (info /= 0) then
                    write (*,*)
                    write (*, '("Error(genlofr): degenerate local-orbital radial functions")')
                    write (*, '(" for species ", I4)') is
                    write (*, '(" atom ", I4)') ia
                    write (*, '(" and local-orbital ", I4)') ilo
                    write (*, '(" ZGESV returned INFO = ", I8)') info
                    write (*,*)
                    call terminate
                end if

                ! generate linear superposition of radial functions
                p0s(1:nr) = 0.0_dp
                do io1 = 1, lorbord(ilo, is)
                    p0s(1:nr) = p0s(1:nr) + b(io1) * p0_temp(1:nr, io1)
                end do
                ! normalize local orbital
                fr(1:nr) = p0s(1:nr) ** 2
                call fderiv (-1, nr, spr(1:nr, is), fr, gr, cf)

                matching_coefficients(ias, ilo, 1:lorbord(ilo, is)) = b(1:lorbord(ilo, is)) / sqrt(abs(gr(nr)))

            end do !ilo
        end do !ia
    end do !is

    end subroutine get_normalized_radial_functions_and_matching_coefficients

    !> Construct the local orbitals from the radial functions and matching coefficients. 
    subroutine compute_local_orbitals(matching_coefficients, p0, p1, local_orbitals)

        !> local orbitals matching coefficients
        real(dp), intent(in) :: matching_coefficients(natmtot, nlomax, maxlorbord)
        !> radial functions times r
        real(dp), intent(in) :: p0(nrmtmax, natmtot, nlomax, maxlorbord)
        !> first radial derivative of p0
        real(dp), intent(in) :: p1(nrmtmax, natmtot, nlomax, maxlorbord)
        !> local orbitals
        real(dp), intent(out) :: local_orbitals(nrmtmax, 2, nlomax, natmtot)

        ! local variable 
        real(dp) :: p1s(nrmtmax)
        real(dp) :: p0s(nrmtmax)
        real(dp) :: r_inv(nrmtmax)
        real(dp) :: t1
        integer(i32) :: is, ia, ias, ilo, io1, ir, nr

        do is = 1, nspecies
            nr = nrmt (is)
            r_inv(1:nr) = 1/spr(1:nr, is)
            do ia = 1, natoms (is)
                ias = idxas (ia, is)
                do ilo = 1, nlorb (is)
                    p0s(1:nr) = 0.0_dp
                    p1s(1:nr) = 0.0_dp
                    do io1 = 1, lorbord(ilo, is)
                        p0s(1:nr) = p0s(1:nr) + matching_coefficients(ias, ilo, io1) * p0(1:nr, ias, ilo, io1)
                        p1s(1:nr) = p1s(1:nr) + matching_coefficients(ias, ilo, io1) * p1(1:nr, ias, ilo, io1)
                    end do !io1

                    local_orbitals(1:nr, 1, ilo, ias) = p0s (1:nr) * r_inv(1:nr)
                    local_orbitals(1:nr, 2, ilo, ias) = (p1s(1:nr)-p0s(1:nr)*r_inv(1:nr)) * r_inv(1:nr)

                end do !ilo
            end do !ia
        end do !is  

    end subroutine compute_local_orbitals

    !> Computes radial functions and matching coefficients, then constructs the local orbitals and saves them in 
    !> the global array lofr. 
    subroutine genlofr
        use mod_APW_LO, only: lofr
        use mod_timing, only: stopwatch

        ! local variables 
        real(dp) :: matching_coefficients(natmtot, nlomax, maxlorbord)
        real(dp) :: p0(nrmtmax, natmtot, nlomax, maxlorbord)
        real(dp) :: p1(nrmtmax, natmtot, nlomax, maxlorbord)

        call stopwatch("exciting:genlofr", 1)

        call get_normalized_radial_functions_and_matching_coefficients(matching_coefficients, p0, p1)
        call compute_local_orbitals(matching_coefficients, p0, p1, lofr)

        call stopwatch("exciting:genlofr", 0)

    end subroutine genlofr

end module mod_gen_lo
