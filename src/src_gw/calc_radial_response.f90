!> This module contains subroutines that compute radial basis response to build the lapw and lo response. 
!> In order to help with the readability of the code, we show here a number of equations, which are referenced
!> via comments in the code. 
!> \begin{align}
!>     \epsilon^{(1)}_{l \alpha, I} 
!>     &= \langle u_{l\alpha \; p=0}| \chi_{I} | u_{l\alpha \; p=0} \rangle \nonumber \\ 
!>     &= \int \big [ u_{l \alpha \; p=0}(r) \big ]^{2} \; \chi_{I}(r) \; r^2 \; dr
!> \end{align} 

!> \begin{align}
!>     \tilde{H}_{l \alpha, l^{\prime}}(r; \omega) &= \frac{l^{\prime}(l^{\prime}+1)}{2r^2} + 
!>      v^{\alpha}_{\text{eff}, 0}(r) - \epsilon_{l \alpha} - \omega
!> \end{align}

!> \begin{equation}
!>     \tilde{H}_{l \alpha, l^{\prime}}(r; \omega)\; r \; u^{\text{hom}}_{l \alpha, l^{\prime}}(r; \omega) = 0
!> \end{equation}

!> \begin{equation}
!>     \mathscr{I}_{l \alpha 0, Il^{\prime}}(r, \omega) = \bigg [ \delta_{ll^{\prime}} \epsilon^{(1)}_{l \alpha, I} 
!>      - \chi_{I}(r) \bigg ] \;  u_{l\alpha 0}(r) \; r 
!> \end{equation}

!> \begin{equation}
!>     \mathscr{I}_{l \alpha p, I l^{\prime}}(r, \omega) = \bigg [ \delta_{ll^{\prime}}\epsilon^{(1)}_{l \alpha, I} 
!>      - \chi_{I}(r) \bigg ]  u_{l\alpha p}(r) \; r + p \; r \; u^{(1)}_{l \alpha p-1, I l^{\prime}}(r; \omega)
!> \end{equation}

!> \begin{equation}
!>     \tilde{H}_{l \alpha, l^{\prime}}(r; \omega)\; r \; u^{(1)}_{l \alpha p, I l^{\prime} } (r; \omega) 
!>     = \bigg [ \delta_{ll^{\prime}} \epsilon^{(1)}_{l \alpha, I} - \chi_{I}(r) \bigg ] 
!>      u_{l\alpha p}(r) \; r  + p \; r \;  u^{(1)}_{l \alpha \;p-1, I l^{\prime}}(r; \omega)
!> \end{equation}

!> \begin{equation}
!>    \begin{bmatrix}
!>     u^{\mathrm{hom}}_{l \alpha, l^{\prime}}(r_{\mathrm{MT}}; \omega) & u_{l^{\prime}\alpha p=0}(r_{\mathrm{MT}}) \\[6pt]
!>     u^{\prime \mathrm{hom}}_{l \alpha, l^{\prime}}(r_{\mathrm{MT}}; \omega) & u^{\prime}_{l^{\prime}\alpha p=0}(r_{\mathrm{MT}})
!>     \end{bmatrix}
!>     \begin{bmatrix}
!>     \alpha \\[6pt] \beta
!>     \end{bmatrix}
!>     = -
!>     \begin{bmatrix}
!>     u^{(1)}_{l \alpha p, I l^{\prime}}(r_{\mathrm{MT}}; \omega) \\[6pt]
!>     u^{\prime (1)}_{l \alpha p, I l^{\prime}}(r_{\mathrm{MT}}; \omega)
!>     \end{bmatrix}
!> \end{equation}

!> \begin{equation}
!>     \tilde{u}^{(1)}_{l \alpha p, I l^{\prime}}(r; \omega)
!>     = u^{(1)}_{l \alpha p, I l^{\prime}}(r; \omega) + 
!>     \alpha \; u^{\mathrm{hom}}_{l \alpha, l^{\prime}}(r; \omega) +
!>     \beta \; u_{l^{\prime}\alpha p=0}(r) 
!> \end{equation}

!> \begin{align}
!>     \tilde{H}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r; \omega) &= \frac{l^{\prime}(l^{\prime}+1)}{2r^2} + 
!>      v^{\alpha}_{\text{eff}, 0}(r) - \epsilon_{l_{\zeta} \alpha \zeta} - \omega
!> \end{align}

!> \begin{align}
!>     \epsilon^{(1)}_{l_{\zeta} \alpha \zeta, I} 
!>     &= \langle u_{l_{\zeta}\alpha \; p=0}(\epsilon_{l_{\zeta}\alpha \zeta})| \chi_{I} | 
!>      u_{l_{\zeta}\alpha \; p=0}(\epsilon_{l_{\zeta}\alpha \zeta}) \rangle \nonumber \\ 
!>     &= \int \big [ u_{l_{\zeta} \alpha \; p=0}(r, \epsilon_{l_{\zeta}\alpha \zeta}) 
!>      \big ]^{2} \; \chi_{I}(r) \; r^2 \; dr
!> \end{align} 

!> \begin{equation}
!>     \tilde{H}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r; \omega)\; r \; 
!>      u^{\text{hom}}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r; \omega) = 0
!> \end{equation}

!> \begin{equation}
!>     \mathscr{I}_{l_{\zeta} \alpha \zeta \;0, Il^{\prime}}(r, \omega) = 
!>      \bigg [ \delta_{l_{\zeta}l^{\prime}} \epsilon^{(1)}_{l_{\zeta} \alpha \zeta, I} - 
!>      \chi_{I}(r) \bigg ] \;  u_{l_{\zeta}\alpha \zeta 0}(r, \epsilon_{l_{\zeta}\alpha\zeta}) \; r 
!> \end{equation}

!> \begin{equation}
!>     \mathscr{I}_{l_{\zeta} \alpha \zeta p, I l^{\prime}}(r, \omega) = 
!>      \bigg [ \delta_{l_{\zeta}l^{\prime}}\epsilon^{(1)}_{l_{\zeta} \alpha \zeta, I} 
!>      - \chi_{I}(r) \bigg ]  u_{l_{\zeta}\alpha \zeta p}(r, \epsilon_{l_{\zeta}\alpha\zeta}) \; r 
!>      + p \; r \; u^{(1)}_{l_{\zeta} \alpha \zeta p-1, I l^{\prime}}(r; \omega)
!> \end{equation}

!> \begin{equation}
!>     \tilde{H}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r; \omega)\; 
!>      r \; u^{(1)}_{l_{\zeta} \alpha \zeta p, I l^{\prime} } (r; \omega) 
!>     = \mathscr{I}_{l_{\zeta} \alpha \zeta p, I l^{\prime}}(r, \omega)
!> \end{equation}

!> \begin{equation}
!>     \begin{bmatrix}
!>     u^{\mathrm{hom}}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r_{\mathrm{MT}}; \omega) 
!>      & u_{l^{\prime}\alpha \zeta p=0}(r_{\mathrm{MT}}, \epsilon_{l_{\zeta}\alpha\zeta}) \\[6pt]
!>     u^{\prime \mathrm{hom}}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r_{\mathrm{MT}}; \omega) 
!>      & u^{\prime}_{l^{\prime}\alpha \zeta p=0}(r_{\mathrm{MT}}, \epsilon_{l_{\zeta}\alpha\zeta})
!>     \end{bmatrix}
!>     \begin{bmatrix}
!>     \alpha \\[6pt] \beta
!>     \end{bmatrix}
!>     = -
!>     \begin{bmatrix}
!>     u^{(1)}_{l_{\zeta} \alpha \zeta p_{\zeta}, I l^{\prime}}(r_{\mathrm{MT}}; \omega) \\[6pt]
!>     u^{\prime (1)}_{l_{\zeta} \alpha \zeta p_{\zeta}, I l^{\prime}}(r_{\mathrm{MT}}; \omega)
!>     \end{bmatrix}
!> \end{equation}

!> \begin{equation}
!>     \tilde{u}^{(1)}_{l_{\zeta} \alpha \zeta p_{\zeta}, I l^{\prime}}(r; \omega)
!>     = u^{(1)}_{l_{\zeta} \alpha \zeta p_{\zeta}, I l^{\prime}}(r; \omega) + 
!>     \alpha \; u^{\mathrm{hom}}_{l_{\zeta} \alpha \zeta, l^{\prime}}(r; \omega) +
!>     \beta \; u_{l^{\prime}\alpha \zeta p=0}(r, \epsilon_{l_{\zeta}\alpha\zeta}) 
!> \end{equation}

!> \begin{equation}
!>     \phi^{(1)}_{\mu, I l^{\prime}}(r;\omega) = \sum_{\zeta} a_{\mu \zeta} \; 
!>      \tilde{u}^{(1)}_{l_{\zeta} \alpha \zeta p_{\zeta}, I l^{\prime}}(r; \omega)
!> \end{equation}
module calc_radial_response
#include "asserts.fpp"
    use precision, only: dp, i32
    use mod_atoms, only: idxas, spr, natmtot, natoms, nspecies
    use mod_muffin_tin, only: nrmt, nrmtmax, rmt
    use mod_APW_LO, only: apwfr, apwe, apwdm, apword, maxapword, lorbe, lorbdm, &
                          lorbord, maxlorbord, nlomax, lorbl, nlorb
    use mod_product_basis, only: umix, nmix, bigl
    use modinput, only: input
    use constants, only: y00, maxlapw
    use mod_potential_and_density, only: veffmt
    use mgga_poteff, only: veffmt_gga
    use mod_convergence, only: iscl
    use errors_warnings, only: terminate_if_false
    use mod_sternheimer_u, only:  build_inhomogeneity, build_h_tilde, build_epsilon_perturbed, &
                                  normalize, compute_sternheimer_u
    use mod_gen_lo, only: get_normalized_radial_functions_and_matching_coefficients
    use linear_algebra_2d, only: solve_2d_cramer
    use modmpi, only: mpiglobal

    implicit none 

    private
    public :: compute_radial_lapw_response, compute_radial_lo_response

contains

    !> Compute the radial response which later builds the lapw response.
    subroutine compute_radial_lapw_response(irm, omega, apwfr_response)

        !> index of the radial part of the muffin tin mixed product basis function
        integer(i32),   intent(in)  :: irm
        !> frequency
        real(dp),       intent(in)  :: omega
        !> radial response and first radial derivative of radial response
        real(dp),       intent(out) :: apwfr_response(nrmtmax, 2, maxapword, 0:maxlapw, 0:maxlapw, natmtot)

        ! local variables
        integer(i32)    :: nr, ir
        integer(i32)    :: is, ia, ias
        integer(i32)    :: l1, l2, bl
        integer(i32)    :: p, iom, io1
        integer(i32)    :: l2_min, l2_max
        integer(i32)    :: info
        real(dp)        :: r_inv(nrmtmax)
        real(dp)        :: epsilon_perturbed
        real(dp)        :: h_tilde(nrmtmax), inhomogeneity(nrmtmax)
        real(dp)        :: vr(nrmtmax), bare_inhomogeneity(nrmtmax)
        real(dp)        :: r_apwfr_homogeneous(nrmtmax, 2)
        real(dp)        :: r_apwfr_inhomogeneous(nrmtmax, 2, maxapword)
        real(dp)        :: r_polynom(2), u_polynom(2)
        real(dp)        :: u_inhom_polynom(2), u_hom_polynom(2)
        real(dp)        :: a(2, 2), b(2), c(2), x(2)

        real(dp), external :: polynom
        character(len=*), parameter :: assert_msg = &
          'The radial functions in an LAPW basis function must be ordered ' // &
          'by the order of energy derivative (p), starting with p=0.'
        
        apwfr_response = 0.0_dp

        do is=1, nspecies
            nr = nrmt(is)
            r_inv(1:nr) = 1.0_dp / spr(1:nr, is) 
            do ia=1, natoms(is)                
                ias = idxas(ia, is)

                if (associated(input%groundstate%mgga) .and. iscl > 1) then 
                    vr(1:nr) = veffmt_gga(1, 1:nr, ias) * y00
                else 
                    vr (1:nr) = veffmt(1, 1:nr, ias) * y00
                end if

                bl =  bigl(irm,ias)
                do l1 = 0, input%groundstate%lmaxapw
                
                    ! Eq.1
                    epsilon_perturbed = build_epsilon_perturbed(umix(1:nr,irm,ias), &
                                                                apwfr(1:nr, 1, 1, l1, ias), &
                                                                spr(1:nr, is), nrmt(is))

                    ! Only loop over combinations of l1, l2, and bl for which non-zero Gaunt coefficients exist.
                    l2_min = abs(bl-l1)
                    l2_max = min(bl+l1,input%groundstate%lmaxGwIbc)        
                    do l2 = l2_min, l2_max
                        ! Eq.2
                        h_tilde(1:nr) = build_h_tilde(l2, apwe(1, l1, ias), omega, nr, vr(1:nr), spr(1:nr, is))

                        ! Eq.3
                        inhomogeneity(1:nr) = 0.0_dp
                        call compute_sternheimer_u(nr, spr(1:nr, is), h_tilde(1:nr), inhomogeneity(1:nr), &
                                                   r_apwfr_homogeneous(1:nr, 1:2))
                        call normalize(r_apwfr_homogeneous(1:nr, 1:2), spr(1:nr, is), nr)

                        do iom = 1, apword(l1, is)
                            ! Eq.4
                            bare_inhomogeneity(1:nr) = build_inhomogeneity(l1, l2, umix(1:nr,irm,ias), &
                                                                           spr(1:nr, is) * apwfr(1:nr, 1, iom, l1, ias), &
                                                                           epsilon_perturbed, nrmt(is))
                            
                            p = apwdm(iom, l1, is)

                            CALL_ASSERT(iom == p + 1, assert_msg )

                            if (p /= 0)  then
                                ! Eq.5
                                ! Stored in apwfr_inhomogeneous(1:nr, 1, iom-1) is the result from the last iom iteration. 
                                ! This should be the solution of the inhomogeneous equation with matching order-1. 
                                ! This is only true for LAPWs. 
                                inhomogeneity(1:nr) = bare_inhomogeneity(1:nr) +  p * r_apwfr_inhomogeneous(1:nr, 1, iom-1)
                            else 
                                inhomogeneity(1:nr) = bare_inhomogeneity(1:nr)
                            end if

                            ! Eq.6
                            call compute_sternheimer_u(nr, spr(1:nr, is), h_tilde(1:nr), inhomogeneity(1:nr), &
                                                       r_apwfr_inhomogeneous(1:nr, 1:2 , iom))

                            r_polynom(1:2) = spr(nr-1:nr, is)
                            u_hom_polynom(1:2) = r_apwfr_homogeneous(nr-1:nr, 1) * r_inv(nr-1:nr)
                            u_polynom(1:2) = apwfr(nr-1:nr, 1, 1, l2, ias)
                            u_inhom_polynom(1:2) = r_apwfr_inhomogeneous(nr-1:nr, 1, iom) * r_inv(nr-1:nr)

                            ! Eq. 7
                            do io1 = 1, 2
                                a(io1, 1) = polynom(io1-1, 2, r_polynom, u_hom_polynom, c, rmt(is))
                                a(io1, 2) = polynom(io1-1, 2, r_polynom, u_polynom, c, rmt(is))
                                b(io1) = - polynom(io1-1, 2, r_polynom, u_inhom_polynom, c, rmt(is))
                            end do ! io1

                            ! Solve system of linear equations
                            call solve_2d_cramer(a, b, x, info)

                            call terminate_if_false(mpiglobal, info == 0, &
                            "Error(compute_radial_lapw_response): degenerate lapw radial response." )
                            
                            ! Eq. 8
                            do ir=1,nr
                                apwfr_response(ir, 1:2, iom, l1, l2, ias) = r_apwfr_inhomogeneous(ir, 1:2, iom) * r_inv(ir) &
                                                                            + x(1) * r_apwfr_homogeneous(ir, 1:2) * r_inv(ir) &
                                                                            + x(2) * apwfr(ir, 1:2, 1, l2, ias)
                            end do 
                        end do !iom
                    end do !l2
                end do !l1
            end do !ia
        end do !is

    end subroutine compute_radial_lapw_response

    !> Compute the local orbital response.
    subroutine compute_radial_lo_response(irm, omega, lofr_response)

        !> index of the radial part of the muffin tin mixed product basis function
        integer(i32),   intent(in)  :: irm
        !> frequency
        real(dp),       intent(in)  :: omega
        !> local orbital response and first radial derivative of local orbital response
        real(dp),       intent(out) :: lofr_response(nrmtmax, 2, nlomax, 0:maxlapw, natmtot)

        ! local variables
        integer(i32)    :: nr, ir
        integer(i32)    :: is, ia, ias
        integer(i32)    :: l1, l2, bl
        integer(i32)    :: ilo
        integer(i32)    :: p, iom, io1
        integer(i32)    :: l2_min, l2_max
        integer(i32)    :: info, nn
        real(dp)        :: r_inv(nrmtmax)
        real(dp)        :: epsilon_perturbed
        real(dp)        :: h_tilde(nrmtmax), inhomogeneity(nrmtmax)
        real(dp)        :: vr(nrmtmax), bare_inhomogeneity(nrmtmax)
        real(dp)        :: lo_u(nrmtmax, 1:2, natmtot, nlomax, maxlorbord)
        real(dp)        :: matching_coefficients(natmtot, nlomax, maxlorbord)
        real(dp)        :: lo_r_u_homogeneous(nrmtmax, 2), lo_r_u_inhomogeneous(nrmtmax, 2)
        real(dp)        :: p0_temp(nrmtmax), p1_temp(nrmtmax), q0_temp(nrmtmax), q1_temp(nrmtmax)
        real(dp)        :: lo_u_response(nrmtmax, 2)
        real(dp)        :: r_polynom(2), u_polynom(2), u_inhom_polynom(2), u_hom_polynom(2)
        real(dp)        :: a(2, 2), b(2), c(2), x(2)

        real(dp), external        :: polynom

        lofr_response = 0.0_dp

        ! Note: This routine outputs u*r, so we have to still divide by r
        call get_normalized_radial_functions_and_matching_coefficients(matching_coefficients, &
                                                                       lo_u(1:nrmtmax, 1, 1:natmtot, 1:nlomax, 1:maxlorbord), &
                                                                       lo_u(1:nrmtmax, 2, 1:natmtot, 1:nlomax, 1:maxlorbord))

        do is=1, nspecies
            nr = nrmt(is)
            r_inv(1:nr) = 1.0_dp / spr(1:nr, is)
            do ia=1, natoms(is)                
                ias = idxas(ia, is)

                ! Note: Even though this is referred to as lo_u, it is still multiplied by r, so we have to divide it here.
                do ir=1, nr
                    lo_u(ir, 1:2, ias, 1:nlomax, 1:maxlorbord) = lo_u(ir, 1:2, ias, 1:nlomax, 1:maxlorbord) * r_inv(ir)
                end do
                
                if (associated(input%groundstate%mgga) .and. iscl > 1) then 
                    vr(1:nr) = veffmt_gga(1, 1:nr, ias) * y00
                else 
                    vr (1:nr) = veffmt(1, 1:nr, ias) * y00
                end if

                bl =  bigl(irm,ias)

                do ilo = 1, nlorb(is)
                    l1 = lorbl(ilo, is)

                    ! Only loop over combinations of l1, l2, and bl for which non-zero Gaunt coefficients exist.
                    l2_min = abs(bl-l1)
                    l2_max = min(bl+l1,input%groundstate%lmaxGwIbc)        
                    do l2 = l2_min, l2_max

                        do iom = 1, lorbord(ilo, is)
                            ! Eq.9
                            h_tilde(1:nr) = build_h_tilde(l2, lorbe(iom, ilo, ias), omega, nr, vr(1:nr), spr(1:nr, is))

                            call rschroddme(0, l1, 0, lorbe(iom, ilo, ias), nr, spr(1:nr, is), vr, nn, &
                                    p0_temp(1:nr), p1_temp(1:nr), q0_temp(1:nr), q1_temp(1:nr))
                            ! Eq.10
                            epsilon_perturbed = build_epsilon_perturbed(umix(1:nr,irm,ias), p0_temp(1:nr) * r_inv(1:nr), &
                                                                        spr(1:nr, is), nrmt(is))

                            ! Eq.11
                            inhomogeneity(1:nr) = 0.0_dp
                            call compute_sternheimer_u(nr, spr(1:nr, is), h_tilde(1:nr), inhomogeneity(1:nr), &
                                                       lo_r_u_homogeneous(1:nr, 1:2))
                            call normalize(lo_r_u_homogeneous(1:nr, 1:2), spr(1:nr, is), nr)
                            
                            ! Eq.12
                            bare_inhomogeneity(1:nr) = build_inhomogeneity(l1, l2, umix(1:nr,irm,ias), p0_temp(1:nr), &
                                                                           epsilon_perturbed, nrmt(is))

                            do p = 0, lorbdm(iom, ilo, is) !ToDo: lots of duplicate calculations. 
                                
                                if (p /= 0)  then                 
                                    ! Eq.13
                                    call rschroddme(p, l1, 0, lorbe(iom, ilo, ias), nr, spr(1:nr, is), vr, nn, &
                                                    p0_temp(1:nr), p1_temp(1:nr), q0_temp(1:nr), q1_temp(1:nr))
                                    bare_inhomogeneity(1:nr) = build_inhomogeneity(l1, l2, umix(1:nr,irm,ias), &
                                                                                   p0_temp(1:nr), &
                                                                                   epsilon_perturbed, nrmt(is))
                                    inhomogeneity(1:nr) = bare_inhomogeneity(1:nr) +  p * lo_r_u_inhomogeneous(1:nr, 1)
                                else 
                                    inhomogeneity(1:nr) = bare_inhomogeneity(1:nr)
                                end if

                                ! Eq.14
                                call compute_sternheimer_u(nr, spr(1:nr, is), h_tilde(1:nr), inhomogeneity(1:nr), &
                                                           lo_r_u_inhomogeneous(1:nr, 1:2))

                            end do !p

                            ! Eq.15
                            r_polynom(1:2) = spr(nr-1:nr, is)
                            u_hom_polynom(1:2) = lo_r_u_homogeneous(nr-1:nr, 1) * r_inv(nr-1:nr)
                            u_polynom(1:2) = lo_u(nr-1:nr, 1, ias, ilo, iom)
                            u_inhom_polynom(1:2) = lo_r_u_inhomogeneous(nr-1:nr, 1) * r_inv(nr-1:nr)

                            do io1 = 1, 2
                                a(io1, 1) = polynom(io1-1, 2, r_polynom, u_hom_polynom, c, rmt(is))
                                a(io1, 2) = polynom(io1-1, 2, r_polynom, u_polynom, c, rmt(is))
                                b(io1) = - polynom(io1-1, 2, r_polynom, u_inhom_polynom, c, rmt(is))
                            end do ! io1

                            ! Solve system of linear equations
                            call solve_2d_cramer(a, b, x, info)

                            call terminate_if_false(mpiglobal, info == 0, &
                            "Error(compute_radial_lo_response): degenerate local orbital response." )

                            ! Eq.16
                            lo_u_response(1:nr, 1) = lo_r_u_inhomogeneous(1:nr, 1) * r_inv(1:nr) &
                                                     + x(1) * lo_r_u_homogeneous(1:nr, 1) * r_inv(1:nr) &
                                                     + x(2) * lo_u(1:nr, 1, ias, ilo, iom)
                            lo_u_response(1:nr, 2) = lo_r_u_inhomogeneous(1:nr, 2) * r_inv(1:nr) &
                                                     + x(1) * lo_r_u_homogeneous(1:nr, 2) * r_inv(1:nr) &
                                                     + x(2) * lo_u(1:nr, 2, ias, ilo, iom)
                            
                            ! Eq.17
                            do ir=1,nr
                                lofr_response(ir, 1:2, ilo, l2, ias) = lofr_response(ir, 1:2, ilo, l2, ias) &
                                                                       + matching_coefficients(ias, ilo, iom) &
                                                                       * lo_u_response(ir, 1:2)
                            end do !ir 
                        end do !iom
                    end do !l2
                end do !ilo
            end do !ia
        end do !is

    end subroutine compute_radial_lo_response
    
end module calc_radial_response
