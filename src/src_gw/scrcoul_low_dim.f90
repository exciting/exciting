!>Module that uses the exact \(q\to 0\) limit of the screened potential in the case where the 2d Coulomb cutoff is applied (PhysRevB.94.155406).
!>
!>The correlation part of the self energy is defined as
!>\begin{equation}
!>\label{eq:selfc}\Sigma^c_{n\mathbf{k}}(\omega) = \frac{1}{N_c}\sum_{\mathbf{q}}\sum_m\frac{i}{2\pi}\int_{-\infty}^{\infty} d\omega^{\prime} \frac{X_{nm}(\mathbf{k},\mathbf{q};\omega^{\prime})}{\omega +\omega^{\prime}-\tilde{\epsilon}_{m\mathbf{k-q}}}
!>\end{equation}
!>with the elements \( X_{nm}\) defined as:
!>\begin{equation}
!>X_{nm}(\mathbf{k},\mathbf{q};\omega)=\sum_{\mu\nu}\left[M^{\mu}_{nm}(\mathbf{k},\mathbf{q})\right]^*W^c_{\mu\nu}(\mathbf{q},\omega)M^{\nu}_{nm}(\mathbf{k},\mathbf{q}).\label{eq:xnm}
!>\end{equation}
!>In the 2D limit the correlation part of the screened Coulomb potential in the coulomb diagonal basis can be written as
!>\begin{align}
!>W_{00}^c(\mathbf{q}\to 0) &= -\left(\frac{4\pi(1-e^{-qL/2})}{q^2}\right)^2\frac{\hat{\mathbf{q}}\cdot A \hat{\mathbf{q}}}{1+4\pi(1-e^{-qL/2})\hat{\mathbf{q}}\cdot A \hat{\mathbf{q}}}\label{eq:whead}\\
!>W_{\mu 0}^c(\mathbf{q}\to 0) &= -\frac{4\pi(1-e^{-qL/2})}{q^2}\frac{\sqrt{v_\mu}\hat{\mathbf{q}}\cdot \mathbf{a}_\mu}{1+4\pi(1-e^{-qL/2})\hat{\mathbf{q}}\cdot A \hat{\mathbf{q}}}\label{eq:wwing}\\
!>W_{0 \mu}^c(\mathbf{q}\to 0) &= -\frac{4\pi(1-e^{-qL/2})}{q^2}\frac{\sqrt{v_\mu}\hat{\mathbf{q}}\cdot \mathbf{b}_\mu}{1+4\pi(1-e^{-qL/2})\hat{\mathbf{q}}\cdot A \hat{\mathbf{q}}}\\
!>W_{\mu\nu}^c(\mathbf{q}\to 0) &= \sqrt{v_\mu v_\nu}\left[B_{\mu\nu}^{-1} -\delta_{\mu\nu} + \frac{4\pi(1-e^{-qL/2})(\hat{\mathbf{q}}\cdot\mathbf{a}_{\mu})(\hat{\mathbf{q}}\cdot\mathbf{b}_{\nu})}{1+4\pi(1-e^{-qL/2})\hat{\mathbf{q}}\cdot A \hat{\mathbf{q}}}\right]\label{eq:wbody}
!>\end{align}
!>with interlayer distance \(L\). For the sake of readability, the frequency dependence is not written explicitly. The vectors \(\mathbf{a}_{\mu}, \mathbf{b}_{\mu}\)
!>and the tensor \(A\) also depend implictly on the frequency. They are defined as
!>
!>\begin{align}
!>\mathbf{a}_{\mu} &= -\sum_{\nu\neq 0}B^{-1}_{\mu\nu}\sqrt{v_\nu}\mathbf{p}_\nu\\
!>\mathbf{b}_{\mu} &= -\sum_{\nu\neq 0}\sqrt{v_\nu}\mathbf{s}_\nu B^{-1}_{\nu\mu}\\
!>A &= -\mathcal{P} - \sum_{\mu\neq 0}\sqrt{v_\mu}\mathbf{s}_\mu\otimes\mathbf{a}_\mu
!>\end{align}
!>
!>and depend on the quantities:
!>\begin{align}
!>\mathbf{p}_\mu &= \frac{1}{\sqrt{\Omega}}\sum_{n,m}\sum_{\mathbf{k}}F_{nm,\mathbf{k}}\frac{\mathbf{p}_{nm,\mathbf{k}}}{\epsilon_{m\mathbf{k}}-\epsilon_{n\mathbf{k}}}\left[M_{nm}^\mu(\mathbf{k},0)\right]^*\\
!>\mathbf{s}_\mu &= \frac{1}{\sqrt{\Omega}}\sum_{n,m}\sum_{\mathbf{k}}F_{nm,\mathbf{k}}\frac{\mathbf{p}_{nm,\mathbf{k}}}{\epsilon_{m\mathbf{k}}-\epsilon_{n\mathbf{k}}}M_{nm}^\mu(\mathbf{k},0)\\
!>\mathcal{P} &= \frac{1}{\Omega}\sum_{n,m}\sum_{\mathbf{k}}F_{nm,\mathbf{k}}\frac{\mathbf{p}_{nm,\mathbf{k}}\otimes \mathbf{p}_{nm,\mathbf{k}}}{(\epsilon_{m\mathbf{k}}-\epsilon_{n\mathbf{k}})^2}.
!>\end{align}
!>
!>The global routine `calcwings` evaluates the vectors \( -\sqrt{4\pi v_\mu}\mathbf{p}_\mu\) and \( -\sqrt{4\pi v_\mu}\mathbf{s}_\mu\) whilst the global routine `calchead`
!>evaluates the tensor \(4\pi\mathcal{P}\). The inverse body \( B_{\mu\nu}^{-1}\) is evaluated within the global routines `calcepsilon_2d` and `calcinveps_2d`.
!>
!>The \( \mathbf{q}=0 \)-term is now treated by replacing \( W^c_{\mu\nu}(0) \) in eq. \ref{eq:xnm} by the average:
!>\begin{equation}
!>\overline{W}_{\mu\nu}^c \equiv \frac{1}{\Omega_0}\int_{\Omega_0}W^c_{\mu\nu}(\mathbf{q})d \mathbf{q}.
!>\end{equation}
!>Where \( \Omega_0 \) is the volume of the Brillouin zone (mini BZ) of the lattice defined by the q-grid. Furthermore we write for the elements  \( M^0_{nm} \):
!>\begin{equation}
!>M^0_{nm}(\mathbf{k},0) = \frac{1}{\sqrt{\Omega}}\delta_{nm}.
!>\end{equation}
!>With this we rewrite eq. \ref{eq:xnm} for \( \mathbf{q}=0 \) as
!>\begin{align}
!>\nonumber X_{nm}(\mathbf{k},\mathbf{q}=0)&=\sum_{\mu\nu\neq 0}\left[M^\mu_{nm}(\mathbf{k},0)\right]^*\overline{W}_{\mu\nu}^cM^\nu_{nm}(\mathbf{k},0)\\
!>\nonumber&+\frac{\delta_{nm}}{\Omega}\overline{W}_{00}^c\\
!>\nonumber&+\frac{\delta_{nm}}{\sqrt{\Omega}} \sum_{\mu\neq 0}\overline{W}_{0\mu}^cM^\mu_{nm}(\mathbf{k},0)\\
!>&+\frac{\delta_{nm}}{\sqrt{\Omega}} \sum_{\mu\neq 0}\overline{W}_{\mu0}^c\left[M^\mu_{nm}(\mathbf{k},0)\right]^*\label{eq:sing_xnm}
!>\end{align}
!>
!>The evaluation of the elements \( X_{nm} \) takes palce in the global routine `calcmwm`. Due to computational reasons a slightly different formula than eq. \ref{eq:xnm} is used:
!>\begin{equation}
!>X_{nm}(\mathbf{k},\mathbf{q};\omega)=\sum_{\mu\nu}\left[\tilde{M}^{\mu}_{nm}(\mathbf{k},\mathbf{q})\right]^*(\epsilon_{\mu\nu}^{-1}(\mathbf{q},\omega)-\delta_{\mu\nu})\tilde{M}^{\nu}_{nm}(\mathbf{k},\mathbf{q})
!>\end{equation}
!>with
!>\begin{equation}
!>\tilde{M}^{\mu}_{nm}(\mathbf{k},\mathbf{q})=\sqrt{v_{\mu}}M^{\mu}_{nm}(\mathbf{k},\mathbf{q}).
!>\end{equation}
!>Therefore, we split of the finite eigenvalues \( \sqrt{v_{\mu}}\) of eqs. \ref{eq:wwing}-\ref{eq:wbody}. The integration can still be performed as they were assumed constant anyways.
!>Then the averaged modified screened Coulomb potential can be stored to the global arrays `epsilon`, `epsw1`, `epsw2` and `epsh`. The last thing that has to be accounted for is that
!>the prefactors in eq. \ref{eq:sing_xnm} are not consistent with the ones in the routine `calcmwm`. This is corrected within the routine `set_singc12`.

module scrcoul_low_dim
    use precision, only: dp
#include "asserts.fpp"
    implicit none
    private
    public :: apply_2d_limit, set_singc12
contains

    !> Calculate the vectors \( \mathbf{a}_\mu \) according to the formula
    !> \[ \mathbf{a}_{\mu} = -\sum_{\nu\neq 0}B^{-1}_{\mu\nu}\sqrt{v_\nu}\mathbf{p}_\nu \]
    !> where the \( \sqrt{v_\mu}\) are already accounted in the \( \mathbf{p}_\nu\)
    subroutine construct_amu(inverse_epsilon_body, wing1, amu)
        use constants, only: zzero, zone
        use general_matrix_multiplication, only: matrix_multiply
        !> Inverse of the body of the dielectric function
        complex(dp), intent(in) :: inverse_epsilon_body(:, :)
        !> 1st wing of the dielectric function
        complex(dp), intent(in) :: wing1(:, :)
        !> Output vector \(\mathbf{a}_\mu\)
        complex(dp), intent(out), allocatable :: amu(:, :)

        integer :: mbsiz

        mbsiz = size(inverse_epsilon_body, 1)
        CALL_ASSERT(size(inverse_epsilon_body, 2) == mbsiz, 'inverse_epsilon_body must be a square matrix')
        CALL_ASSERT(size(wing1, 1) == mbsiz, 'wing1 must have size mbsiz along 1st dim.')
        CALL_ASSERT(size(wing1, 2) == 3, 'wing1 must have size 3 along 2nd dim.')

        allocate (amu(mbsiz, 3))
        call matrix_multiply(inverse_epsilon_body, wing1, amu)
        amu = -amu

    end subroutine

    !> Calculate the vectors \( \mathbf{b}_\mu \) according to the formula
    !> \[ \mathbf{b}_{\mu} = -\sum_{\nu\neq 0}\sqrt{v_\nu}\mathbf{s}_\nu B^{-1}_{\nu\mu} \]
    !> where the \( \sqrt{v_\mu}\) are already accounted in the \( \mathbf{s}_\nu\)
    subroutine construct_bmu(inverse_epsilon_body, wing2, bmu)
        use constants, only: zzero, zone
        use general_matrix_multiplication, only: matrix_multiply
        !> Inverse of the body of the dielectric function
        complex(dp), intent(in) :: inverse_epsilon_body(:, :)
        !> 2nd wing of the dielectric function
        complex(dp), intent(in) :: wing2(:, :)
        !> Output vector \(\mathbf{b}_\mu\)
        complex(dp), intent(out), allocatable :: bmu(:, :)

        integer :: mbsiz

        mbsiz = size(inverse_epsilon_body, 1)
        CALL_ASSERT(size(inverse_epsilon_body, 2) == mbsiz, 'inverse_epsilon_body must be a square matrix')
        CALL_ASSERT(size(wing2, 1) == mbsiz, 'wing2 must have size mbsiz along 1st dim.')
        CALL_ASSERT(size(wing2, 2) == 3, 'wing2 must have size 3 along 2nd dim.')

        allocate (bmu(mbsiz, 3))
        call matrix_multiply(inverse_epsilon_body, wing2, bmu, 't', 'n')
        bmu = -bmu

    end subroutine

!> \[ A = -\mathcal{P} - \sum_{\mu\neq 0}\sqrt{v_\mu}\mathbf{s}_\mu\otimes\mathbf{a}_\mu \]
    subroutine construct_a(head, wing2, amu, a)
        use constants, only: zzero, zone
        use general_matrix_multiplication, only: matrix_multiply
        !> Head of the polarizability \( \mathcal{P} \)
        complex(dp), intent(in) :: head(3, 3)
        !> 2nd wing of the dielectric function \( \mathbf{s}_\mu \)
        complex(dp), intent(in) :: wing2(:, :)
        !> vector \(\mathbf{a}_\mu\)
        complex(dp), intent(in) :: amu(:, :)
        !> Output Tensor \( A \)
        complex(dp), intent(out) :: a(3, 3)

        integer :: mbsiz

        mbsiz = size(wing2, 1)
        CALL_ASSERT(size(wing2, 2) == 3, 'wing2 must have size 3 along 2nd dim.')
        CALL_ASSERT(size(amu, 1) == mbsiz, 'amu must have size mbsiz along 1st dim.')
        CALL_ASSERT(size(amu, 2) == 3, 'amu must have size 3 along 2nd dim.')

        call matrix_multiply(wing2, amu, a, 't', 'n')
        a = -a - head

    end subroutine

    !>Evaluate the head analytically for \( \mathbf{q}=0\) on subgrid. The radius \(r_0\) is
    !>chosen such that \( \pi r_0^2=\frac{\Omega_0}{N_{q0}}\)
    !> \[ -2\pi\int_0^{r_0}\frac{x}{1+(1+4\pi A)x}dx = \frac{-2\pi\left(4\pi Ar_0+r_0-\ln(4\pi r_0+r_0+1)\right)}{(4\pi A+1)^2} \]
    pure function analytic_head(a_avg, r0) result(w_out)
        use constants, only: twopi
        !>Spherical average of the tensor \( A\)
        complex(dp), intent(in) :: a_avg
        !>Radius for analytical integral
        real(dp), intent(in) :: r0
        !>Result of the integral
        complex(dp) :: w_out

        w_out = -twopi*(2*twopi*a_avg*r0 + r0 - log(2*twopi*a_avg*r0 + r0 + 1))/((2*twopi*a_avg + 1)**2)

    end function

    !>Calculate the area each q-point of the subgrid has.
    pure function calc_area_q(bvec_sublattice, ngrid_q0) result(area_q)
        use linear_algebra_3d, only: cross_product
        !> Basis vectors defining the sub lattice
        real(dp), intent(in) :: bvec_sublattice(3, 3)
        !>Number of q-points in each direction for sampling of small region around q=0
        integer, intent(in) :: ngrid_q0(3)
        !> Area around each q-point of the subgrid.
        real(dp) :: area_q

        real(dp) :: ab_norm(3)

        ab_norm = cross_product(bvec_sublattice(:, 1), bvec_sublattice(:, 2))
        area_q = norm2(ab_norm)/dble(product(ngrid_q0))

    end function

    !>Calculate the head of the screened interaction according to \ref{eq:whead}.
    subroutine calculate_screened_interaction_head(averaging_type, q, a, r_cut, area_q, screened_interaction_head)
        use constants, only: pi, fourpi, zzero, zone
        use mod_misc_gw, only: gammapoint
        use modmpi, only: terminate
        !>Specify, whether 1D or 2D averaging should be applied
        character(2), intent(in) :: averaging_type
        !>q-point in cartesian coordinates for which the screened interaction is evaluated
        real(dp), intent(in) :: q(3)
        !>Tensor \(A\)
        complex(dp), intent(in) :: a(3, 3)
        !>Cutoff length of coulomb potential
        real(dp), intent(in) :: r_cut
        !>2D volume each q-point of the subgrid has
        real(dp), intent(in) :: area_q
        !>Head of screened interaction
        complex(dp), intent(out) :: screened_interaction_head

        !local variables
        logical :: gamma
        complex(dp) :: a_avg, analytic_int
        real(dp) :: r0, vcoul_times_qsq, denominator, qhat(3)

        gamma = gammapoint(q)
        if (gamma) then !special treatment for gamma point on sub grid
            select case (averaging_type)
            case ('2d')
                a_avg = 0.5_dp*(a(1, 1) + a(2, 2))
                r0 = sqrt(area_q/pi)*r_cut
                analytic_int = analytic_head(a_avg, r0)
                screened_interaction_head = 16.0_dp*pi**2*a_avg*analytic_int/area_q
                !TODO(Ben) case('1d')
            case default
                call terminate('Only 2D case implemented')

            end select
        else

            select case (averaging_type)
            case ('2d')
                vcoul_times_qsq = fourpi*(1.0_dp - exp(-norm2(q)*r_cut))
                qhat = q/norm2(q)
                a_avg = dot_product(qhat, matmul(a(:, :), qhat))

                !TODO(Ben) case('1d')
            case default
                call terminate('Only 2D case implemented')

            end select
            denominator = zone + vcoul_times_qsq*a_avg
            screened_interaction_head = -(vcoul_times_qsq/norm2(q))**2*a_avg/denominator
        end if

    end subroutine

    !> Calculate analytically the expressions for the q-dependent part of the body of the screened Coulomb potential
    !>\[ \frac{4\pi(1-e^{-qL/2})(\hat{\mathbf{q}}\cdot\mathbf{a}_{\mu})(\hat{\mathbf{q}}\cdot\mathbf{b}_{\nu})}{1+4\pi(1-e^{-qL/2})\hat{\mathbf{q}}\cdot A \hat{\mathbf{q}}}\]
    !>and adds it to the content of the array `screened_interaction_body`. For \( \mu\neq0 \) the eigenvalues of the coulomb potential are split of to the
    !>elements \( \tilde{M}^\mu_{nm}\).
    subroutine calculate_screened_interaction_body(averaging_type, q, amu, bmu, a, r_cut, screened_interaction_body)
        use constants, only: fourpi, zone, zzero
        use mod_misc_gw, only: gammapoint
        use modmpi, only: terminate
        use general_matrix_multiplication, only: matrix_multiply
        use vector_multiplication, only: outer_product

        !>Specify whether 1D or 2D averaging should be applied
        character(2), intent(in) :: averaging_type
        !>q-point in cartesian coordinates for which the screened interaction is evaluated
        real(dp), intent(in) :: q(3)
        !> vector \(\mathbf{a}_\mu\)
        complex(dp), intent(in) :: amu(:, :)
        !> vector \(\mathbf{b}_\mu\)
        complex(dp), intent(in) :: bmu(:, :)
        !> vector \(A\)
        complex(dp), intent(in) :: a(3, 3)
        !>Cutoff length of coulomb potential
        real(dp), intent(in) :: r_cut
        !> Out = In + Body of the screened interaction
        complex(dp), intent(inout) :: screened_interaction_body(:, :)

        integer :: mbsiz
        real(dp) :: vcoul_times_qsq
        complex(dp) :: prefactor_body
        complex(dp), allocatable :: amu_dot_qhat(:), bmu_dot_qhat(:)
        logical :: gamma
        complex(dp) :: a_avg
        complex(dp):: aq_tmp(3), qhat(3)

        mbsiz = size(amu, 1)

        gamma = gammapoint(q)
        if (.not. gamma) then !for gamma point on the sub gird the q-dependent term goes to zero
            select case (averaging_type)
            case ('2d')
                vcoul_times_qsq = fourpi*(1.0_dp - exp(-norm2(q)*r_cut))
                qhat = q/norm2(q)
                aq_tmp = matmul(a(:, :), qhat)
                a_avg = dot_product(qhat, aq_tmp)

                !TODO(Ben) case('1d')
            case default
                call terminate('Only 2D case implemented')

            end select

            allocate (amu_dot_qhat(mbsiz), bmu_dot_qhat(mbsiz))

            call matrix_multiply(amu, qhat, amu_dot_qhat)
            call matrix_multiply(bmu, qhat, bmu_dot_qhat)

            prefactor_body = vcoul_times_qsq/(zone + vcoul_times_qsq*a_avg)
            call outer_product(prefactor_body*amu_dot_qhat(:), bmu_dot_qhat(:), screened_interaction_body) !adds outer product to screened_interaction_body

        end if
    end subroutine

    !>Evaluate the averaged screened coulomb interaction
    !>\[ \overline{W}_{\mu\nu}^c = \frac{1}{\Omega_0}\int_{\Omega_0}W^c_{\mu\nu}(\mathbf{q})dq = \frac{1}{N_{q_0}}\sum_{\mathbf{q}}W^c_{\mu\nu}(\mathbf{q}) \]
    !>Where \( N_{q_0} \) is the number of q-points in the subgrid
    subroutine integrate_screened_interaction(averaging_type, vqc, area_q, inverse_epsilon_body, amu, bmu, a, r_cut, integrated_screened_interaction_head, integrated_screened_interaction_body)
        use constants, only: zzero
        !>Specify, whether 1D or 2D averaging should be applied
        character(2), intent(in) :: averaging_type
        !>set of q-points for the numeric integration
        real(dp), intent(in) :: vqc(:, :)
        !>2D volume each q-point of the subgrid has
        real(dp), intent(in) :: area_q
        !>Inverse of the body of the dielectric function
        complex(dp), intent(in) :: inverse_epsilon_body(:, :)
        !> vector \(\mathbf{a}_\mu\)
        complex(dp), intent(in) :: amu(:, :)
        !> vector \(\mathbf{b}_\mu\)
        complex(dp), intent(in) :: bmu(:, :)
        !> vector \(A\)
        complex(dp), intent(in) :: a(3, 3)
        !>Cutoff length of coulomb potential
        real(dp), intent(in) :: r_cut
        !>Integrated head of screened interaction
        complex(dp), intent(out) :: integrated_screened_interaction_head
        !>Integrated body of screened interaction
        complex(dp), intent(out), allocatable :: integrated_screened_interaction_body(:, :)

        !local variables
        integer :: iq, mbsiz, nqpt
        complex(dp) :: screened_interaction_head
        complex(dp) ::  wqp

        mbsiz = size(inverse_epsilon_body, 1)
        nqpt = size(vqc, 2)

        wqp = 1/dble(nqpt)

        allocate (integrated_screened_interaction_body(mbsiz, mbsiz), source=zzero)

        integrated_screened_interaction_head = zzero

        do iq = 1, nqpt
            call calculate_screened_interaction_head(averaging_type, vqc(:, iq), a, r_cut, area_q, screened_interaction_head)
            call calculate_screened_interaction_body(averaging_type, vqc(:, iq), amu, bmu, a, r_cut, integrated_screened_interaction_body)

            integrated_screened_interaction_head = integrated_screened_interaction_head + wqp*screened_interaction_head

        end do

        integrated_screened_interaction_body(:, :) = wqp*integrated_screened_interaction_body(:, :) + inverse_epsilon_body

    end subroutine

    !>Small hack that redefines the globals `singc1` and `singc2` such that the variables `coefs1` and `coefs2` in routine `calcmwm`
    !>have the right value
    subroutine set_singc12
        use modgw, only: kqset, singc1, singc2
        use constants, only: fourpi

        real(dp) :: wkq
        wkq = 1/dble(kqset%nkpt)
        singc1 = wkq/sqrt(fourpi)
        singc2 = wkq/fourpi
    end subroutine

    !>Apply analytic formulas to calculate the averaged screened Coulomb potential in the limit \( q\to 0 \).
    subroutine apply_2d_limit(bvec, r_cut, symt2, eps, epsh, epsw1, epsw2)
        use mod_kpointset, only: kq_set, generate_kq_vectors
        use modinput, only: input
        use constants, only: fourpi, zone, zzero, pi

        !>Reciprocal basis vectors
        real(dp), intent(in) :: bvec(3, 3)
        !>Cutoff length of coulomb potential
        real(dp), intent(in) :: r_cut
        !>Symmetrization tensor
        real(dp), intent(in) :: symt2(3, 3, 3, 3)
        !>In: \( B_{\mu\nu}^{-1}\), Out: \( \overline{W}_{\mu\nu}^c\)
        complex(dp), intent(inout) :: eps(:, :)
        !>In: \( 4\pi\mathcal{P} \), Out: \( \overline{W}_{00}^c\) (stored in `epsh(1, 1)`)
        complex(dp), intent(inout) :: epsh(:, :)
        !>In: \( -\sqrt{4\pi v_\mu}\mathbf{p}_\mu\), Out: \( \overline{W}_{\mu 0}^c\) (stored in `epsw1(:, 1)`)
        complex(dp), intent(inout) :: epsw1(:, :)
        !>In: \( -\sqrt{4\pi v_\mu}\mathbf{s}_\mu\), Out: \( \overline{W}_{0\mu}^c\) (stored in `epsw2(:, 1)`)
        complex(dp), intent(inout) :: epsw2(:, :)

        complex(dp), allocatable :: amu(:, :), bmu(:, :), integrated_screened_interaction_body(:, :)
        complex(dp) :: integrated_screened_interaction_head
        complex(dp) :: a(3, 3), a_tmp(3, 3)
        type(kq_set) :: kqset_q0
        real(dp) :: bvec_sublattice(3, 3), area_q
        integer :: i, iq, iv(3), iop, jop

        !normalize head and wings to convention used in this module
        do iop = 1, 3
            epsh(iop, iop) = zone - epsh(iop, iop)
        end do
        epsh = epsh/fourpi
        epsw1 = -epsw1/sqrt(fourpi)
        epsw2 = -epsw2/sqrt(fourpi)

        !construct needed vectors and tensors
        call construct_amu(eps(:, :), epsw1(:, :), amu)
        call construct_bmu(eps(:, :), epsw2(:, :), bmu)
        call construct_a(epsh(:, :), epsw2(:, :), amu, a)

        !symmetrize A
        a_tmp(:, :) = a(:, :)
        do iop = 1, 3
            do jop = 1, 3
                call symt2app(iop, jop, 1, symt2, a_tmp, a(iop, jop))
            end do
        end do

        !Define basis vectors of subgrid around q = 0
        do i = 1, 3
            bvec_sublattice(:, i) = bvec(:, i)/dble(input%gw%ngridq(i))
        end do

        call generate_kq_vectors(kqset_q0, &
        &                        bvec_sublattice, &
        &                        input%gw%scrcoul%subgrid_q0, &
        &                        [0.0_dp, 0.0_dp, 0.0_dp], &
        &                        .false.)

        ! Map grid to mini-BZ around q=0
        do iq = 1, kqset_q0%nkpt
            call vecfbz(input%structure%epslat, bvec_sublattice, kqset_q0%vql(:, iq), iv)
            call r3mv(bvec_sublattice, kqset_q0%vql(:, iq), kqset_q0%vqc(:, iq))
        end do

        area_q = calc_area_q(bvec_sublattice, input%gw%scrcoul%subgrid_q0)
        call integrate_screened_interaction(input%gw%scrcoul%averaging, kqset_q0%vqc, area_q, eps(:, :), amu, bmu, a, r_cut, integrated_screened_interaction_head, integrated_screened_interaction_body)

        !adjust globals for further calculation
        call set_singc12
        epsh(1, 1) = integrated_screened_interaction_head
        !wings are zero because analytic expression is odd
        epsw1(:, 1) = zzero
        epsw2(:, 1) = zzero
        eps(:, :) = integrated_screened_interaction_body(:, :)

    end subroutine

end module
