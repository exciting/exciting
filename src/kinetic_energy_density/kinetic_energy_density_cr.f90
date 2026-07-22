!> Module provides subroutines to calculate the kinetic energy density
!> in the muffin-tin region for the core electrons `ked_cr`. 
module kinetic_energy_density_cr

    use precision, only: dp
    use errors_warnings, only: terminate_if_true, terminate_if_false
    use modmpi, only: mpiglobal
    use modinput, only: input
    use mod_muffin_tin, only: nrmt, nrmtmax, idxlm
    use mod_atoms, only: natoms, idxas, nspecies, &
        & spnst, spcore, spr, spk, spocc, spl, natmtot, spnstmax
    use mod_SHT, only: zbshtvr
    use kinetic_energy_density_vars
    use constants, only: zzero, zone
    use mod_corestate, only: rhocr, rwfcr

    implicit none 
    private

    public :: gen_ked_cr, get_clebsch_gordan_coeffs

    contains
    
    !> Returns the kinetic energy density for the core electrons in the muffin-tin region.
    !> This is done by generating the four components of the core wavefunction and then taking the 
    !> gradient for each component. 
    !> From the gradient of the wavefunction, the total kinetic energy density can be 
    !> constructed.
    subroutine gen_ked_cr(ked_cr)
        !> Kinetic energy density for the core electrons
        real(dp), intent(inout) :: ked_cr(:,:,:)

        integer :: ist, nr_core_state, wf_component
        integer :: is, ia, ias
        integer :: nr, m, lm, ir

        complex(dp), allocatable :: ked_cr_per_state(:,:,:,:)
        real(dp), allocatable :: wf_cr_per_state(:, :)

        integer, allocatable :: core_states(:,:), core_states_count(:)

        allocate(ked_cr_per_state(ked_lmmaxvr, nrmtmax, 4, natmtot), source=zzero)
        allocate(wf_cr_per_state(nrmtmax, 4), source=0.0_dp)

        allocate(core_states(spnstmax, nspecies))  
        allocate(core_states_count(nspecies))

        ! For each species, get the number of core states and their indices
        do is = 1, nspecies
            core_states_count(is) = COUNT(spcore(:, is))
            core_states(1:core_states_count(is), is) = PACK([(ist, ist=1, spnst(is))], spcore(1:spnst(is), is))
        end do
    
        ked_cr = 0.0_dp
        do is = 1, nspecies
            nr = nrmt(is)
            do ia = 1, natoms(is)
                ias = idxas(ia, is)
                do nr_core_state = 1, core_states_count(is)
                    ist = core_states(nr_core_state, is) 
                    do m = -spk(ist, is), spk(ist, is)-1
                        call gen_wf_cr_per_state(is, ia, ist, m, wf_cr_per_state) 
                        call gen_ked_cr_per_state(is, ia, ist, m, wf_cr_per_state, ked_cr_per_state(:, :, :, ias))
                        do wf_component = 1, 4
                            do lm = 1, ked_lmmaxvr
                                do ir = 1, nrmt(is)
                                    ked_cr(lm, ir, ias) = ked_cr(lm, ir, ias) + dble(ked_cr_per_state(lm, ir, wf_component, ias)) 
                                end do 
                            end do 
                        end do 
                    end do 
                end do
            end do
        end do

        deallocate(ked_cr_per_state, wf_cr_per_state)
    end subroutine  

    !> Returns the four-component relativistic Dirac wavefunction in spherical harmonic representation  
    !> for atom \(\alpha\) and for state `ist`, defined by the quantum numbers \(\kappa\) and \(m\).  
    !>  
    !> The full wavefunction is given by:  
    !> \[  
    !> \psi^{\alpha}_{\kappa m}(\mathbf{r}) =  
    !>      \begin{pmatrix}  
    !>          g^{\alpha}_{\kappa}(r) \, \Omega_{\kappa m}(\hat{\mathbf{r}}) \\  
    !>          -i f^{\alpha}_{\kappa}(r) \, \Omega_{-\kappa m}(\hat{\mathbf{r}}),  
    !>      \end{pmatrix},  
    !> \]  
    !> where \(g^{\alpha}_{\kappa}(r)\) and \(-i f^{\alpha}_{\kappa}(r)\) are the radial functions  
    !> for the large and small components, respectively, obtained by solving the radial Dirac equations (see gencore.f90).  
    !>  
    !> The quantum number \(\kappa\) is defined as:  
    !> \[  
    !> \kappa =  
    !>      \begin{cases}  
    !>         -l - 1 & \text{for } j = l + \frac{1}{2}, \\  
    !>         l & \text{for } j = l - \frac{1}{2},  
    !>      \end{cases}  
    !> \]  
    !> and the spin spherical harmonics \(\Omega_{lsjm}(\hat{\mathbf{r}})\) with \(s=\frac{1}{2}\) are defined as:  
    !> \[  
    !> \Omega_{ljm}(\theta, \phi) = \sum_{m_{s} = \pm \frac{1}{2}} C^{j m}_{l, m-m_{s}, \frac{1}{2}, m_{s}} Y_{l, m - m_{s}}(\theta, \phi) \, \chi_{m_{s}},  
    !> \]  
    !> where \(C^{j m_{j}}_{l m s m_s}\) are the Clebsch-Gordan coefficients and \(\chi_{\sigma}\) for \(m_{s}= \pm \frac{1}{2}\) is given by:  
    !> \[  
    !> \chi_{\frac{1}{2}} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}, \quad \chi_{-\frac{1}{2}} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}.  
    !> \]  
    !>  
    !> The spin spherical harmonics can then be expressed as:  
    !> \[  
    !> \Omega_{l,l\pm\frac{1}{2},m}(\hat{\mathbf{r}}) =  
    !>    \begin{pmatrix}  
    !>      \pm \sqrt{\frac{l\pm m+\frac{1}{2}}{2l+1}} \, Y_{l,m\mp\frac{1}{2}}(\hat{\mathbf{r}}) \\  
    !>      \sqrt{\frac{l\mp m+\frac{1}{2}}{2l+1}} \, Y_{l,m\pm\frac{1}{2}}(\hat{\mathbf{r}})  
    !>    \end{pmatrix}.  
    !> \]  
    !> In practice, the \(m\) index is provided through `m`, which is defined as \(m - \frac{1}{2}\).  
    !>  
    !> The returned wavefunction is represented in terms of spherical harmonics \(Y_{l,m}(\hat{\mathbf{r}})\),  
    !> taking the following form:  
    !> \[  
    !> \psi^{\alpha}_{\kappa m}(r) =  
    !>      \begin{pmatrix}  
    !>          g^{\alpha}_{\kappa}(r) \, \vec{C}_{l,l\pm\frac{1}{2},m} \\  
    !>          -i f^{\alpha}_{\kappa}(r) \, \vec{C}_{l,l\pm\frac{1}{2},m},  
    !>      \end{pmatrix},  
    !> \]  
    !> where the explicit angular dependence via \(Y_{l,m}(\hat{\mathbf{r}})\) is omitted  
    !> and \(\vec{C}_{l,l\pm\frac{1}{2},m}\) represents the Clebsch-Gordan coefficients of the  
    !> spin spherical harmonics:  
    !> \[  
    !> \vec{C}_{l,l\pm\frac{1}{2},m} =  
    !>    \begin{pmatrix}  
    !>      \pm \sqrt{\frac{l\pm m+\frac{1}{2}}{2l+1}}  \\  
    !>      \sqrt{\frac{l\mp m+\frac{1}{2}}{2l+1}}  
    !>    \end{pmatrix},  
    !> \]  
    !> obtained through the function [[get_diracwf_cg_coeffs(function)]].
    subroutine gen_wf_cr_per_state(is, ia, ist, m, wf_cr_per_state)
        !> species index
        integer, intent(in) :: is
        !> atom index
        integer, intent(in) :: ia
        !> core state
        integer, intent(in) :: ist
        !> magnetic quantum number (passed in as m-1/2)
        integer, intent(in) :: m
        !> core wavefunction 
        real(dp), intent(out) :: wf_cr_per_state(:, :)

        ! local variables
        integer :: nr, ias, k, l
        real(dp) :: wf_cg(2)
        character(:), allocatable :: error_message

        nr = nrmt(is)
        ias = idxas(ia, is)

        l = spl(ist, is)
        k = spk(ist, is)
                
        call terminate_if_false(mpiglobal, input%groundstate%CoreRelativity == "dirac", "currently only works with dirac eq.")
        call wf_cr_error_msg(error_message, l, k, m, is, ia, ist)
        call terminate_if_true(mpiglobal, (m < -k) .or. (m > (k-1)), error_message)
        call terminate_if_false(mpiglobal, (k == (l+1) .or. k == l), error_message)
        
        wf_cg = get_diracwf_cg_coeffs(l, k, m)
        
        wf_cr_per_state = 0.0_dp
        if (abs(m) <= l) then
            wf_cr_per_state(1:nr, 1) = (wf_cg(1) * rwfcr (1:nr, 1, ist, ias)) / spr (1:nr, is)
            wf_cr_per_state(1:nr, 3) = (wf_cg(1) * rwfcr (1:nr, 2, ist, ias)) / spr (1:nr, is)        
        end if
        
        if (abs(m+1) <= l) then
            wf_cr_per_state(1:nr, 2) = (wf_cg(2) * rwfcr (1:nr, 1, ist, ias)) / spr (1:nr, is)
            wf_cr_per_state(1:nr, 4) = (wf_cg(2) * rwfcr (1:nr, 2, ist, ias)) / spr (1:nr, is)
        end if 
    end subroutine

    !> Returns the core kinetic energy density for state `ist`, defined by the quantum numbers \(\kappa\) and \(m\).
    !> Given the relativistic Dirac wavefunction `wf_cr_per_state`, obtained through the subroutine  
    !> [[gen_wf_cr_per_state(subroutine)]], the kinetic energy density `ked_cr_per_state` is returned. 
    !>
    !> The following kinetic energy density is evaluated for both \(i = 1, 2\) corresponding to the major and minor component:
    !>
    !>  \[
    !>     \tau^{\alpha}_{i, L M} = \frac{1}{2} \left( \vec{C}_{l, l \pm \frac{1}{2}, m}\right)^* 
    !>      \vec{C}_{l, l \pm \frac{1}{2}, m} \sum_{\pm} \sum_{\pm'} \vec{P}^{\pm \pm'}_{L M; l, m} 
    !>      u^{\alpha \pm^{*}}_{i, \kappa} u^{\alpha \pm'}_{i, \kappa},
    !> \]
    !> where \(\vec{C}_{l, l \pm \frac{1}{2}, m}\) are the cofficients of the Dirac wavefunctions as defined in 
    !> [[get_diracwf_cg_coeffs(function)]] and \(P^{\pm \pm'}_{L M; l, m}\) is given by 
    !> \[
    !> P^{\pm \pm'}_{L M; l, m} = 
    !>      \sum_{\substack{m'_l = \\ - (l \pm 1)}}^{l \pm 1}   \sum_{\substack{m''_l = \\ - (l \pm 1)}}^{l \pm 1} \sum_{\gamma} 
    !>      \left( \vec{C}^{l, m}_{l \pm 1, m'_{l}, \gamma} \right)^{*} \vec{C}^{l, m}_{l \pm' 1, m''_{l}, \gamma} 
    !>      G^{L M}_{l \pm 1, m'_{l}, l \pm' 1, m''_{l}} \hat{\mathbf{e}}_{\gamma}.
    !> \]
    !> Here \(\vec{C}_{l \pm 1, m'_{l}, \gamma}^{l, m}\) is a vector consisting of Clebsch-Gordan coefficients: 
    !> \[
    !> \vec{C}_{l \pm 1, m'_{l}, \gamma}^{l, m} =  
    !>      \begin{pmatrix}
    !>          C^{l, m-\frac{1}{2}}_{l \pm 1, m_{l}, 1, \gamma}  \\
    !>          C^{l, m+\frac{1}{2}}_{l \pm 1, m_{l}, 1, \gamma} 
    !>      \end{pmatrix},
    !> \] 
    !> which is evaluated through the subroutine [[get_clebsch_gordan_coeffs(subroutine)]].
    !>
    !> \(G^{L M}_{l \pm 1, m'_{l}, l \pm' 1, m''_{l}}\) are the Gaunt cofficients. 
    !>
    !> The radial functions \(u^{\alpha \pm}_{i, \kappa}\) are obtained through the function: [[get_grad_rad_fun(function)]] 
    !>
    !> If \(|m-\frac{1}{2}| <= l\), then the upper part of the spin spherical harmonics is evaluated for the major and  
    !> minor component of the Dirac wavefunction:
    !>
    !> If \(|m+\frac{1}{2}| <= l\), then the lower part of the spin spherical harmonics is evaluated for the major and  
    !> minor component of the Dirac wavefunction.
    !>
    !> In practice, the \(m\) index is provided through `m`, which is defined as \(m - \frac{1}{2}\).  
    subroutine gen_ked_cr_per_state(is, ia, ist, m, wf_cr_per_state, ked_cr_per_state)
        !> species index
        integer, intent(in) :: is
        !> atom index
        integer, intent(in) :: ia
        !> core state
        integer, intent(in) :: ist
        !> magnetic quantum number (passed in as m-1/2)
        integer, intent(in) :: m
        !> core wavefunction 
        real(dp), intent(in) :: wf_cr_per_state( :, :)
        !> kinetic energy density for core state
        complex(dp), intent(inout) :: ked_cr_per_state(:,:,:)

        integer :: l 
        l = spl(ist, is)

        ked_cr_per_state = zzero   
        if (abs(m) <= l) then
            call gen_ked_cr_per_comp(is, nrmt(is), spr(:, is), wf_cr_per_state( :, 1), ked_cr_per_state(:, :, 1), l, m)
            call gen_ked_cr_per_comp(is, nrmt(is), spr(:, is), wf_cr_per_state( :, 3), ked_cr_per_state(:, :, 3), l, m)
        end if
        if (abs(m+1) <= l) then
            call gen_ked_cr_per_comp(is, nrmt(is), spr(:, is), wf_cr_per_state( :, 2), ked_cr_per_state(: ,:, 2), l, m+1)
            call gen_ked_cr_per_comp(is, nrmt(is), spr(:, is), wf_cr_per_state( :, 4), ked_cr_per_state(:, :, 4), l, m+1)
        end if 
    end subroutine 

    !> Returns the kinetic energy density for one component in spherical harmonics representation.
    !> Given the relativistic Dirac wavefunction `wf_per_comp` for one component, obtained through the routine  
    !> [[gen_wf_cr_per_state(subroutine)]], the kinetic energy density `ked_cr_per_comp`for one component is returned.
    subroutine gen_ked_cr_per_comp(is, nr, r, wf_cr_per_comp, ked_cr_per_comp, current_l, current_m) 
        use gaunt
        !> species index
        integer, intent(in) :: is 
        !> number of radial grid points 
        integer, intent(in) :: nr
        !> radial mesh array in the MT-region
        real(dp), intent(in) :: r(nr)
        !> cofficients of the wavefunction in spherical harmonics representation
        real(dp), intent(in) :: wf_cr_per_comp(nr)
        !> kinetic energy density per state
        complex(dp), intent(out) :: ked_cr_per_comp(:, :)
        !> angular quantum number
        integer, intent(in) :: current_l
        !> magnetic quantum number 
        integer, intent(in) :: current_m
        
        ! local variables
        integer :: ir, l1, l2, m1, m2, lm1, lm2, i, k, lm_cr, l, m, lm, grad_array(2), l1_pm, l2_pm
        complex(dp) :: cg1, cg2, cg3
        complex(dp), allocatable :: tmp(:), res(:,:), gzfmt1(:,:,:), gzfmt2(:,:,:), zt1(:), zt2(:)
        real(dp), allocatable :: f1(:), f2(:), f3(:), rzfmt(:)

        allocate(tmp(nr), res(nr,ked_lmmaxvr), f3(nr))
        allocate(gzfmt1(ked_lmmaxvr, nrmtmax, 3), source = zzero)
        allocate(gzfmt2(ked_lmmaxvr, nrmtmax, 3), source = zzero)
        allocate(zt1(nrmtmax), zt2(nrmtmax), source = zzero)
        allocate(rzfmt(ked_lmmaxvr), source = 0.0_dp)

        ked_cr_per_comp = zzero
        grad_array = [-1, 1]
        l = current_l
        m = current_m
        lm = idxlm(l, m)

        do l1_pm = 1, 2 
            l1 = l + grad_array(l1_pm)
            If ((l1 >= 0) .and. (l1 <= ked_lmaxvr)) then
                do l2_pm = 1, 2
                    l2 = l + grad_array(l2_pm)
                    if ((l2 >= 0) .and. (l2 <= ked_lmaxvr)) then 
                        f1 = get_grad_rad_fun(r, nr, l, grad_array(l1_pm), wf_cr_per_comp(1:nr))
                        f2 = get_grad_rad_fun(r, nr, l, grad_array(l2_pm), wf_cr_per_comp(1:nr))

                        do m1 = - l1, l1
                            lm1 = idxlm(l1, m1)
                            do m2 = -l2, l2
                                lm2 = idxlm(l2, m2)
                                do i = 1, 3 
                                    call get_clebsch_gordan_coeffs(l, m, l1, m1, i, cg1)
                                    call get_clebsch_gordan_coeffs(l, m, l2, m2, i, cg2)
                                    
                                    zt1 = cmplx(f1, 0.0_dp, dp)
                                    zt2 = cmplx(f2, 0.0_dp, dp)

                                    gzfmt1(lm1, 1:nr, i) = cg1 * zt1(1:nr)
                                    gzfmt2(lm2, 1:nr, i) = cg2 * zt2(1:nr)

                                    tmp = 0.5_dp * conjg(gzfmt2( lm2, 1:nr, i)) * gzfmt1(lm1, 1:nr, i)
                                    do k = 1, gaunt_coeff_yry%num(lm1, lm2)
                                        lm_cr = gaunt_coeff_yry%lm2(k, lm1, lm2)
                                        if (lm_cr > ked_lmmaxvr) exit
                                        call zaxpy( nr, gaunt_coeff_yry%val(k, lm1, lm2), tmp, 1, ked_cr_per_comp(lm_cr, 1:nr), 1 )
                                    end do 
                                end do
                            end do 
                        end do
                
                    end if 
                end do 
            end if 
        end do 

    end subroutine
    
    !> Evaluates the Clebsch-Gordan coefficients \(C^{l, m}_{l \pm 1, m'_{l}, 1, \gamma}\), where 
    !> where \(l\) is given by `l`, \(m\) by `m`, \(l \pm 1\) by `ll`, \(m'_{l}, 1\) by `mm`. 
    !> Then returns the corresponding Cartesian-like representation for the specified direction \(i = 1, 2\) or \(3\).
    !>
    !> The Cartesian-like components \(\hat{\mathbf{e}}_{i}\) are derived from the spherical basis vectors 
    !> \(\hat{\mathbf{e}}_{\gamma}\) with \(\gamma = -1, 0, 1\), defined as:
    !> \[
    !>  \hat{\mathbf{e}}_{x} = - \frac{1}{\sqrt{2}} (\hat{\mathbf{e}}_{1} + \hat{\mathbf{e}}_{-1}),
    !>  \hat{\mathbf{e}}_{y} = \frac{i}{\sqrt{2}} (\hat{\mathbf{e}}_{1} - \hat{\mathbf{e}}_{-1}),
    !>  \hat{\mathbf{e}}_{z} = \hat{\mathbf{e}}_{0}.
    !> \]
    subroutine get_clebsch_gordan_coeffs(l, m, ll, mm, i, cg)
        use wigner3j_symbol, only: clebsch_gordan
        use constants, only: sqrt_two, zi
        !> angular quantum number 
        integer, intent(in) :: l
        !> magnetic quantum number 
        integer, intent(in) :: m
        !> angular quantum number with \(l \pm 1\)
        integer, intent(in) :: ll
        !> magnetic quantum number corresponding to the angular quantum number `ll` 
        integer, intent(in) :: mm
        !> direction of gradient (possible values are \(i = 1, 2, 3\)
        integer, intent(in) :: i
        !> clebsch-gordan coefficient
        complex(dp), intent(inout) :: cg

        real(dp) :: t1, t2, t3, g

        t1 = clebsch_gordan(ll, 1, l, mm, -1, m)
        t2 = clebsch_gordan(ll, 1, l, mm,  0, m)
        t3 = clebsch_gordan(ll, 1, l, mm,  1, m)

        if (i == 1) then 
            g = t1 - t3
            cg = cmplx( g / sqrt_two, 0, dp )
        elseif (i == 2) then 
            g = t1 + t3
            cg = cmplx( 0, - g / sqrt_two, dp )
        elseif (i == 3) then 
            cg = cmplx( t2, 0, dp )
        end if 

    end subroutine 

    !> Returns the Clebsch-Gordan coefficients of the Dirac wavefunction.
    !> If \(|\kappa| = l + 1\) they are given by: 
    !> \[
    !> \vec{C}_{l, l+\frac{1}{2}, m} = 
    !>    \begin{pmatrix}
    !>       \sqrt{\frac{l + m + \frac{1}{2}}{2l+1}}  \\
    !>      \sqrt{\frac{l - m + \frac{1}{2}}{2l+1}} 
    !>    \end{pmatrix}.
    !> \]
    !> If \(|\kappa| = l \)
    !> \[
    !> \vec{C}_{l,l-\frac{1}{2}, m} = 
    !>    \begin{pmatrix}
    !>       - \sqrt{\frac{l - m + \frac{1}{2}}{2l+1}}  \\
    !>      \sqrt{\frac{l + m + \frac{1}{2}}{2l+1}} 
    !>    \end{pmatrix},
    !> \]
    !>
    !> In practice the \(m\) index is given by `m`, which is defined as \(m - \frac{1}{2}\).  
    function get_diracwf_cg_coeffs(l, k, m) result(cg)
        !> angular quantum number
        integer, intent(in) :: l
        !> quantum number \(\kappa\)
        integer, intent(in) :: k
        !> magnetic quantum number (passed in as m-1/2)
        integer, intent(in) :: m
        !> clebsch-gordan cofficients for \(s=\pm \frac{1}{2}\)
        real(dp), allocatable :: cg(:)

        if (allocated(cg)) deallocate(cg)
        allocate(cg(2), source=0.0_dp)

        if (k == l+1) then
            cg(1) = sqrt(dble(l+m+1)/dble(2*l+1))
            cg(2) = sqrt(dble(l-m)/dble(2*l+1))
        else if (k == l) then
            cg(1) = - sqrt(dble(l-m)/dble(2*l+1))
            cg(2) = sqrt(dble(l+m+1)/dble(2*l+1))
        end if
    end function

    !> Returns the input radial function \(u^\alpha_{\kappa}(r)\) with a prefactor based on the specified gradient component.
    !> The input radial function \(u^\alpha_{\kappa}(r)\) is modified as follows:
    !>
    !> If `grad=-1`:
    !> \[ u^{\alpha,-}_{\kappa}(r) = \left[\frac{l}{2l+1}\right]^\frac{1}{2} \left[ \frac{l+1}{r} + \frac{{\rm d}}{{\rm d}r} \right] 
    !>    u^\alpha_{\kappa}(r). \]
    !>
    !> If `grad=1`: 
    !> \[ u^{\alpha,+}_{\kappa}(r) = \left[\frac{l+1}{2l+1}\right]^\frac{1}{2} \left[ \frac{l}{r} - \frac{{\rm d}}{{\rm d}r} \right] 
    !>    u^\alpha_{\kappa}(r). \]
    function get_grad_rad_fun(r, nr, l, grad, rad_fun) result(grad_rad_fun)
        !> radial mesh array in the MT-region
        real(dp), intent(in) :: r(:) 
        !> number of radial grid points 
        integer, intent(in) ::  nr
        !> angular quantum number
        integer, intent(in) :: l
        !> gradient component (-1 or 1) for radial functions
        integer, intent(in) :: grad 
        !> radial function
        real(dp), intent(in) :: rad_fun(:)
        !> modified radial function  
        real(dp), allocatable :: grad_rad_fun(:)

        real(dp) :: lfac
        real(dp), allocatable :: deriv_rad_func(:), cf(:,:)
        
        if( allocated( grad_rad_fun ) ) deallocate( grad_rad_fun )
        allocate( grad_rad_fun(nr), source=0.0_dp )

        if( grad < 0 ) then
            lfac = sqrt( dble( l ) / dble( 2 * l + 1 ) )
        else
            lfac = sqrt( dble( l + 1 ) / dble( 2 * l + 1 ) )
        end if

        allocate(deriv_rad_func(nr), cf(3, nr), source=0.0_dp)
        call fderiv(1, nr, r(1:nr), rad_fun(1:nr), deriv_rad_func(1:nr), cf)

        if( grad < 0 ) then
            grad_rad_fun(1:nr) = lfac * ( dble(l + 1) * rad_fun(1:nr) / r(1:nr)  + deriv_rad_func(1:nr) )
        else
            grad_rad_fun(1:nr) = lfac * (dble(l) * rad_fun(1:nr) /  r(1:nr)  - deriv_rad_func(1:nr) )
        end if
        deallocate(deriv_rad_func, cf)
    end function 

    !> Constructs and returns a formatted error message indicating a mismatch
    !> in quantum numbers or indices for wavefunction reconstruction. 
    subroutine wf_cr_error_msg(error_message, l, k, m, is, ia, ist)
        !> error message
        character(:), allocatable, intent(out) :: error_message
        !> angular quantum number
        integer, intent(in) :: l
        !> quantum number \(\kappa\)
        integer, intent(in) :: k 
        !> magnetic quantum number (passed in as m-1/2)
        integer, intent(in) :: m 
        !> species index
        integer, intent(in) :: is
        !> atom index
        integer, intent(in) :: ia
        !> core state
        integer, intent(in) :: ist

        character(len=3) :: l_str, k_str, m_str, is_str, ia_str, ist_str
        character(len=500) :: buffer_string

        write(l_str, '(I0)') l
        write(k_str, '(I0)') k
        write(m_str, '(I0)') m
        write(is_str, '(I0)') is
        write(ia_str, '(I0)') ia
        write(ist_str, '(I0)') ist

        write(buffer_string, '(A, A, A, A, A, A, A, A, A, A, A, A, A, A)') &
            "Error(wavefcr): mismatched l, k or m : ", trim(l_str), ", ", trim(k_str), ", ", trim(m_str), ", for species ", trim(is_str), ", atom ", trim(ia_str), " and state ", trim(ist_str)

        error_message = trim(buffer_string)
    end subroutine 

end module 