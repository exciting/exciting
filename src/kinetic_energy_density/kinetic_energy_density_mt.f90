!> Module provides subroutines to calculate the kinetic energy density in the muffin-tin region for the valence electrons `ked_mt`. 
module kinetic_energy_density_mt
    use kinetic_energy_density_vars
    use mod_atoms, only: nspecies, natoms, idxas, natmtot         
#include "asserts.fpp"
    use precision, only: dp
    use constants, only: zzero, zone
    use mod_APW_LO, only: nlotot, nlorb
    use mod_eigensystem, only: idxlo
    use mod_spin, only: ncmag
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use general_matrix_multiplication, only: matrix_multiply

    implicit none 
    private

    public :: gen_denmat_k,& 
              gen_ked_mt

    interface gen_denmat_k
        procedure gen_denmat_k_spin_unpolarised
        procedure gen_denmat_k_spin_polarised
    end interface

    interface gen_ked_mt
        procedure gen_ked_mt_polarised
        procedure gen_ked_mt_unpolarised
    end interface

    contains 

    !> Main routine to calculate the spin-polarised kinetic energy density `ked_mt`
    !> and to calculate the spin kinetic energy density `ked_magmt`.
    subroutine gen_ked_mt_polarised(ked_mt, ked_magmt, ked_mat_alpha, ked_mat_beta, ked_mat_ab)
        !> scalar kinetic energy density
        real(dp), intent(inout) :: ked_mt(:, :, :)
        !> spin kinetic energy density  
        real(dp), intent(inout) :: ked_magmt(:, :, :, :)
        !> diagonal components of the density matrix
        complex(dp), intent(in) :: ked_mat_alpha(:, :, :), ked_mat_beta(:, :, :)
        !> off-diagonal components of the density matrix 
        complex(dp), optional, intent(in) :: ked_mat_ab(:, :, :)

        complex(dp), allocatable:: ked_mat_alpha_plus_beta(:, :, :), ked_mat_alpha_min_beta(:, :, :)
        
        if( allocated( ked_mat_alpha_plus_beta ) ) deallocate( ked_mat_alpha_plus_beta )
        allocate(ked_mat_alpha_plus_beta(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max, natmtot), source=zzero)

        if( allocated( ked_mat_alpha_min_beta ) ) deallocate( ked_mat_alpha_min_beta )
        allocate(ked_mat_alpha_min_beta(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max, natmtot), source=zzero)

        ! Calculate the kinetic energy density: ked_mt
        ked_mat_alpha_plus_beta = ked_mat_alpha + ked_mat_beta 
        call gen_ked_or_ked_magmt(zone, ked_mat_alpha_plus_beta, zzero, ked_mt) 

        ! In the spin polarised case, calculate: ked_magmt
        ked_mat_alpha_min_beta = ked_mat_alpha - ked_mat_beta
        if (ncmag) then
            CALL_ASSERT( present(ked_mat_ab), message='if ncmag=true the array ked_mat_ab needs to be given.')
            ! noncollinear case
            call gen_ked_or_ked_magmt(zone, ked_mat_alpha_min_beta, zzero, ked_magmt(:, :, :, 3))
            call gen_ked_or_ked_magmt(cmplx(2.0, 0.0, dp), ked_mat_ab, zzero, ked_magmt(:, :, :, 1))
            call gen_ked_or_ked_magmt(cmplx(0.0, 2.0, dp), ked_mat_ab, zzero, ked_magmt(:, :, :, 2))
        else
            ! collinear case
            call gen_ked_or_ked_magmt(zone, ked_mat_alpha_min_beta, zzero, ked_magmt(:, :, :, 1))
        end if
    end subroutine 

    !> Main routine to calculate the spin-unpolarised kinetic energy density `ked_mt`.
    subroutine gen_ked_mt_unpolarised(ked_mt, ked_mat_alpha)
        !> kinetic energy density
        real(dp), intent(inout) :: ked_mt(:, :, :)
        !> component of the density matrix
        complex(dp), intent(in) :: ked_mat_alpha(:, :, :)

        call gen_ked_or_ked_magmt(zone, ked_mat_alpha, zzero, ked_mt) 
    end subroutine 

    !> Calculates the kinetic energy density (or spin kinetic energy density) in the muffin-tin region 
    !> for the valence electrons for all species. 
    !> `ked_mt` is returned in spherical harmonics. 
    !> For more details see [[gen_ked_mt_ias(subroutine)]].
    subroutine gen_ked_or_ked_magmt(alpha, ked_mat, beta, ked_mt)
        !> prefactor \(a\)
        complex(dp), intent(in) :: alpha
        !> density matrix for the kinetic energy density
        complex(dp), intent(in) :: ked_mat(:, :, :)
        !> prefactor \(b\)
        complex(dp), intent(in) :: beta
        !> kinetic energy density (could also be u-array - to be changed)
        real(dp), intent(inout) :: ked_mt(:, :, :) 
                
        integer :: is, ia, ias

        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia, is)
                call gen_ked_mt_ias(is, ias, alpha, ked_mat(:, :, ias), beta, ked_mt(:, :, ias)) 
            end do 
        end do 
    end subroutine 

    !> Returns the kinetic energy density for the muffin-tin region and valence electrons `ked_mt_ias`
    !> for given atom \(\alpha\).
    !> This routine computes, for every species \(\alpha\), kinetic energy density for the muffin-tin region and valence electrons `ked_mt_ias`
    !> using the following formula: 
    !>  \[ \tau^{\alpha,{\rm out}}_{\rm lm} = b\, \tau^{\alpha,{\rm in}}_{\rm lm} + 
    !>  a\, \sum_{\lambda_1 \lambda_2} D_{\lambda_1 \lambda_2}^{\alpha} \sum_{\pm} \sum_{\pm'} P^{\pm \pm'}_{lm; l_1 m_1, l_2 m_2} 
    !>      g_{\lambda_1}^{\alpha \pm} g_{\lambda_2}^{\alpha \pm'} \]
    !> where \(\lambda\) counts all muffin-tin basis functions and is a combined \((l,m,\xi)\) index for (L)APWs and \(L\) for all LOs. For more details on the 
    !> construction of the basis and indexing, see [[muffin_tin_basis(module)]].
    !>
    !> \(D_{\lambda \lambda^{'}}^{\alpha}\) is the density matrix `ked_mat` for the species \(\alpha\). 
    !> For a detailed derivation of the calculation of \(D_{\lambda_1 \lambda_2}^{\alpha}\), see routine: [[gen_denmat_k(subroutine)]].
    !>
    !> \(P^{\pm \pm'}_{lm; l_1 m_1, l_2 m_2}\) is a product of Clebsch-Gordan and Gaunt coefficients, see routine: [[clebsch_gordan_gaunt_product(subroutine)]].
    !>
    !> \(g_{\lambda_1}^{\alpha \pm'} g_{\lambda_2}^{\alpha \pm'}\) is a radial product of the basis functions, see routine: [[radial_product(subroutine)]].
    subroutine gen_ked_mt_ias(is, ias, alpha, ked_mat_ias, beta, ked_mt_ias) 
        use mod_muffin_tin, only: idxlm
        !> index of the species of the MT
        integer, intent(in) :: is
        !> index of the atom (within the species) of the MT
        integer, intent(in) :: ias
        !> prefactor \(a\)
        complex(dp), intent(in) :: alpha
        !> density matrix for the kinetic energy density
        complex(dp), intent(in) :: ked_mat_ias(:, :)
        !> prefactor \(a\)
        complex(dp), intent(in) :: beta
        !> kinetic energy density in the MT region
        real(dp), intent(inout) :: ked_mt_ias(:, :) 

        ! radial function product
        real(dp), allocatable :: radprod(:)
        ! product of clebsch-gordan and gaunt
        real(dp), allocatable :: cg_gaunt_prod(:)
        ! temporary array of kinetic energy density
        real(dp), allocatable :: ked_temp(:, :)
        
        ! loop variables
        integer :: l1, l2, lam1, lam2, m1, m2, lm1, lm2, idx1, idx2, lm, nr
        integer :: l1_pm, l2_pm, llstart1, llstart2, llstop1, llstop2, ll1, ll2
        integer :: grad_array(2), gradient_dir
        complex(dp) :: ked_mat_val

        if( beta == zzero ) then
           ked_mt_ias = zzero 
        else if( beta /= zone ) then
           ked_mt_ias = beta * ked_mt_ias
        end if
        if( alpha == zzero ) return

        nr = ked_mt_basis%n_rad_grid(is)
        grad_array = [-1, 1]

        allocate( cg_gaunt_prod(ked_lmmaxvr) )
        allocate( ked_temp(nr, ked_lmmaxvr), source=0.0_dp )
        allocate( radprod(nr) )

        do l1_pm = 1, 2 
            llstart1 = max(0, grad_array(l1_pm))
            llstop1 = ked_lmaxapw + grad_array(l1_pm) 
            do l2_pm = 1, 2
                llstart2 = max(0, grad_array(l2_pm))
                llstop2 = ked_lmaxapw + grad_array(l2_pm)
!$omp parallel default( shared ) private( ll1, ll2, l1, l2, m1, m2, lm, lm1, lm2, lam1, lam2, idx1, idx2, gradient_dir, radprod, ked_mat_val, cg_gaunt_prod ) reduction( +:ked_temp )
!$omp do collapse(2)
                do ll1 = llstart1, llstop1
                    do ll2 = llstart2, llstop2
                        l1 = ll1 - grad_array(l1_pm)
                        l2 = ll2 - grad_array(l2_pm)                
                        do lam1 = 1, ked_mt_basis%n_rad_fun(l1, is)
                            do lam2 = 1, ked_mt_basis%n_rad_fun(l2, is)
                                call radial_product(is, ias, l1, l2, lam1, lam2, grad_array(l1_pm), grad_array(l2_pm), radprod)                             

                                cg_gaunt_prod = 0.0_dp
                                do m1 = -l1, l1 
                                    lm1 = idxlm(l1, m1)
                                    idx1 = ked_mt_basis%idx_basis_fun(lm1, lam1, is)
                                    do m2 = -l2, l2
                                        lm2 = idxlm(l2, m2)
                                        idx2 = ked_mt_basis%idx_basis_fun(lm2, lam2, is)

                                        ked_mat_val = ked_mat_ias(idx1, idx2)
                                        if( alpha /= zone ) ked_mat_val = alpha * ked_mat_val 
                                        do gradient_dir = 1, 3
                                            call clebsch_gordan_gaunt_product(lm1, lm2, grad_array(l1_pm), grad_array(l2_pm), gradient_dir, ked_mat_val, cg_gaunt_prod) 
                                        end do 
                                    end do 
                                end do 

                                do lm = 1, ked_lmmaxvr
                                   if (cg_gaunt_prod(lm) .ne. 0.0_dp) then
                                        ked_temp(1:nr, lm) = ked_temp(1:nr, lm) +  cg_gaunt_prod(lm) * radprod(1:nr) 
                                   end if 
                                end do
                            end do 
                        end do 
                     end do
                 end do 
!$omp end do
!$omp end parallel
            end do 
        end do 
         
        ked_mt_ias(1:ked_lmmaxvr, 1:nr) = transpose(ked_temp(1:nr, 1:ked_lmmaxvr))

        deallocate(ked_temp, cg_gaunt_prod, radprod)    
    end subroutine 

    !> Returns the sum of Clebsch-Gordan coefficient times Gaunt coefficients with a given pre-factor `alpha`
    !> for a given `l1`, `m1` and `l2`, `m2`, i.e., 
    !> \[ 
    !>  P^{\pm \pm', i}_{l m; l_{\lambda_1} m_{\lambda_1} l_{\lambda_2} m_{\lambda_2}} = \alpha
    !>  \sum_{\substack{m' = \\ -(l_{\lambda_1} \pm 1)}}^{l_{\lambda_1} \pm 1} \sum_{\substack{m'' = \\ -(l_{\lambda_2} \pm' 1)}}^{l_{\lambda_2} \pm' 1} \sum_{i} 
    !>   \left( C^{l_{\lambda_1}, m_{\lambda_1}}_{l_{\lambda_1} \pm 1, m^{'}, 1, i}\right)^{*} C^{l_{\lambda_2}, m_{\lambda_2}}_{l_{\lambda_2} \pm' 1, m^{\prime \prime}, 1, i} 
    !>  G^{l m}_{l_{\lambda_1} \pm 1, m^{ \prime}, l_{\lambda_2} \pm' 1, m^{ \prime  \prime}}.
    !> \]
    !> for a specified direction \(i = 1, 2, 3\).
    subroutine clebsch_gordan_gaunt_product(lm1, lm2, lgrad, rgrad, gradient_dir, alpha, cg_gaunt_prod) 
        use gaunt
        !> Combined index of \(l_1,m_1\)
        integer, intent(in) :: lm1
        !> Combined index of \(l_2,m_2\)
        integer, intent(in) :: lm2
        !> gradient component (-1 or 1) for left basis functions
        integer, intent(in) :: lgrad
        !> gradient component (-1 or 1) for right basis functions
        integer, intent(in) :: rgrad
        !> Direction of gradient (1, 2 or 3)
        integer, intent(in):: gradient_dir
        !> pre-factor 
        complex(dp), intent(in) :: alpha
        !> Product of Clebsch-Gordan and Gaunt coefficients
        real(dp), intent(inout) :: cg_gaunt_prod(:)
        
        integer :: llmm1, llmm2, lm, i, j, k
        complex(dp) :: cg1, cg2, cg3
        type(non_zero_gaunt_complex), pointer :: gntz

        do i = 1, ked_cg_num(gradient_dir, lgrad, lm1)
            llmm1 = ked_cg_lm(i, gradient_dir, lgrad, lm1)
            cg1 = ked_cg_val(i, gradient_dir, lgrad, lm1) 
            do j = 1, ked_cg_num(gradient_dir, rgrad, lm2)
                llmm2 = ked_cg_lm(j, gradient_dir, rgrad, lm2)
                cg2 = ked_cg_val(j, gradient_dir, rgrad, lm2) 
                cg3 = conjg(cg1) * cg2
                gntz => gaunt_coeff_yry
                do k = 1, gntz%num(llmm1, llmm2) 
                    lm = gntz%lm2(k, llmm1, llmm2)
                    if( lm > ked_lmmaxvr) exit
                    cg_gaunt_prod(lm) = cg_gaunt_prod(lm) + dble(cg3 *  gntz%val(k, llmm1, llmm2) * alpha)
                end do
            end do 
        end do 
    end subroutine 


    !> Computes product of two radial MT basis functions with angular momentum
    !> \(l_1\) and \(l_2\), respectively, 
    !> \[ g_{\tilde{\lambda}_1}(r) \, g_{\tilde{\lambda}_2}(r) \;. \]
    !> Depending on the input of `lgrad` and `rgrad`, will calculate the product of the following radial basis functions:
    !> When `lgrad=-1` or `rgrad=-1`: 
    !> \[ g^{\alpha,-}_{\tilde{\lambda}}(r) = \left[\frac{l}{2l+1}\right]^\frac{1}{2} \left[ \frac{l+1}{r} + \frac{{\rm d}}{{\rm d}r} \right] 
    !>    g^\alpha_{\tilde{\lambda}}(r) \]
    !> When `lgrad=1` or `rgrad=1`: 
    !> \[ g^{\alpha,+}_{\tilde{\lambda}}(r) = \left[\frac{l+1}{2l+1}\right]^\frac{1}{2} \left[ \frac{l}{r} - \frac{{\rm d}}{{\rm d}r} \right] 
    !>    g^\alpha_{\tilde{\lambda}}(r) \]
    !> with \(g^\alpha_{\tilde{\lambda}}(r)\) being the bare basis function.
    subroutine radial_product(is, ias, l1, l2, lam1, lam2, lgrad, rgrad, f3)
        !> index of the species of the MT
        integer, intent(in) :: is
        !> index of the atom (within the species) of the MT
        integer, intent(in) :: ias
        !> angular momentum of left basis functions
        integer, intent(in) :: l1
        !> angular momentum of right basis functions
        integer, intent(in) :: l2
        !> \(\tilde{\lambda}\) on the left 
        integer, intent(in) :: lam1
        !> \(\tilde{\lambda}\) on the right 
        integer, intent(in) :: lam2
        !> gradient component (-1 or 1) for left basis functions // (`>0` for \(g^{\alpha,+}_\tilde{\lambda}\), `<0` for \(g^{\alpha,-}_\tilde{\lambda}\)\) )
        integer, intent(in) :: lgrad
        !> gradient component (-1 or 1) for right basis functions
        integer, intent(in) :: rgrad
        !> radial product
        real(dp), intent(inout) :: f3(:)

        integer ::  nr
        real(dp), allocatable :: f1(:), f2(:)
        
        nr = ked_mt_basis%n_rad_grid(is)
        
        f1 = ked_mt_basis%get_gradient_rad_fun(l1, is, ias, lam1, lgrad, radial_derivative=0)
        f2 = ked_mt_basis%get_gradient_rad_fun(l2, is, ias, lam2, rgrad, radial_derivative=0)

        f3(1:nr) =  f1(1:nr) * f2(1:nr)
    end subroutine

    !> Returns the matrices \(D^{\mathbf{k} \alpha'}_{\alpha \alpha}\), \(D^{\mathbf{k} \alpha'}_{\alpha \beta}\) and \(D^{\mathbf{k} \alpha'}_{\beta \beta}\) 
    !> of the total density matrix for the kinetic energy density for the \(i\)-th \(\mathbf{k}\)-point `ik` and for each species \(\alpha'\):
    !> \[ D^{\mathbf{k} \alpha'} =
    !>    \begin{pmatrix}
    !>    D^{\mathbf{k} \alpha'}_{\alpha \alpha} & D^{\mathbf{k} \alpha'}_{\beta \alpha}\\
    !>    D^{\mathbf{k} \alpha'}_{\alpha \beta} & D^{\mathbf{k} \alpha'}_{\beta \beta}
    !>     \end{pmatrix} \]
    !> where each matrix is constructed as follows: 
    !> \[
    !>   D^{\mathbf{k} \alpha'}_{\lambda_1 \lambda_2, \sigma \sigma'} = 
    !>      \sum_{n} w_{\mathbf{k}} f_{n \mathbf{k}} \left( \sum_{j} C^{n \mathbf{k} \alpha, \text{SV}}_{j \sigma} 
    !>      C^{j \mathbf{k} \alpha, \text{FV}}_{\lambda_1 \sigma} \right)^{*} \left( \sum_{i} C^{n \mathbf{k} \alpha, \text{SV}}_{i \sigma'} 
    !>      C^{i \mathbf{k} \alpha, \text{FV}}_{\lambda_2 \sigma'} \right).
    !> \] 
    !> where \(\sigma\) and \(\sigma'\) can correspond to \(\alpha\) or \(\beta\).
    subroutine gen_denmat_k_spin_polarised(ik, occsvk, evecfv_k, evecsv_k, apwalmk, ked_mat_alpha, ked_mat_beta, ked_mat_ab)
        !> i-th k-point
        integer, intent(in) :: ik
        !> occupation numbers at \({\bf k}\)
        real(dp), intent(in) :: occsvk(:)
        !> 1st var eigenvectors at \({\bf k}\)
        complex(dp), intent(in) :: evecfv_k(:, :)
        !> 2nd var eigenvectors at \({\bf k}\)
        complex(dp), intent(in) :: evecsv_k(:, :)
        !> (L)APW matching coefficients \(A^\alpha'_{{\bf G+k},lm,\xi}\) 
        complex(dp), intent(in) :: apwalmk(:, :, :, :)
        !> muffin-tin KED matrix alpha
        complex(dp), intent(inout) :: ked_mat_alpha(:, :, :)
        !> muffin-tin KED matrix beta
        complex(dp), intent(inout) :: ked_mat_beta(:, :, :)
        !> muffin-tin KED matrix alphabeta
        complex(dp), optional, intent(inout) :: ked_mat_ab(:, :, :)
        
        integer :: i, ist, nst, is, ia, ias
        logical :: lo_dependent(nspecies)

        complex(dp), allocatable :: evecfv_mt(:, :)
        complex(dp), allocatable :: wfalpha(:, :), wfalpha2(:, :)
        complex(dp), allocatable :: wfbeta(:, :), wfbeta2(:, :)
        complex(dp), allocatable :: ked_mat_ab_temp(:, :)

        allocate( ked_mat_ab_temp(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max), source = zzero )

        lo_dependent = [ (merge(.true., .false., nlorb(is) /= 0), is = 1, nspecies) ]

        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia, is)
                    ! generate MT eigenvector for first variational eigenvector at current k-point
                    call ked_mt_basis%transform_evec( is, ked_Gkset%ngk(1, ik), nlotot, idxlo(:, :, ias), apwalmk(:, :, :, ias), evecfv_k(:, :), evecfv_mt, lo_dependent(is))
                    ! Calculate ked_mat_alpha
                    call calculate_ked_mat(is, ik, occsvk, evecfv_mt, evecsv_k(1:nstfv, 1:nstsv), wfalpha, wfalpha2, ked_mat_alpha(:, :, ias))

                    ! Calculate ked_mat_beta
                    call calculate_ked_mat(is, ik, occsvk, evecfv_mt, evecsv_k(1+nstfv:nstsv, 1:nstsv), wfbeta, wfbeta2, ked_mat_beta(:, :, ias))

                    ! Calculate ked_mat_ab
                    if (ncmag) then
                        CALL_ASSERT( present(ked_mat_ab), "ked_mat_ab needs to be given.")
                        ked_mat_ab_temp = zzero
                        call matrix_multiply(wfalpha, wfbeta2(1:ked_mt_basis%n_basis_fun_max, 1:nstsv), ked_mat_ab_temp(:, :), 'n', 'c')
                        ked_mat_ab(:, :, ias) = ked_mat_ab(:, :, ias) + ked_mat_ab_temp(:, :)
                    end if 
            end do 
        end do

        deallocate( wfalpha, wfbeta, evecfv_mt, ked_mat_ab_temp)
    end subroutine 
    
    !> Only returns \(D^{\mathbf{k} \alpha'}_{\alpha \alpha}\), for more see routine: [[gen_denmat_k_spin_polarised(subroutine)]].
    subroutine gen_denmat_k_spin_unpolarised(ik, occsvk, evecfv_k, evecsv_k, apwalmk, ked_mat_alpha)
        !> i-th k-point
        integer, intent(in) :: ik
        !> occupation numbers at \({\bf k}\)
        real(dp), intent(in) :: occsvk(:)
        !> 1st var eigenvectors at \(k\)
        complex(dp), intent(in) :: evecfv_k(:, :)
        !> 2nd var eigenvectors at \(k\)
        complex(dp), intent(in) :: evecsv_k(:, :)
        !> (L)APW matching coefficients \(A^\alpha'_{{\bf G+k},lm,\xi}\) 
        complex(dp), intent(in) :: apwalmk(:, :, :, :)
        !> muffin-tin density response matrix alpha
        complex(dp), intent(inout) :: ked_mat_alpha(:, :, :)

        integer :: i, ist, nst
        integer :: is, ia, ias
        logical :: lo_dependent(nspecies)

        complex(dp), allocatable :: evecfv_mt(:, :)
        complex(dp), allocatable :: wfalpha(:, :), wfalpha2(:, :)

        lo_dependent = [ (merge(.true., .false., nlorb(is) /= 0), is = 1, nspecies) ]

        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia, is)
                call ked_mt_basis%transform_evec( is, ked_Gkset%ngk(1, ik), nlotot, idxlo(:, :, ias), apwalmk(:, :, :, ias), evecfv_k(:, :), evecfv_mt, lo_dependent(is))
                call calculate_ked_mat(is, ik, occsvk, evecfv_mt, evecsv_k(1:nstfv, 1:nstsv), wfalpha, wfalpha2, ked_mat_alpha(:, :, ias))
 
                deallocate( evecfv_mt)
            end do 
        end do              

        deallocate( wfalpha )
    end subroutine 

    !> Returns a component of the density matrix for the kinetic energy density for a given \(\mathbf{k}\)-point,
    !> \[
    !>   D^{\mathbf{k}}_{\lambda_1 \lambda_2} = 
    !>      \sum_{n} w_{\mathbf{k}} f_{n \mathbf{k}} \left( \sum_{j} C^{n \mathbf{k}, \text{SV}}_{j} 
    !>      C^{j \mathbf{k}, \text{FV}}_{\lambda_1} \right)^{*} \left( \sum_{i} C^{n \mathbf{k}, \text{SV}}_{i} 
    !>      C^{i \mathbf{k}, \text{FV}}_{\lambda_2} \right).
    !> \] 
    subroutine calculate_ked_mat(is, ik, occsvk, evecfv_mt, evecsv_k, wf, wf2, ked_mat)
        use general_matrix_multiplication, only: matrix_multiply
        !> i-th \({\bf k}\)-point
        integer, intent(in) :: ik
        !> index of the species of the MT
        integer, intent(in) :: is 
        !> occupation numbers at \({\bf k}\)
        real(dp), intent(in) :: occsvk(:)
        !> first-variational eigenvector \(C^{j \mathbf{k}, \text{FV}}_{\lambda_1}\) in muffin-tin basis
        complex(dp), intent(in) :: evecfv_mt(:, :)
        !> second-variational eigenvector \(C^{n \mathbf{k}, \text{SV}}_{j}\)
        complex(dp), intent(in) :: evecsv_k(:, :)
        !> Product of the first- and second-variational eigenvectors: \(C^{j \mathbf{k}, \text{FV}}_{\lambda_1} C^{n \mathbf{k}, \text{SV}}_{j}\)
        complex(dp), allocatable, intent(out) :: wf(:, :)
        !> Product of the first- and second-variational eigenvectors and \(w_{\mathbf{k}} f_{n \mathbf{k}}\)
        complex(dp), allocatable, intent(out) :: wf2(:, :)
        !> muffin-tin density response matrix alpha
        complex(dp), intent(inout) :: ked_mat(:, :)

        complex(dp), allocatable :: ked_mat_temp(:, :)
        real(dp), allocatable :: wgt(:)
        integer, allocatable :: bands(:)
        integer :: i, nst, ist

        if (allocated(wf)) deallocate(wf) 
        allocate( wf(ked_mt_basis%n_basis_fun(is), nstsv), source = zzero )
        allocate( ked_mat_temp(ked_mt_basis%n_basis_fun(is), ked_mt_basis%n_basis_fun(is)), source = zzero)

        CALL_ASSERT(size(evecsv_k,1) == nstfv, message='Dim 1 of second-variational eigenvector has to be equal to nstfv.')
        
        allocate( wgt(nstsv) )
        wgt = [(ked_kset%wkpt(ik)*occsvk(i), i=1, nstsv)]
        bands = pack( [(i, i=1, nstsv)], &
                      [(wgt(i) == wgt(i), i=1, nstsv)] )
        nst = size( bands )
        
        ! Calculate ked_mat 
        call matrix_multiply(evecfv_mt, evecsv_k(1:nstfv, 1:nstsv), wf)

        if (allocated(wf2)) deallocate(wf2) 
        allocate(wf2(ked_mt_basis%n_basis_fun(is), nstsv), source=zzero)
        do i = 1, nstsv
            wf2(:, i) = 0.5_dp * ked_kset%wkpt(ik)*occsvk(i) * wf(:, i) 
        end do

        ked_mat_temp = zzero
        call matrix_multiply(wf, wf2, ked_mat_temp, 'n', 'c')
        ked_mat(1:ked_mt_basis%n_basis_fun(is), 1:ked_mt_basis%n_basis_fun(is)) = ked_mat(1:ked_mt_basis%n_basis_fun(is), 1:ked_mt_basis%n_basis_fun(is)) + ked_mat_temp
    end subroutine

end module 
