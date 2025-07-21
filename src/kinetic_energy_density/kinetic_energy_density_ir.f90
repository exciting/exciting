!> Module provides subroutines to calculate the kinetic energy density `ked_ir` for the interstitial region.
module kinetic_energy_density_ir
    use precision, only: dp
    use mod_spin, only : nspinor
    use modinput, only: input
    use mod_eigenvalue_occupancy, only: nstsv, nstfv
    use mod_lattice, only: omega
    use mod_eigenvalue_occupancy, only: occsv
    use constants, only: zzero, zi
    use kinetic_energy_density_vars, only: ked_Gkset, ked_kset, ked_Gset
    use general_matrix_multiplication, only: matrix_multiply
    use m_zfftifc, only: zfftifc

    implicit none 
    private
    public :: gen_ked_ir
    
    interface gen_ked_ir
        procedure gen_ked_ir_polarised
        procedure gen_ked_ir_unpolarised
    end interface
    
    contains

        !> Given the first-variational eigenvector `evecfv` at the current k-point `ik`, 
        !> calculates the kinetic energy density \( \tau(\mathbf{r}) \) for the interstitial region `kinetic_energy_density_ir`.
        !> 
        !> First, the gradient of the wavefunction is calculated in the reciprocal space. 
        !> For a detailed derivation see [[gen_grad_wavefunction_ir_unpolarised(subroutine)]]. 
        !>
        !> The returned gradient \(\nabla \phi_{\mathbf{G}j}^{\mathbf{k}}\) in reciprocal space is transformed to real space.
        !> The kinetic energy density in the spin-unpolarised case is then obtained in real space as:
        !> \[ \tau(\mathbf{r}) = \frac{1}{2} \sum_{j \mathbf{k}} \frac{w_{\mathbf{k}} f_{j\mathbf{k}}}{\Omega} 
        !>                                      \sum_{\mathbf{G} \mathbf{G}^{'}} 
        !>                                      ( \nabla \phi_{\mathbf{G}+\mathbf{k}}^{j}(\mathbf{r}) ) \cdot 
        !>                                      ( \nabla \phi_{\mathbf{G^{'}}+\mathbf{k}}^{j}(\mathbf{r}) )^{*} \]
        !> where \(\nabla \phi_{\mathbf{G}+\mathbf{k}}^{j}(\mathbf{r})\) denotes the real-space gradient obtained 
        !> after transforming from reciprocal space. The second-variational occupation number is \(f_{j\mathbf{k}}\), 
        !> \(w_{\mathbf{k}}\) are the k-point weights, and \(\Omega\) is the unit cell volume.
        subroutine gen_ked_ir_unpolarised(ik, evecfv, ked_ir_k)
            !> i-th k-point
            integer, intent(in) :: ik
            !> first variational eigenvectors
            complex(dp), intent(in) :: evecfv (:, :) 
            !> k-dependent kinetic energy density in the interstitial region
            real(dp), intent(inout) :: ked_ir_k(:)

            real(dp) :: wo, w
            integer :: jst, grad_dir
            complex(dp), allocatable :: zfft_grad(:, :), zfft(:)

            allocate( zfft_grad(ked_Gset%ngrtot, 3) )
            allocate( zfft(ked_Gset%ngrtot) )

            ked_ir_k = 0.0_dp
            do jst = 1, nstsv
                wo = ked_kset%wkpt(ik) * occsv(jst, ik)
                if (abs(wo) .gt. input%groundstate%epsocc) then
                    w = wo / omega

                    call gen_grad_wavefunction_ir_unpolarised(ik, evecfv(:, jst), zfft_grad)

                    do grad_dir = 1, 3
                        call zfftifc(3, ked_Gset%ngrid, 1, zfft_grad(:, grad_dir))
                        ked_ir_k(:) = ked_ir_k(:) + w * 0.5_dp * ( dble( zfft_grad(:, grad_dir) )**2 + aimag( zfft_grad(:, grad_dir) )**2 ) 
                    end do

                end if 
            end do 
        end subroutine 

        !> Given the first-variational eigenvector `evecfv` and second-variational eigenvector `evecsv` at the current k-point `ik`, 
        !> calculates the kinetic energy density \( \tau(\mathbf{r}) \) for the interstitial region `kinetic_energy_density_ir`.
        !> In the spin-polarised case the spin kinetic energy density vector `ked_magir` is also computed.
        !> 
        !> First, the gradients of the wavefunction are calculated in the reciprocal space. 
        !> For a detailed derivation see [[gen_grad_wavefunction_ir_polarised(subroutine)]]. 
        !>
        !> The returned gradients \(\nabla \phi_{\alpha,\mathbf{G}j}^{\mathbf{k}}\) and \(\nabla \phi_{\beta,\mathbf{G}j}^{\mathbf{k}}\)
        !> are then transformed to real space and from that the components of kinetic energy density matrix are calculated:
        !> \[ D = \begin{pmatrix}
        !>      \tau_{\alpha \alpha}(\mathbf{r}) & \tau_{\beta \alpha}(\mathbf{r}) \\
        !>      \tau_{\alpha \beta}(\mathbf{r}) & \tau_{\beta \beta}(\mathbf{r})
        !>    \end{pmatrix} \]
        !> where:
        !> \[ \tau_{\alpha \alpha}(\mathbf{r}) = \frac{1}{2} \sum_{\mathbf{G} \mathbf{G}^{'}} \sum_{j \mathbf{k}} \frac{w_{\mathbf{k}} f_{j\mathbf{k}}}{\Omega} 
        !>                                      ( \nabla \phi_{\alpha,\mathbf{G}+\mathbf{k}}^{j} (\mathbf{r}) ) \cdot 
        !>                                      ( \nabla \phi_{\alpha,\mathbf{G^{'}}+\mathbf{k}}^{j} (\mathbf{r}) )^{*} \]
        !> \[ \tau_{\beta \beta}(\mathbf{r}) = \frac{1}{2} \sum_{\mathbf{G} \mathbf{G}^{'}} \sum_{j \mathbf{k}} \frac{w_{\mathbf{k}} f_{j\mathbf{k}}}{\Omega} 
        !>                                      ( \nabla \phi_{\beta,\mathbf{G}+\mathbf{k}}^{j} (\mathbf{r}) ) \cdot 
        !>                                      ( \nabla \phi_{\beta,\mathbf{G^{'}}+\mathbf{k}}^{j} (\mathbf{r}) )^{*} \]
        !> \[ \tau_{\alpha \beta}(\mathbf{r}) = \frac{1}{2} \sum_{\mathbf{G} \mathbf{G}^{'}} \sum_{j \mathbf{k}} \frac{w_{\mathbf{k}} f_{j\mathbf{k}}}{\Omega} 
        !>                                      ( \nabla \phi_{\alpha,\mathbf{G}+\mathbf{k}}^{j} (\mathbf{r}) ) \cdot 
        !>                                      ( \nabla \phi_{\beta,\mathbf{G^{'}}+\mathbf{k}}^{j} (\mathbf{r}) )^{*} \]
        !>
        !> with the second-variational occupation number \(f_{j\mathbf{k}}\), the k-point weights \(w_{\mathbf{k}}\) and the unit cell volume \(\Omega\).
        !> In the spin-polarised case one needs to further differentiate between non-collinearism and collinearism.
        !> In the non-collinear case \(\tau(\mathbf{r})\) and the components of the spin kinetic energy density vector are calculated as:
        !> \[ \tau(\mathbf{r}) = \tau_{\alpha \alpha}(\mathbf{r}) + \tau_{\beta \beta}(\mathbf{r}) \]
        !> \[ u_{x}(\mathbf{r}) = \tau_{\alpha \beta}(\mathbf{r}) + \tau_{\beta \alpha}(\mathbf{r}) \]
        !> \[ u_{y}(\mathbf{r}) = i(\tau_{\alpha \beta}(\mathbf{r}) - \tau_{\beta \alpha}(\mathbf{r})) \]
        !> \[ u_{z}(\mathbf{r}) = \tau_{\alpha \alpha}(\mathbf{r}) - \tau_{\beta \beta}(\mathbf{r})  \]
        !> In the collinear case, \(u_{x}(\mathbf{r})\) and \(u_{y}(\mathbf{r})\) are zero. 
        subroutine gen_ked_ir_polarised(ik, evecfv, evecsv, ked_ir_k, uir_k)
            use mod_spin, only: ncmag
            !> i-th k-point
            integer, intent(in) :: ik
            !> first variational eigenvectors
            complex(dp), intent(in) :: evecfv (:, :) !currently spinorb true ignored
            !> second variational eigenvectors
            complex(dp), intent(in) :: evecsv (:, :)
            !> k-dependent kinetic energy density in the interstitial region
            real(dp), intent(inout) :: ked_ir_k(:)
            !> k-dependent spin kinetic energy density vector in the interstitial region
            real(dp), intent(inout) :: uir_k(:, :)

            !> Weights
            real(dp) :: wo, w
            integer :: ispn, jst, ir, grad_dir

            real(dp) :: up2(ked_Gset%ngrtot), dn2(ked_Gset%ngrtot)
            complex(dp) :: up(ked_Gset%ngrtot), dn(ked_Gset%ngrtot), updn(ked_Gset%ngrtot)

            complex(dp), allocatable :: zfft_grad(:, :, :), zfft(:, :)

            allocate( zfft_grad(ked_Gset%ngrtot, 3, nspinor) )
            allocate( zfft(ked_Gset%ngrtot, nspinor) )
            
            ked_ir_k = 0.0_dp
            uir_k = 0.0_dp
            do jst = 1, nstsv
                wo = 0.5_dp * ked_kset%wkpt(ik) * occsv(jst, ik)
                if (abs(wo) .gt. input%groundstate%epsocc) then
                    w = wo / omega

                    call gen_grad_wavefunction_ir_polarised(ik, evecfv, evecsv(:, jst), zfft_grad)

                    up2(:) = 0.0_dp
                    dn2(:) = 0.0_dp
                    updn(:) = 0.0_dp

                    do grad_dir = 1, 3
                        do ispn = 1, 2
                            call zfftifc(3, ked_Gset%ngrid, 1, zfft_grad(:, grad_dir, ispn))
                        end do

                        up2(:) = up2(:) + ( dble(zfft_grad(:, grad_dir, 1))**2 + aimag(zfft_grad(:, grad_dir, 1))**2 )
                        dn2(:) = dn2(:) + ( dble(zfft_grad(:, grad_dir, 2))**2 + aimag(zfft_grad(:, grad_dir, 2))**2 )
                        updn(:) = updn(:) + zfft_grad(:, grad_dir, 1) * conjg(zfft_grad(:, grad_dir, 2))
                    end do

                    ked_ir_k(:) = ked_ir_k(:) + w * ( up2(:) + dn2(:) )
                    if (ncmag) then 
                        ! non-collinear
                        uir_k(:, 1) = uir_k(:, 1) + 2.0_dp * w * dble( updn(:))
                        uir_k(:, 2) = uir_k(:, 2) - 2.0_dp * w * aimag( updn(:) )
                        uir_k(:, 3) = uir_k(:, 3) +  w * ( up2(:) - dn2(:) )
                    else
                        ! collinear
                        uir_k(:, 1) = uir_k(:, 1) + w * ( up2(:) - dn2(:) )
                    end if 
                end if
            end do 

            deallocate(zfft_grad, zfft)
        end subroutine 

        !> Given the first-variational eigenvector `evecfv_j` for the j-th state and the i-th \(\mathbf{k}\) point `ik`, 
        !> returns the gradient of the interstitial wavefunction. The gradient of the wavefunction is calculated and returned in reciprocal space. 
        !>
        !> In the spin-unpolarised case the following gradient is calculated: 
        !> \[ \nabla \phi_{\mathbf{G}j}^{\mathbf{k}} = i (\mathbf{G} + \mathbf{k}) C_{\mathbf{G} j}^{\mathbf{k}} \]
        !> where \(C_{\mathbf{G} j}^{\mathbf{k}}\) represents the first-variational eigenvector.
        subroutine gen_grad_wavefunction_ir_unpolarised(ik, evecfv_j, zfft_grad)
            !> i-th \(\mathbf{k}\)-point
            integer, intent(in) :: ik
            !> first variational eigenvector for the j-th state
            complex(dp), intent(in) :: evecfv_j(:) 
            !> gradient of wavefunction
            complex(dp), intent(inout) :: zfft_grad(:, :)

            integer :: grad_dir

            zfft_grad(:, :) = zzero
            do grad_dir = 1, 3
                zfft_grad(ked_Gset%igfft(ked_Gkset%igkig(1:ked_Gkset%ngk(1, ik), 1, ik)), grad_dir) = zi * ked_Gkset%vgkc(grad_dir, 1:ked_Gkset%ngk(1, ik), 1, ik) * evecfv_j(1:ked_Gkset%ngk(1, ik)) 
            end do 
        end subroutine 

        !> Given the first-variational eigenvector `evecfv` and the second-variational eigenvector `evecsv_j` for
        !> the j-th second variational state and the i-th k point `ik`, returns the gradient of the interstitial wavefunction.
        !> The gradient of the wavefunction is calculated and returned in reciprocal space. 
        !>
        !> In the spin-polarised case the following gradients are calculated:
        !> \[ \nabla \phi_{\alpha,\mathbf{G}j}^{\mathbf{k}} = \sum_{n=1}^{\text{nstfv}} i (\mathbf{G} + \mathbf{k}) C_{\mathbf{G}n}^{\mathbf{k}} B_{nj}^{\mathbf{k}} \]
        !> \[ \nabla \phi_{\beta,\mathbf{G}j}^{\mathbf{k}} = \sum_{n=\text{nstfv}+1}^{\text{nstsv}} i (\mathbf{G} + \mathbf{k}) C_{\mathbf{G}n}^{\mathbf{k}} B_{nj}^{\mathbf{k}} \]
        !> where \(C_{\mathbf{G}n}^{\mathbf{k}}\) is the first-variational eigenvector and \(B_{nl}^{\mathbf{k}}\) is the second-variational eigenvector.
        !> `nstsv` denotes the number of second-variational states and `nstfv` represents the number of first-variational states.
        subroutine gen_grad_wavefunction_ir_polarised(ik, evecfv, evecsv_j, zfft_grad)
            !> i-th k-point
            integer, intent(in) :: ik
            !> first variational eigenvectors
            complex(dp), intent(in) :: evecfv(:, :) !currently spinorb=true ignored
            !> second variational eigenvectors
            complex(dp), intent(in) :: evecsv_j(:)
            !> gradient of wave function
            complex(dp), intent(inout) :: zfft_grad(:, :, :)

            integer :: i, ispn, ist, grad_dir
            complex(dp) :: evecsv_ij
            complex(dp), allocatable :: evec_tmp(:)

            allocate(evec_tmp(ked_Gkset%ngk(1, ik)))

            zfft_grad(:, :, :) = zzero

            do ispn = 1, nspinor
                i = (ispn-1)*nstfv + 1
                call matrix_multiply(evecfv(1:ked_Gkset%ngk(1, ik), 1:nstfv), evecsv_j(i:i+nstfv-1), evec_tmp)
                call gen_grad_wavefunction_ir_unpolarised( ik, evec_tmp, zfft_grad(:, :, ispn) )
            end do
            deallocate ( evec_tmp )
        end subroutine 


end module
