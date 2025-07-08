!> Routines to calculate the mgga or gga potential to construct the basis.
module mgga_potxc
    use precision, only: dp
    use kinetic_energy_density_vars
    use mod_atoms, only: nspecies, natoms, idxas, natmtot
    use mod_SHT, only: rfshtvr, rbshtvr
    use modxcifc, only: xcifc
    use mod_spin, only: ndmag
    use modinput
    use mod_Gvector, only: ngrtot
    use mod_muffin_tin, only: nrmtmax, nrmt, lmmaxvr

    implicit none 

    private 

    !> interstitial non-multiplicative meta-GGA potential
    real(dp), public, allocatable :: vxcir_mgga_nonmult(:)
    !> muffin-tin non-multiplicative meta-GGA potential 
    real(dp), public, allocatable :: vxcmt_mgga_nonmult(:, :, :)

    public :: potxc_ir_spinunpolarised, potxc_mt_spinunpolarised
    public :: init_potxc_non_mult_mgga


    contains 

        !> Allocate all arrays needed to calculate the non-multiplicative 
        !> meta-GGA exchange-correlation potential
        subroutine init_potxc_non_mult_mgga()
            if ( allocated( vxcir_mgga_nonmult) ) deallocate(vxcir_mgga_nonmult )
            allocate(vxcir_mgga_nonmult(ngrtot))
            if ( allocated(vxcmt_mgga_nonmult) ) deallocate(vxcmt_mgga_nonmult)
            allocate(vxcmt_mgga_nonmult(lmmaxvr, nrmtmax, natmtot))
        end subroutine 

        !> In the case of a GGA potential, only calculates the exchange-correlation potential `vxcir` 
        !> and exchange and correlation energy densities `exir` and `ecir`. 
        !> In the case of a meta-GGA potential, calculates the multiplicative exchange-correlation potential `vxcir` 
        !> and the non-multiplicative potential `vxcir_non_mult' and the exchange and correlation energy densities `exir` and `ecir`. 
        !> `vxcir` is defined as follows:
        !> \[ 
        !>  v_{\text{xc}}^{\text{mult}}(\mathbf{r}) = \frac{\partial \varepsilon^{\text{mGGA}}_{\text{xc}}}{\partial n(\mathbf{r})} 
        !>                  - 2 \left[ \nabla \left(\frac{\partial \varepsilon^{\text{mGGA}}_{\text{xc}}}{\partial \gamma(\mathbf{r})}\right) \cdot \nabla n(\mathbf{r}) 
        !>                  + \frac{\partial \varepsilon^{\text{mGGA}}_{\text{xc}}}{\partial \gamma(\mathbf{r})} \nabla^2 n(\mathbf{r}) \right]. 
        !> \] 
        !> and `vxcir_non_mult` as follows:
        !> \[ 
        !>  v^{\text{non-mult}}_{\text{xc}, \tau}(\mathbf{r}) =  \frac{\partial \varepsilon^{\text{mGGA}}_{\text{xc}}}{\partial \tau(\mathbf{r})} 
        !> \]
        subroutine potxc_ir_spinunpolarised(xcgrad, xctype, rhoir, vxcir, exir, ecir, tauir, vxcir_non_mult)
            !> degree of exchange-correlation potential (can be 2 - for GGA or 3 - for meta-GGA)
            integer, intent(in) :: xcgrad
            !> exchange-correlation type
            integer, intent(in) :: xctype(:)
            !> density 
            real(dp), intent(in) :: rhoir(:)
            !> exchange-correlation potential 
            real(dp), intent(inout) :: vxcir(:)
            !> exchange energy density 
            real(dp), intent(inout) :: exir(:)
            !> correlation energy density 
            real(dp), intent(inout) :: ecir(:)
            !> kinetic energy density 
            real(dp), optional, intent(in) :: tauir(:)
            !> non-multiplicative exchange-correlation potential 
            real(dp), optional, intent(inout) :: vxcir_non_mult(:)

            integer :: n, nr, ir

            ! density
            real(dp), allocatable :: rho(:)
            ! gradient of density
            real(dp), allocatable :: gvrho(:)
            ! grad^2 of density 
            real(dp), allocatable :: g2rho(:)
            ! (grad of density)^2
            real(dp), allocatable :: grho2(:)
            ! exchange and correlation energy densities
            real(dp), allocatable :: ex(:), ec(:)
            ! exchange-correlation potentials
            real(dp), allocatable :: vxc(:), vx(:), vc(:)
            ! de_x/d(|grad rho|^2), de_c/d(|grad rho|^2)
            real(dp), allocatable :: dxdg2(:), dcdg2(:)
            ! non-multiplicative xc potential 
            real(dp), allocatable :: vxc_non_mult(:)
            ! kinetic energy densities
            real(dp), allocatable :: tau(:)
            ! de_xc/d(laplacian of rho)
            real(dp), allocatable :: dxdl(:), dcdl(:)
            ! de_xc/d(tau)
            real(dp), allocatable :: dxdtau(:), dcdtau(:)

            n = ngrtot
            allocate(rho(n), ex(n), ec(n), vxc(n), vx(n), vc(n))
            if (xcgrad == 2 .or. xcgrad == 3) then
                ! allocate everything for gga (libxc)
                allocate(g2rho(n),gvrho(3*n),grho2(n))
                allocate(dxdg2(n),dcdg2(n))
            end if 
            if (xcgrad == 3) then 
                ! allocate everything for mgga (libxc)
                allocate(tau(n), dxdl(n), dcdl(n), dxdtau(n), dcdtau(n))
            end if

            ! GGA 
            if (xcgrad == 2) then
                call ggair_2a(g2rho,gvrho,grho2)
                call xcifc(xctype,n=ngrtot,rho=rhoir,grho2=grho2,ex=exir,ec=ecir,vx=vx, vc=vc,dxdg2=dxdg2,dcdg2=dcdg2)
                 call ggair_2b(g2rho,gvrho,vx,vc,dxdg2,dcdg2)
                vxcir(1:ngrtot) = vx(1:ngrtot) + vc(1:ngrtot)
            
            ! meta-GGA
            else if (xcgrad == 3) then 
               ! multiplicative part of meta-GGA
                call ggair_2a(g2rho,gvrho,grho2)
                
                ! TASK functional uses LDA for correlation
                if (xctype(2) == 707) then 
                    call xcifc(xctype=xctype, n=ngrtot, rho=rhoir, grho2=grho2, g2rho=g2rho, tau=tauir, ex=exir, ec=ecir, vx=vx, & 
                              vc=vc, dxdg2=dxdg2, dxdl=dxdl, dxdtau=dxdtau)
                    dcdg2(:) = 0.0_dp
                    call ggair_2b(g2rho,gvrho,vx,vc,dxdg2,dcdg2)
                    
                    ! non-multiplicative part of meta-GGA
                    vxcir_non_mult(1:ngrtot) = dxdtau(1:ngrtot)
                else 
                    call xcifc(xctype=xctype, n=ngrtot, rho=rhoir, grho2=grho2, g2rho=g2rho, tau=tauir, ex=exir, ec=ecir, vx=vx, & 
                                vc=vc, dxdg2=dxdg2, dcdg2=dcdg2, dxdl=dxdl, dcdl=dcdl, dxdtau=dxdtau, dcdtau=dcdtau)
                    call ggair_2b(g2rho,gvrho,vx,vc,dxdg2,dcdg2)
                    
                    ! non-multiplicative part of meta-GGA
                    vxcir_non_mult(1:ngrtot) = (dcdtau(1:ngrtot) + dxdtau(1:ngrtot))
                end if 
                
                ! multiplicative potential 
                vxcir(1:ngrtot) = vx(1:ngrtot) + vc(1:ngrtot)

                deallocate(tau, dxdl, dcdl, dxdtau, dcdtau)
            end if

            deallocate(dxdg2, dcdg2, vx, vc)
            deallocate(g2rho, gvrho, grho2)
        end subroutine 

        !> In the case of a GGA potential, only calculates the exchange-correlation potential `vxcmt` 
        !> and exchange and correlation energy densities `exmt` and `ecmt`. 
        !> In the case of a meta-GGA potential, calculates the multiplicative exchange-correlation potential `vxcmt` 
        !> and the non-multiplicative potential `vxcmt_non_mult' and the exchange and correlation energy densities `exmt` and `ecmt`. 
        !> The multiplicative and non-multiplicative potentials are defined as in [[potxc_ir_spinunpolarised(subroutine)]].
        subroutine potxc_mt_spinunpolarised(xcgrad, xctype, rhomt, vxcmt, exmt, ecmt, taumt, vxcmt_non_mult)
            !> degree of exchange-correlation potential (can be 2 - for GGA or 3 - for meta-GGA)
            integer, intent(in) :: xcgrad
            !> exchange-correlation type
            integer, intent(in) :: xctype(:)
            !> density 
            real(dp), intent(in) :: rhomt(:, :, :)
            !> exchange-correlation potential 
            real(dp), intent(inout) :: vxcmt(:, :, :)
            !> exchange energy density 
            real(dp), intent(inout) :: exmt(:, :, :)
            !> correlation energy density 
            real(dp), intent(inout) :: ecmt(:, :, :)
            !> kinetic energy density 
            real(dp), optional, intent(in) :: taumt(:, :, :)
            !> non-multiplicative exchange-correlation potential 
            real(dp), optional, intent(inout) :: vxcmt_non_mult(:, :, :)

            integer :: is, ia, ias, lm, n, nr

            ! density
            real(dp), allocatable :: rho(:)
            ! gradient of density
            real(dp), allocatable :: gvrho(:)
            ! grad^2 of density 
            real(dp), allocatable :: g2rho(:)
            ! (grad of density)^2
            real(dp), allocatable :: grho2(:)
            ! exchange and correlation energy densities
            real(dp), allocatable :: ex(:), ec(:)
            ! exchange-correlation potentials
            real(dp), allocatable :: vxc(:), vx(:), vc(:)
            ! de_x/d(|grad rho|^2), de_c/d(|grad rho|^2)
            real(dp), allocatable :: dxdg2(:), dcdg2(:)
            ! non-multiplicative xc potential 
            real(dp), allocatable :: vxc_non_mult(:)
            ! kinetic energy densities
            real(dp), allocatable :: tau(:)
            ! de_xc/d(laplacian of rho)
            real(dp), allocatable :: dxdl(:), dcdl(:)
            ! de_xc/d(tau)
            real(dp), allocatable :: dxdtau(:), dcdtau(:)

            n = lmmaxvr*nrmtmax
            
            allocate(rho(n), ex(n), ec(n), vxc(n), vx(n), vc(n))
            if (xcgrad == 2 .or. xcgrad == 3) then
                ! allocate everything for gga 
                allocate(g2rho(n),gvrho(3*n),grho2(n))
                allocate(dxdg2(n),dcdg2(n))
            end if 
            if (xcgrad == 3) then 
                ! allocate everything for mgga
                allocate(tau(n), dxdl(n), dcdl(n), dxdtau(n), dcdtau(n), vxc_non_mult(n))
            end if

            do is = 1, nspecies
                nr = nrmt(is)
                n = lmmaxvr*nr
                do ia = 1, natoms(is)
                    ias = idxas(ia,is)
                    ! compute the density in spherical coordinates
                    call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, rhomt(:,:,ias), &
                        lmmaxvr, 0.d0, rho, lmmaxvr)

                    ! GGA potential
                    if (xcgrad == 2) then
                        call ggamt_2a(is, ia, g2rho, gvrho, grho2)
                        call xcifc(xctype,n=n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
                                    dxdg2=dxdg2,dcdg2=dcdg2)
                        call ggamt_2b(is, g2rho, gvrho, vx, vc, dxdg2, dcdg2)
                        
                        vxc(1:n) = vx(1:n)+ vc(1:n) 
                        call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, vxc, lmmaxvr, &
                                0.d0, vxcmt(:,:,ias), lmmaxvr)
                    
                    ! meta-GGA potential
                    else if (xcgrad == 3) then 
                        call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, taumt(:,:,ias), &
                        lmmaxvr, 0.d0, tau, lmmaxvr)
                        
                        ! multiplicative part of meta-GGA
                        call ggamt_2a(is, ia, g2rho, gvrho, grho2)

                        if (xctype(2) == 707) then 
                            call xcifc(xctype=xctype, n=n, rho=rho, grho2=grho2, g2rho=g2rho, tau=tau, ex=ex, ec=ec, vx=vx, & 
                                    vc=vc, dxdg2=dxdg2, dxdl=dxdl, dxdtau=dxdtau)
                            dcdg2(:) = 0.0_dp
                            call ggamt_2b(is, g2rho, gvrho, vx, vc, dxdg2, dcdg2) 
                            
                            ! non-multiplicative part of meta-GGA
                            vxc_non_mult(1:n) = dxdtau(1:n)
                        else
                            call xcifc(xctype=xctype, n=n, rho=rho, grho2=grho2, g2rho=g2rho, tau=tau, ex=ex, ec=ec, vx=vx, & 
                                        vc=vc, dxdg2=dxdg2, dcdg2=dcdg2, dxdl=dxdl, dcdl=dcdl, dxdtau=dxdtau, dcdtau=dcdtau)
                            call ggamt_2b(is, g2rho, gvrho, vx, vc, dxdg2, dcdg2) 
                        
                            ! non-multiplicative part of meta-GGA
                            vxc_non_mult(1:n) = (dcdtau(1:n) + dxdtau(1:n))
                        end if 
                        ! multiplicative part of meta_GGA
                        vxc(1:n) = vx(1:n)+ vc(1:n) 

                        call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, vxc, lmmaxvr, &
                                0.d0, vxcmt(:,:,ias), lmmaxvr)
                        call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, vxc_non_mult, lmmaxvr, &
                               0.d0, vxcmt_non_mult(:,:,ias), lmmaxvr)
                    end if

                    ! convert exchange and correlation energy densities to spherical harmonics
                    call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, ex, lmmaxvr, &
                            0.d0, exmt(:,:,ias), lmmaxvr)
                    call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, ec, lmmaxvr, &
                            0.d0, ecmt(:,:,ias), lmmaxvr)
                end do 
            end do 
            
            if (xcgrad == 2 .or. xcgrad == 3) then
                ! deallocate everything for gga (libxc)
                deallocate(g2rho,gvrho,grho2)
                deallocate(dxdg2,dcdg2)
            end if 
            if (xcgrad == 3) then 
                ! deallocate everything for mgga
                deallocate(tau, dxdl, dcdl, dxdtau, dcdtau, vxc_non_mult)
            end if

        end subroutine 

end module 