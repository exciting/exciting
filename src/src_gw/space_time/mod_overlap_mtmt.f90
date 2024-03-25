module mod_overlap_mtmt
    use precision, only: wp
    ! use constants, only: pi
    use modinput, only: input
    ! use mod_misc_gw, only: atposl
    ! use modgw, only: kqset
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_APW_LO, only: apword, apwordmax, nlorb, lorbl, lolmmax
    use mod_muffin_tin, only: idxlm, lmmaxapw
    use mod_gaunt_coefficients, only: getgauntcoef
    use mod_product_basis, only: bradketa, bradketlo

    !! test and then DELETE later
    use mod_product_basis,     only: matsiz, mbsiz, locmatsiz
    
    implicit none
    private
    public :: overlap_mtmt

    contains

    !!------------------------------------------------------------------------------------------------------------
    !> Overlap matrix in MT-MT region 
    !>
    !> \begin{equation}
		!>   S^{\alpha}_{I,pp'} = 
    !>    \left\langle\chi_I^\alpha \mid f_p^\alpha f_{p^{\prime}}^\alpha\right\rangle 
		!>    = \left\langle \gamma^{\mathbf{q}}_{\alpha NLM}(\mathbf{r}) \mid 
		!>    f_{\xi \ell m}^\alpha f_{{\xi' \ell' m'}}^\alpha\right\rangle 
    !> \end{equation} 
    !> 
    !> where, $$ \gamma^{\mathbf{q}}_{\alpha NLM}(\mathbf{r}) = e^{i \mathbf{q} \mathbf {r}_{\alpha}} 
    !>                                       {v}_{\alpha NL}(r_\alpha) Y_{LM}(\hat{\mathbf r}^{\alpha}) $$ 
    !> and $$ f_{\xi \ell m}^\alpha =  u^{\alpha}_{\xi \ell}(r^{\alpha}) Y_{\ell m}(\hat{\mathbf r}^\alpha) $$ 
    !> Therefore, 
    !> \begin{equation} 
		!>   S^{\alpha I}_{\xi \ell m ; \xi' \ell' m'} = 
		!>    \left\langle e^{i \mathbf{q} \mathbf {r}_{\alpha}} {v}_{\alpha NL}(r_\alpha) 
		!>    Y_{LM}(\hat{\mathbf r}^{\alpha}) \mid u^{\alpha}_{\xi \ell}(r^{\alpha}) 
		!>    Y_{\ell m}(\hat{\mathbf r}^\alpha) 
		!>    u^{\alpha}_{\xi' \ell'}(r^{\alpha}) Y_{\ell' m'}(\hat{\mathbf r}^\alpha) 
		!>    \right\rangle 
    !> \end{equation} 
    !> Finally 
    !> \begin{equation} 
		!>   S^{\alpha I}_{\xi \ell m ; \xi' \ell' m'} = 
		!>    e^{- i \mathbf{q} \mathbf {r}_{\alpha}} 
    !>    \int_{\Omega} d\hat{\mathbf r} 
    !>    Y_{LM}^{*}(\hat{\mathbf r}^{\alpha}) 
    !>    Y_{\ell m}(\hat{\mathbf r}^\alpha) 
    !>    Y_{\ell' m'}(\hat{\mathbf r}^\alpha)
    !>    \int_{0}^{R^{\alpha}_{MT}} r^2 dr
    !>    {v}_{\alpha NL}(r_\alpha) 
		!>    u^{\alpha}_{\xi \ell}(r^{\alpha}) 
		!>    u^{\alpha}_{\xi' \ell'}(r^{\alpha}) 
    !> \end{equation}
    !> Second integration of the last equation can be APW and LO. Depending on this G will have 4 combinations 
    !> (1) APW-APW  (2) APW-LO  (3) LO-APW  (4) LO-LO
    !!
    !! TODO(Manoar): Do we need phase factor $$ e^{i \mathbf{q} \mathbf {r}_{\alpha}} $$ ?????
    !!               If needed then implement it
    !!------------------------------------------------------------------------------------------------------------

    subroutine overlap_mtmt(nlm_indices, s_mtmt)
        !> Combined N, L, M index
        integer, intent(in) :: nlm_indices(:, :, :)
        !> Overlap matrix: MT-MT case
        real(wp), allocatable, intent(out) :: s_mtmt(:,:,:,:,:,:,:,:)

        !! Local variables
        !> Max of lm including APW and LO
        integer :: lm_max
        integer :: ia, is, ias
        integer :: mbn, mbl, mbm
        integer :: l1, m1, lm1, l2, m2, lm2
        integer :: l2min, l2max
        integer :: io1, io2
        integer :: ilo1, ilo2
        integer :: nlm, nlm_tot
        !> Gaunt coefficient 
        real(wp) :: gaunt_coefficient
        ! !> Phase 
        ! real(wp) :: phase
        ! !> Phase factor 
        ! complex(wp) :: phase_factor 
        !--------------------------------------------------------------

        nlm_tot = size(nlm_indices(:, 1, 1))
        lm_max = max(lmmaxapw, lolmmax)

        allocate(s_mtmt(apwordmax, lm_max, apwordmax, lm_max, nlm_tot, natmtot, 3, 3))

        s_mtmt = 0.0_wp
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            ! ! phase = atposl(1,ia,is)*kqset%vql(1, iq)+ atposl(2,ia,is)*kqset%vql(2, iq)+ atposl(3,ia,is)*kqset%vql(3, iq)
            ! ! ! phase = atposl(1,ia,is)*kqset%vql(1, 5)+ atposl(2,ia,is)*kqset%vql(2, 5)+ atposl(3,ia,is)*kqset%vql(3, 5)
            ! phase = dot_product(atposl(:, ia, is), kqset%vql(:, iq))
            ! phase = 2.0_wp * pi * phase
            ! phase_factor = cmplx(cos(phase), -sin(phase), kind=wp)
            ! print*, 'phase, phase_factor', phase, phase_factor
            do nlm = 1, nlm_tot
              mbn = nlm_indices(nlm, ias, 3)
              mbl = nlm_indices(nlm, ias, 4)
              mbm = nlm_indices(nlm, ias, 5)
              do l1 = 0, input%groundstate%lmaxapw
                do m1 = -l1, l1
                  lm1 = idxlm(l1,m1)
                  m2 = - mbm + m1
                  l2min = abs(m2)
                  l2max = min(mbl+l1, input%groundstate%lmaxapw)
                  do l2 = l2min, l2max
                    lm2 = idxlm(l2,m2)
                    gaunt_coefficient = getgauntcoef(l2, mbl, l1, m2, mbm)

                    !! APW - APW/LO ---------------------------------------------------------
                    do io1 = 1, apword(l1, is)
                      !! APW-APW
                      do io2 = 1, apword(l2, is)
                        s_mtmt(io2, lm2, io1, lm1, nlm, ias, 3, 3) = & 
                        gaunt_coefficient * bradketa(2, mbn, l1, io1, l2, io2, ias)
                      enddo   !! io2
                      !! APW-LO
                      io2 = 1
                      do ilo2 = 1, nlorb(is)
                        if(lorbl(ilo2, is)==l2) then
                          s_mtmt(io2, lm2, io1, lm1, nlm, ias, 3, 2) = & 
                          gaunt_coefficient * bradketa(3, mbn, l1, io1, l2, io2, ias)
                        endif 
                      enddo   !! ilo2
                    enddo   !! io1

                    !! LO - APW/LO ----------------------------------------------------------
                    io1 = 1
                    do ilo1 = 1, nlorb(is)
                      if(lorbl(ilo1, is)==l1) then 
                        !! LO-APW
                        do io2 = 1, apword(l2, is)
                          s_mtmt(io2, lm2, io1, lm1, nlm, ias, 2, 3) = & 
                          gaunt_coefficient * bradketlo(2, mbn, ilo1, l2, io2, ias) 
                        enddo   !! io2
                        !! LO-LO
                        io2 = 1
                        do ilo2 = 1, nlorb(is)
                          if(lorbl(ilo2, is)==l2) then
                            s_mtmt(io2, lm2, io1, lm1, nlm, ias, 2, 2) = & 
                            gaunt_coefficient * bradketlo(2, mbn, ilo1, l2, io2, ias) 
                          endif 
                        enddo   !! ilo2
                      endif 
                    enddo   !! ilo1

                  enddo   !! l2
                enddo   !! m1
              enddo   !! l1
            enddo   !! nlm
          enddo   !! ia
        enddo   !! is

        print*, 'min, max overlap new ', minval(s_mtmt(:,:,:,:,:,:, 3, 3)), maxval(s_mtmt(:,:,:,:,:,:, 3, 3))
        print*, 'min, max overlap new ', minval(s_mtmt(:,:,:,:,:,:, 3, 2)), maxval(s_mtmt(:,:,:,:,:,:, 3, 2))
        print*, 'min, max overlap new ', minval(s_mtmt(:,:,:,:,:,:, 2, 3)), maxval(s_mtmt(:,:,:,:,:,:, 2, 3))
        print*, 'min, max overlap new ', minval(s_mtmt(:,:,:,:,:,:, 2, 2)), maxval(s_mtmt(:,:,:,:,:,:, 2, 2))


    end subroutine overlap_mtmt

end module mod_overlap_mtmt
