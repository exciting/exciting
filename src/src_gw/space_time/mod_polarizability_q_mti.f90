module mod_polarizability_q_mti
    use precision, only: wp
    use constants, only: zzero, zone
    use mod_atoms, only: natmtot
    use mod_Gvector, only: igfft
    use modgw, only: kqset, Gset, Gqset
    use mod_muffin_tin, only: lmmaxapw
    use modinput, only: input   !! This is for "input%gw%ngridq"
    use mod_product_basis, only: locmatsiz
    use mod_bravais_lattice, only: bravais_lattice
    use mod_polarizability_R_mti, only: polarizability_R_mti
    use mod_overlap_ii, only: overlap_ii
    use modmain, only : ngkmax   !! max value of G-vectors in all k-points
    !! TEST AND DELETE
    use modmain, only: ngrtot
    
    implicit none
    private
    public :: polarizability_q_mti, polarizability_tau_to_omega_mti

    contains

    !!-----------------------------------------------------------------------------------------------
    !! Polarizabiliity in MT-I region in q-space from R-sapce (in 3 steps)
    !! First: FFT from R to q-space 
    !> \begin{equation}
    !>   P_{\alpha I, \mathbf{r}^{\prime} }^{{\mathbf{q}}}\left(\tau\right)
    !>   =\sum_{\mathbf{R}} e^{-i \mathbf{q} \mathbf{R}}
    !>   P_{\alpha I, \mathbf{r}^{\prime}}^{{\mathbf{R}}}\left(\tau\right)
    !> \end{equation}
    !> 
    !> Second: FFT from r to G-space
    !> \begin{equation}
    !>   P_{\alpha I, \mathbf{G}^{\prime} }^{{\mathbf{q}}}\left(\tau\right)
    !>   = \frac{1}{N_{\mathbf{r}}}\sum_{\mathbf{r}'} e^{i (\mathbf{q} + \mathbf{G}') \mathbf{r}'}
    !>   P_{\alpha I, \mathbf{r}^{\prime}}^{{\mathbf{q}}}\left(\tau\right)
    !> \end{equation}
    !> 
    !> Third: from G-space to mixed basis representation by multiplying with I-I overlap matrix
    !> \begin{equation}
    !>   P_{\alpha I, I^{\prime}}^{{\mathbf{q}}} \left(\tau\right)=\sum_{\mathbf{G}^{\prime}}  
    !>   P_{\alpha I, \mathbf{G}^{\prime}}^{{\mathbf{q}}}\left(\tau\right) \left\langle  
    !>   e^{i(\mathbf{q}+ \mathbf{G}^{\prime})\mathbf{r}^{\prime} } 
    !>   \mid \chi^{\mathbf{q}}_{I^{\prime}}\right\rangle_{IR}
    !> \end{equation}
    !> 
    !> where, 
    !> \begin{equation}
    !>   \left\langle  e^{i(\mathbf{q}+ \mathbf{G}^{\prime})\mathbf{r}^{\prime} } \mid 
    !>   \chi^{\mathbf{q}}_{I^{\prime}}\right\rangle_{IR} = S_{\mathbf{G}',I'}^{\mathbf{q}}
    !> \end{equation}
    !!-----------------------------------------------------------------------------------------------

    subroutine polarizability_q_mti(tau, pola_q_mti)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Polarizability in imaginary time and q representation in MT-I
        complex(wp), allocatable, intent(out) :: pola_q_mti(:,:,:,:)

        ! Local variables 
        integer :: ias, nlm
        !> Index for the Bravias vector R 
        integer :: ir_grid, iq
        !> Total number of Bravias vector
        integer :: nbigR
        !> Total number of n, l, m
        integer :: nlm_tot
        !> Total no of G-vectors per k-point
        integer :: ngq
        !> Index for G-vector
        integer :: ig, igq1, igq2
        !> Polarizability in imaginary time and R representation in MT-I
        complex(wp), allocatable :: pola_R_mti(:,:,:,:)
        !> Polarizability in imaginary time and q representation in MT-I with G' vector
        complex(wp), allocatable :: pola_q_mti_g(:,:,:,:)
        !> Temporary array for FFT 
        complex(wp), allocatable :: ztmp(:)
        !> Overlap matrix (between two plane waves) in I-I region
        complex(wp), allocatable :: s_ii(:,:)
        !> Overlap matrix in fft grid
        complex(wp), allocatable :: s_ii_fft(:,:)

        ! Test and delete
        real :: t_i, t_f

        !> External routines for matrix matrix multiplication
        external :: zgemm

        
        nbigR = kqset%nkpt
        nlm_tot = locmatsiz/natmtot

        print*, 'Gset%ngrtot, ngrtot, ngkmax', Gset%ngrtot, ngrtot, ngkmax
        print*, 'shape(pola_q_mti)', shape(pola_q_mti)
        
        allocate(pola_q_mti_g(Gset%ngrtot, nlm_tot, natmtot, nbigR))
        allocate(pola_q_mti(nlm_tot, natmtot, ngkmax, nbigR))
        allocate(s_ii_fft(Gset%ngrtot, ngkmax))
        
        ! pola_q_mti_g = cmplx(0.0_wp, 0.0_wp, kind=wp)
        ! pola_q_mti = cmplx(0.0_wp, 0.0_wp, kind=wp)
        pola_q_mti_g = zzero
        pola_q_mti = zzero

        !> Polarizability in R-space in MT-MT
        call polarizability_R_mti(tau, pola_R_mti)

        allocate(ztmp(nbigR))

        !! FFT from R to q space
        call CPU_TIME(t_i)
        do ias = 1, natmtot
          do nlm = 1, nlm_tot
            do ir_grid = 1, Gset%ngrtot

              ztmp(:) = pola_R_mti(ir_grid, nlm, ias, :)
              call zfftifc(3, input%gw%ngridq, -1, ztmp)
              ! ! call zfftifc(3, input%gw%ngridq, -1, pola_R_mti(ir_grid, nlm, ias, :))
              pola_q_mti_g(ir_grid, nlm, ias, :) = ztmp(:)

            enddo
          enddo
        enddo
        write(4011,*) sum(real(pola_q_mti_g)), sum(imag(pola_q_mti_g))

        !! FFT from r' to G' sapce  !!! TODO: Do I need Gset or Gqset in *%ngrid???
        do iq = 1, kqset%nkpt
          do ias = 1, natmtot
            do nlm = 1, nlm_tot

              call zfftifc(3, Gset%ngrid, 1, pola_q_mti_g(:, nlm, ias, iq))   !!! Do I need -1, or 1 (see eqn.)
              ! call zfftifc(3, Gset%ngrid, -1, pola_q_mti_g(:, nlm, ias, iq))   !!! Do I need -1, or 1 (see eqn.)
              !!! TODO: do we need to re-arrange the G'(r' to G') data ?
              !!! as we did for the reverse (G' to r') (see "mod_green_R_mti")

            enddo
          enddo
          pola_q_mti_g(:,:,:, iq) = pola_q_mti_g(:,:,:, iq) / real(Gset%ngrtot, kind=wp)

          write(4012,*) sum(real(pola_q_mti_g)), sum(imag(pola_q_mti_g))


          !! Polarizability in mixed basis representation

          call overlap_ii(iq, s_ii)

          ngq = Gqset%ngk(1, iq)
          s_ii_fft(:, :) = zzero
          do igq2 = 1, ngq!size(s_ii(1, :))
            ! s_ii_fft(:, igq2) = zzero
            do igq1 = 1, ngq!size(s_ii(:, 1))
              ig = igfft(Gqset%igkig(igq1, 1, iq))
              s_ii_fft(ig, igq2) = s_ii(igq1, igq2)
            enddo
          enddo

          call zgemm('t', 'n', nlm_tot*natmtot, ngq, Gset%ngrtot, &
                     zone, pola_q_mti_g(1,1,1, iq), Gset%ngrtot,   &
                     s_ii_fft(1,1), Gset%ngrtot, zzero, pola_q_mti(1, 1, 1, iq), nlm_tot*natmtot)


        enddo   !! iq loop ends
        call CPU_TIME(t_f)

        print*, 'sum pola_q_mti_g (real and imag), time', sum(real(pola_q_mti_g)), sum(imag(pola_q_mti_g)), t_f-t_i
        print*, 'sum pola_q_mti   (real and imag), time', sum(real(pola_q_mti)), sum(imag(pola_q_mti)), t_f-t_i
        write(4013,*) sum(real(pola_q_mti)), sum(imag(pola_q_mti))


    end subroutine polarizability_q_mti

    !!---------------------------------------------------------------------------
    !! Transformation of the polarizability from $\omega$ to $\tau$ domain
    !> \begin{equation}
    !>   P_{\alpha I, I^{\prime}}^{{\mathbf{q}}}\left(\omega_i\right)
    !>   =\sum_{j=1}^{n} \delta_{ij} cos(\omega_i \tau_j) 
    !>   P_{\alpha I, I^{\prime}}^{{\mathbf{q}}}\left(\tau_j\right)
    !> \end{equation}
    !!---------------------------------------------------------------------------

    subroutine polarizability_tau_to_omega_mti(jt, weight_cos, pola_q_mti, pola_q_mti_omega)
        !> Index for $\tau$
        integer :: jt
        !> Weight and cos(w*t) for tau to omega transformation in MT-I
        real(wp), intent(in) :: weight_cos(:,:)
        !> Polarizability in imaginary time and R representation in MT-I
        complex(wp), intent(in) :: pola_q_mti(:,:,:,:)
        !> Polarizability in R representation and in omega domain (in MT-I) 
        complex(wp), intent(out) :: pola_q_mti_omega(:,:,:,:,:)

        !! Local variables
        integer :: iq, igq, ias, nlm, it, n_tau
        !> Total number of n, l, m
        integer :: nlm_tot

        nlm_tot = locmatsiz/natmtot
        n_tau = size(weight_cos, 1)
        !allocate(pola_q_mti(nlm_tot, natmtot, ngkmax, nbigR))

        do it = 1, n_tau
          do iq = 1, kqset%nkpt
            do igq = 1, ngkmax
              do ias = 1, natmtot
                do nlm = 1, nlm_tot
                  pola_q_mti_omega(nlm, ias, igq, iq, it) = pola_q_mti_omega(nlm, ias, igq, iq, it) + &
                  weight_cos(it, jt) * pola_q_mti(nlm, ias, igq, iq)
                enddo
              enddo
            enddo
          enddo
        enddo 

    end subroutine polarizability_tau_to_omega_mti

end module mod_polarizability_q_mti
