module mod_polarizability_q_ii
    use precision, only: wp
    use constants, only: zzero, zone
    use mod_Gvector, only: igfft
    use modgw, only: Gset, kqset, Gqset
    use modinput, only: input
    use modmain, only : ngkmax   !! max value of G-vectors in all k-points
    use mod_overlap_ii, only: overlap_ii
    use mod_polarizability_R_ii, only: polarizability_R_ii

    implicit none
    private
    public :: polarizability_q_ii, polarizability_tau_to_omega_ii

    contains

    !!------------------------------------------------------------------------------------------------
    !! FFT of polarizability from R to q-space
    !> \begin{equation}
    !>   P^{\mathbf q}_{\mathbf r, \mathbf r'}(\tau) = \sum_{\mathbf R} e^{-i \mathbf{q} \mathbf{R}} 
    !>   P^{\mathbf R}_{\mathbf r, \mathbf r'}(\tau)
    !> \end{equation}
    !>
    !> Fourier transform of polarizability from r' to G'
    !> \begin{equation}
    !>   P^{\mathbf q}_{\mathbf r, \mathbf G'}(\tau) = \frac{1}{N_{\mathbf r}} \sum_{\mathbf r'} 
    !>   e^{-i (\mathbf{q} + \mathbf{G}') \mathbf{r}'} P^{\mathbf q}_{\mathbf r, \mathbf r'}(\tau)
    !> \end{equation}
    !>
    !> Fourier transform of polarizability from r to G
    !> \begin{equation}
    !>   P^{\mathbf q}_{\mathbf G, \mathbf G'}(\tau) = \frac{1}{N_{\mathbf r}} \sum_{\mathbf r} 
    !>   e^{i (\mathbf{q} + \mathbf{G}) \mathbf{r}} P^{\mathbf q}_{\mathbf r, \mathbf G'}(\tau)
    !> \end{equation}
    !>
    !> Polarizability in mixed basis representation by multiplying with I-I overlap matrix
    !> from both side (left and right)
    !> \begin{equation}
    !>   P_{\mathbf{G}, I}^{{\mathbf{q}}} \left(\tau\right) = \sum_{\mathbf{G}^{\prime}}  
    !>   P_{\mathbf{G}, \mathbf{G}^{\prime}}^{{\mathbf{q}}}\left(\tau\right) \left\langle  
    !>   e^{i(\mathbf{q}+ \mathbf{G}^{\prime})\mathbf{r}^{\prime} } \mid 
    !>   \chi^{\mathbf{q}}_{I^{\prime}}\right\rangle_{IR}
    !> \end{equation}
    !> \begin{equation}
    !>   P_{I, I}^{{\mathbf{q}}} \left(\tau\right) = \sum_{\mathbf{G}}  \left\langle  
    !>   e^{i(\mathbf{q}+ \mathbf{G})\mathbf{r}} \mid \chi^{\mathbf{q}}_{I}\right\rangle_{IR}^{*} 
    !>   P_{\mathbf{G}, I^{\prime}}^{{\mathbf{q}}} \left(\tau\right)
    !> \end{equation}
    !>
    !> where, 
    !> \begin{equation}
    !>   \left\langle  e^{i(\mathbf{q}+ \mathbf{G}^{\prime})\mathbf{r}^{\prime} } \mid 
    !>   \chi^{\mathbf{q}}_{I^{\prime}}\right\rangle_{IR} = S_{\mathbf{G}',I'}^{\mathbf{q}}
    !> \end{equation}
    !!------------------------------------------------------------------------------------------------

    subroutine polarizability_q_ii(tau, pola_q_ii)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Polarizability in imaginary time and q representation in I-I
        complex(wp), allocatable, intent(out) :: pola_q_ii(:,:,:)

        !> Local variables
        integer :: iq, ir_grid1, ir_grid2
        !> Total no of G-vectors per k-point
        integer :: ngq
        !> Index for G-vector
        integer :: ig, igq1, igq2
        !> Polarizability in imaginary time and R/q representation in I-I
        complex(wp), allocatable :: pola_q_ii_tmp(:,:,:)
        !> Temporary matrix for polarizability 
        complex(wp), allocatable :: tmp(:,:)
        !> Overlap matrix (between two plane waves) in I-I region
        complex(wp), allocatable :: s_ii(:,:)
        !> Overlap matrix in fft grid
        complex(wp), allocatable :: s_ii_fft(:,:)

        !!! USE and delete 
        real(wp) :: t_i, t_f

        print*, 'from polarizability_q_ii input%gw%ngridq', input%gw%ngridq
        print*, 'from polarizability_q_ii Gset%ngrtot', Gset%ngrtot
        print*, 'ngkmax', ngkmax

        allocate(pola_q_ii(ngkmax, ngkmax, kqset%nkpt))
        pola_q_ii = zzero

        call polarizability_R_ii(tau, pola_q_ii_tmp)

        call CPU_TIME(t_i)
        ! FFT of polarizability function from R to q-space
        do ir_grid1 = 1, Gset%ngrtot
          do ir_grid2 = 1, Gset%ngrtot
            call zfftifc(3, input%gw%ngridq, -1, pola_q_ii_tmp(ir_grid1, ir_grid2, :))
          enddo
        enddo
        call CPU_TIME(t_f)
        print*, 'time, sum pola_q_ii_tmp', t_f - t_i, sum(pola_q_ii_tmp)
        print*, 'test pola_q_ii_tmp', minval(real(pola_q_ii_tmp)), maxval(imag(pola_q_ii_tmp))


        ! FFT from r' to G'-space of polarizability
        call CPU_TIME(t_i)
        do iq = 1, kqset%nkpt
          do ir_grid1 = 1, Gset%ngrtot
            call zfftifc(3, Gset%ngrid, -1, pola_q_ii_tmp(ir_grid1, :, iq))
            !!! TODO: do I need Gqset%ngrid or Gset%ngrid and -1 or 1
          enddo
        enddo 
        pola_q_ii_tmp = pola_q_ii_tmp / real(Gset%ngrtot, kind=wp)
        call CPU_TIME(t_f)
        print*, "pola_q_ii_tmp(r,G') and time", sum(pola_q_ii_tmp), t_f - t_i
        print*, 'test pola_q_ii_tmp', minval(real(pola_q_ii_tmp)), maxval(imag(pola_q_ii_tmp))

        ! FFT from r to G-space of polarizability
        call CPU_TIME(t_i)
        do iq = 1, kqset%nkpt
          do ir_grid2 = 1, Gset%ngrtot
            call zfftifc(3, Gset%ngrid, 1, pola_q_ii_tmp(:, ir_grid2, iq))
            ! call zfftifc(3, Gset%ngrid, -1, pola_q_ii_tmp(:, ir_grid2, iq))
            ! !! TODO: do I need Gqset%ngrid or Gset%ngrid and 1 or -1
          enddo 
        enddo 
        pola_q_ii_tmp = pola_q_ii_tmp / real(Gset%ngrtot, kind=wp)
        call CPU_TIME(t_f)
        print*, "pola_q_ii_tmp(G,G') and time", sum(pola_q_ii_tmp), t_f - t_i
        print*, 'test pola_q_ii_tmp', minval(real(pola_q_ii_tmp)), maxval(imag(pola_q_ii_tmp))


        !! Polarizability in mixed basis representation

        allocate(s_ii_fft(Gset%ngrtot, ngkmax))
        allocate(tmp(Gset%ngrtot, ngkmax))
        ! allocate(tmp(Gset%ngrtot, ngq))

        call CPU_TIME(t_i)
        do iq = 1, kqset%nkpt

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
          !allocate(pola_q_ii_tmp(Gset%ngrtot, Gset%ngrtot, nbigR))

          ! call zgemm('n', 'n', Gset%ngrtot, ngq, Gset%ngrtot,  zone, pola_q_ii_tmp(1,1, iq), Gset%ngrtot, &
          call zgemm('t', 'n', Gset%ngrtot, ngq, Gset%ngrtot,  zone, pola_q_ii_tmp(1,1, iq), Gset%ngrtot, &
                     s_ii_fft(1,1), Gset%ngrtot,  zzero, tmp(1,1), Gset%ngrtot)

          ! call zgemm('t', 'n', ngq, ngq, Gset%ngrtot,  zone, s_ii_fft(1,1), Gset%ngrtot, & 
          call zgemm('c', 'n', ngq, ngq, Gset%ngrtot,  zone, s_ii_fft(1,1), Gset%ngrtot, & 
                     tmp(1,1), Gset%ngrtot,  zzero, pola_q_ii(1,1, iq), ngkmax)
        enddo   !! iq loop ends
        call CPU_TIME(t_f)
        print*, "pola_q_ii(I,I') and time", sum(pola_q_ii), t_f - t_i
        print*, 'test pola_q_ii', minval(real(pola_q_ii)), maxval(imag(pola_q_ii))

    end subroutine polarizability_q_ii

    subroutine polarizability_tau_to_omega_ii(jt, weight_cos, pola_q_ii, pola_q_ii_omega)
        !> Index for $\tau$
        integer :: jt
        !> Weight and cos(w*t) for tau to omega transformation
        real(wp), intent(in) :: weight_cos(:,:)
        !> Polarizability in imaginary time and R representation in I-I
        ! complex(wp), allocatable, intent(in) :: pola_q_ii(:,:,:,:,:)
        complex(wp), intent(in) :: pola_q_ii(:,:,:)
        !> Polarizability in R representation and in omega domain (in I-I)
        complex(wp), intent(out) :: pola_q_ii_omega(:,:,:,:)

        !! Local variables
        integer :: iq, igq1, igq2, it, n_tau


        print*, 'shape of weight_cos', shape(weight_cos)
        print*, 'shape of pola_q_ii', shape(pola_q_ii)

        ! allocate(pola_q_ii_omega, mold=pola_q_ii)

        print*, 'shape of pola_q_ii_omega', shape(pola_q_ii_omega)
        print*, 'size of pola_q_ii_omega', size(pola_q_ii_omega, 1)

        n_tau = size(weight_cos, 1)
        print*, ' no. of tau points, nlm_tot', n_tau

        do it = 1, n_tau
          do iq = 1, kqset%nkpt
            do igq1 = 1, ngkmax  !! need to change with actual dimension
              do igq2 = 1, ngkmax  !! need to change with actual dimension
                pola_q_ii_omega(igq2, igq1, iq, it) = & 
                pola_q_ii_omega(igq2, igq1, iq, it) + &
                weight_cos(it, jt) * pola_q_ii(igq2, igq1, iq)
              enddo 
            enddo 
          enddo 
        enddo 

    end subroutine polarizability_tau_to_omega_ii

end module mod_polarizability_q_ii