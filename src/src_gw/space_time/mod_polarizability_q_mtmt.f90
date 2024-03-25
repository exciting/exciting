module mod_polarizability_q_mtmt
    use precision, only: wp
    use constants, only: zzero
    use modinput, only: input   !! This is for "input%gw%ngridq"
    use mod_atoms, only: natmtot
    use modgw, only: kqset, qsetd
    use mod_product_basis, only: locmatsiz
    use mod_bravais_lattice, only: bravais_lattice
    use mod_polarizability_R_mtmt, only: polarizability_R_mtmt
    
    implicit none
    private
    public :: polarizability_q_mtmt, polarizability_tau_to_omega_mtmt

    contains

    !!----------------------------------------------------------------------------------
    !! Polarizabiliity in MT-MT region in q-space from R-sapce (using FFT)
    !> \begin{equation}
    !>   P_{\alpha I, \alpha^{\prime} I^{\prime}}^{{\mathbf{q}}}\left(\tau\right)
    !>   =\sum_{\mathbf{R}} e^{-i \mathbf{q} \mathbf{R}}
    !>   P_{\alpha I, \alpha^{\prime} I^{\prime}}^{{\mathbf{R}}}\left(\tau\right)
    !> \end{equation}
    !!----------------------------------------------------------------------------------

    subroutine polarizability_q_mtmt(tau, pola_q_mtmt) 
        !> Imaginary time 
        real(wp), intent(in) :: tau 
        !> Polarizability in imaginary time and R representation in MT-MT
        complex(wp), allocatable, intent(out) :: pola_q_mtmt(:,:,:,:,:) 
        
        !> Local variables 
        integer :: ias, ias1, nlm, nlm1 
        !> Index for the Bravias vector R 
        integer :: ir 
        !> Total number of Bravias vector
        integer :: nbigR 
        !> Total number of n, l, m
        integer :: nlm_tot 
        !> Polarizability in imaginary time and R representation in MT-MT
        complex(wp), allocatable :: pola_R_mtmt(:,:,:,:,:) 
        !> Temporary array for FFT 
        complex(wp), allocatable :: zfft(:) 

        
        nbigR = kqset%nkpt 
        nlm_tot = locmatsiz/natmtot 
        
        allocate(pola_q_mtmt(nlm_tot, natmtot, nlm_tot, natmtot, kqset%nkpt)) 

        ! pola_q_mtmt = cmplx(0.0_wp, 0.0_wp, kind=wp) 
        pola_q_mtmt = zzero

        !> Polarizability in R-space in MT-MT
        call polarizability_R_mtmt(tau, pola_R_mtmt) 
        ! STOP   !! DELETE after TEST
          
        ! write(42,*) sum(real(pola_q_mtmt)), sum(imag(pola_q_mtmt))
        
        allocate(zfft(nbigR)) 
        print*, 'q grid', input%gw%ngridq 
        !print*, 'q grid', ngridq

        !! FFT from R to q space
        do ias = 1, natmtot
          do nlm = 1, nlm_tot
            do ias1 = 1, natmtot
              do nlm1 = 1, nlm_tot
                zfft(:) = pola_R_mtmt(nlm1, ias1, nlm, ias, :) 
                !!call zfftifc(3, input%gw%ngridq, -1, pola_q_mtmt(nlm1,ias1,nlm,ias,:))
                ! call zfftifc(3, input%gw%ngridq, -1, pola_q_mtmt(nlm1,ias1,nlm,ias,:))
                call zfftifc(3, input%gw%ngridq, -1, zfft)
                ! pola_q_mtmt(nlm1, ias1, nlm, ias, :) = zfft(:)
                pola_q_mtmt(nlm1, ias1, nlm, ias, :) = zfft(:) * kqset%nkpt   !! TODO:(Manoar) see below
                !! multiplying by "kqset%nkpt" gives the compareable polarizability with quartic 
                !! (but don't know why this so? because it is not present in the eqn)
              enddo
            enddo
          enddo
        enddo
        print*, 'sum pola_q mt-mt', sum(pola_q_mtmt)
        print*, 'test pola_q_mtmt', minval(real(pola_q_mtmt)), maxval(imag(pola_q_mtmt))
        
        write(43,*) sum(real(pola_q_mtmt)), sum(imag(pola_q_mtmt))
        do ir = 1, nbigR
          write(44,*) sum(real(pola_q_mtmt(:,:,:,:,ir))), sum(imag(pola_q_mtmt(:,:,:,:,ir)))
        enddo
        
        
        deallocate(pola_R_mtmt)
        deallocate(zfft)
        
    end subroutine polarizability_q_mtmt


    !!------------------------------------------------------------------------------------
    !! Transforamtion of the polarizability from $\tau$ to $\omega$ domain
    !> \begin{equation}
    !>   P_{\alpha I, \alpha^{\prime} I^{\prime}}^{{\mathbf{q}}}\left(\omega_i\right) = 
    !>   \sum_{j=1}^{n} \delta_{ij} cos(\omega_i \tau_j) P_{\alpha I, \alpha^{\prime} 
    !>   I^{\prime}}^{{\mathbf{q}}}\left(\tau_j\right)
    !> \end{equation}
    !!------------------------------------------------------------------------------------

    subroutine polarizability_tau_to_omega_mtmt(jt, weight_cos, pola_q_mtmt, pola_q_mtmt_omega)
        !> Index for $\tau$
        integer :: jt
        !> Weight and cos(w*t) for tau to omega transformation
        real(wp), intent(in) :: weight_cos(:,:)
        !> Polarizability in imaginary time and R representation in MT-MT
        ! complex(wp), allocatable, intent(in) :: pola_q_mtmt(:,:,:,:,:)
        complex(wp), intent(in) :: pola_q_mtmt(:,:,:,:,:)
        !> Polarizability in R representation and in omega domain (in MT-MT)
        complex(wp), intent(out) :: pola_q_mtmt_omega(:,:,:,:,:,:)

        !! Local variables
        integer :: iq, ias1, ias2, nlm1, nlm2, it, n_tau
        !> Total number of n, l, m
        integer :: nlm_tot

        ! complex(wp), allocatable :: pola_q_mtmt_omega(:,:,:,:,:)

        ! print*, 'shape of weight_cos', shape(weight_cos)
        ! print*, 'shape of pola_q_mtmt', shape(pola_q_mtmt)

        ! allocate(pola_q_mtmt_omega, mold=pola_q_mtmt)

        ! print*, 'shape of pola_q_mtmt_omega', shape(pola_q_mtmt_omega)
        ! print*, 'size of pola_q_mtmt_omega', size(pola_q_mtmt_omega, 1)
        ! ! allocate(pola_q_mtmt(nlm_tot, natmtot, nlm_tot, natmtot, kqset%nkpt))

        nlm_tot = locmatsiz/natmtot
        n_tau = size(weight_cos, 1)
        ! print*, ' no. of tau points, nlm_tot', n_tau, nlm_tot

        do it = 1, n_tau
          do iq = 1, kqset%nkpt
            do ias2 = 1, natmtot
              do nlm2 = 1, nlm_tot
                do ias1 = 1, natmtot
                  do nlm1 = 1, nlm_tot
                    pola_q_mtmt_omega(nlm1, ias1, nlm2, ias2, iq, it) = & 
                    pola_q_mtmt_omega(nlm1, ias1, nlm2, ias2, iq, it) + &
                    weight_cos(it, jt) * pola_q_mtmt(nlm1, ias1, nlm2, ias2, iq)
                    if(abs(pola_q_mtmt_omega(nlm1, ias1, nlm2, ias2, iq, it))>10.0_wp) then 
                      write(51,*) weight_cos(it, jt), pola_q_mtmt(nlm1, ias1, nlm2, ias2, iq)
                    endif 
                  enddo 
                enddo 
              enddo 
            enddo 
          enddo 
        enddo 

    end subroutine polarizability_tau_to_omega_mtmt


end module mod_polarizability_q_mtmt