module mod_green_R_ii
    use precision, only: wp
    use constants, only: zzero, zone, real_one
    use mod_lattice, only: omega
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_bands, only : nomax, numin, nstdf, eveck
    use mod_kpoint, only: nkpt
    use modgw, only: kset, kqset, Gkqset, Gset, Gkset
    use mod_green_band, only: green_band
    use modinput, only: input   !! This is for "input%gw%ngridq"

    !!! TEST and DELETE
    use modmain, only: ngrid, ngrtot
    use mod_kpointset, only: G_set 
    use modxs, only: fftmap_type
    
    implicit none
    type(fftmap_type) :: fftmap
    private
    public :: green_R_ii

    contains

    !!-------------------------------------------------------------------------------------------
    !!
    !> \begin{equation}
    !>   G^{\mathbf k}_{\mathbf G, \mathbf G'}(\tau) = \frac{1}{\Omega} \sum_{n} 
    !>   C^{n \mathbf k}_{\mathbf G} G^{\mathbf k}_{n}(\tau) C^{* n \mathbf k}_{\mathbf G'}
    !> \end{equation}
    !>
    !> Let's say    L = C^{n \mathbf k}_{\mathbf G}
    !>              M = G^{\mathbf k}_{n}(\tau)
    !>              R = C^{* n \mathbf k}_{\mathbf G'}
    !> Therefore, G = L * M * R
    !>
    !> First GF is Fourier trnasformed from G' to r'-space
    !> \begin{equation}
    !>   G^{\mathbf k}_{\mathbf G, \mathbf r'}(\tau) = \sum_{\mathbf G'} 
    !>   G^{\mathbf k}_{\mathbf G, \mathbf G'}(\tau) e^{-i(\mathbf k + \mathbf G')\mathbf r'} 
    !> \end{equation}
    !>
    !> Then again Fourier trnasformed from G to r-space
    !> \begin{equation}
    !>   G^{\mathbf k}_{\mathbf r, \mathbf r'}(\tau) = \sum_{\mathbf G} 
    !<   e^{i(\mathbf k + \mathbf G)\mathbf r} G^{\mathbf k}_{\mathbf G, \mathbf r'}(\tau) 
    !> \end{equation}
    !>
    !> And finally GF in R-space from k-space by FFT 
    !> \begin{equation}
    !>   G^{\mathbf R}_{\mathbf r, \mathbf r'}(\tau) = \frac{1}{N_{\mathbf k}} \sum_{\mathbf k} 
    !>   e^{i \mathbf{k} \mathbf{R}} G^{\mathbf k}_{\mathbf r, \mathbf r'}(\tau)
    !> \end{equation}
    !!-------------------------------------------------------------------------------------------

    subroutine green_R_ii(tau, green_R_occ, green_R_uno)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Green function in imaginary time and R representation occupied part
        complex(wp), allocatable, intent(out) :: green_R_occ(:,:,:)
        !> Green function in imaginary time and R representation unoccupied part
        complex(wp), allocatable, intent(out) :: green_R_uno(:,:,:)
        
        !> Local variables
        integer :: ik, ik_irreduced, ig, igk, igk1, igf, ngk
        integer :: ir_grid, ir_grid1
        !> One by omega(volume)
        real(wp) :: omega_inv
        real(wp), allocatable :: green_n(:)
        complex(wp), allocatable :: green_gg_occ(:,:)
        complex(wp), allocatable :: green_gg_uno(:,:)
        complex(wp), allocatable :: green_gr_occ(:,:)
        complex(wp), allocatable :: green_gr_uno(:,:)
        complex(wp), allocatable :: green_n_g(:,:)
        complex(wp), allocatable :: fft_occ(:), fft_uno(:)

        !> External routines for matrix matrix multiplication
        external :: zgemm

        ! Test and delete
        real :: t_i, t_f

        print*, 'real space grid from I-I ', Gset%ngrid, Gset%ngrtot, ngrid, ngrtot
        print*, 'fftmap_type ngrtot, ngrid', fftmap%ngrtot, fftmap%ngrid
        
        !> nstdf = # states used to calculate the dielectric function
        allocate(green_n(nstdf))
        if(allocated(eveck)) deallocate(eveck)
        !> nmatmax = maximum nmat over all k-points, nstfv = # first-variational states
        allocate(eveck(nmatmax, nstfv))
        !> Temporary array for FFT
        allocate(fft_occ(Gset%ngrtot))
        allocate(fft_uno(Gset%ngrtot))
        !> "green_R_occ/uno" allocation fails for large "Gset%ngrtot"
        allocate(green_R_occ(Gset%ngrtot, Gset%ngrtot, nkpt))
        allocate(green_R_uno(Gset%ngrtot, Gset%ngrtot, nkpt))

        
        do ik = 1, nkpt !kqset%nkpt
          
          ngk = Gkqset%ngk(1,ik)                           !! Number of G-vectors 
          ik_irreduced = kset%ik2ikp(ik)                   !! If symmetry is considered in ground-state
          if(input%groundstate%nosym) ik_irreduced = ik    !! If symmetry is not considered in ground-state
          !! TODO:(Manoar) nscf needs to be included here too

          ! print*, 'ngk, nmatmax ', ngk, nmatmax
          write(63,*) ik, ngk, nkpt
          
          allocate(green_n_g(nstfv, ngk))
          allocate(green_gg_occ(ngk, ngk))
          allocate(green_gg_uno(ngk, ngk))
          allocate(green_gr_occ(ngk, Gset%ngrtot))
          allocate(green_gr_uno(ngk, Gset%ngrtot))
          
          call green_band(tau, ik_irreduced, green_n)

          !> Collecting KS eigenvectors to perform expand_evec
          call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)

          
          !! M * R = MR = green_n_g
          do ig = 1, ngk
            green_n_g(1:nstdf, ig) = green_n(1:nstdf) * conjg(eveck(ig, 1:nstdf))
          enddo

          omega_inv = real_one / omega

          !-------------------------------------------------------------------------
          !! L * MR = G = green_gg_occ/uno
          !-------------------------------------------------------------------------
          !> For occupied part
          ! call zgemm('n', 'n', ngk, ngk, nomax, zone, eveck(1,1), nstdf, &
          call zgemm('n', 'n', ngk, ngk, nomax,  omega_inv, eveck(1,1), nstdf, &
                     green_n_g(1,1), nstdf,  zzero, green_gg_occ(1,1), ngk)

          !> For unoccupied part
          ! ! call zgemm('n', 'n', ngk, ngk, nstdf-nomax, zone, eveck(1,1), nstdf, &
          ! call zgemm('n', 'n', ngk, ngk, nstdf-nomax, zone, eveck(numin,1), nstdf, &
          call zgemm('n', 'n', ngk, ngk, nstdf-nomax,  omega_inv, eveck(numin,1), nstdf, &
                     green_n_g(numin,1), nstdf,  zzero, green_gg_uno(1,1), ngk)
          !-------------------------------------------------------------------------

          ! green_gg_occ = green_gg_occ / omega
          ! green_gg_uno = green_gg_uno / omega

          print*, 'min, max real gg occ', minval(real(green_gg_occ)), maxval(real(green_gg_occ))
          print*, 'min, max real gg uno', minval(real(green_gg_uno)), maxval(real(green_gg_uno))
          print*, 'sum gg', sum(real(green_gg_occ)), sum(imag(green_gg_occ))
          
          write(31,*) 'green_gg_occ', minval(real(green_gg_occ)), maxval(real(green_gg_occ))
          write(31,*) 'green_gg_occ', minval(imag(green_gg_occ)), maxval(imag(green_gg_occ))
          write(32,*) 'green_gg_uno', minval(real(green_gg_uno)), maxval(real(green_gg_uno))
          write(32,*) 'green_gg_uno', minval(imag(green_gg_uno)), maxval(imag(green_gg_uno))

          !! FFT from G' to r'-space of Green's function
          do igk = 1, ngk
            fft_occ = zzero
            fft_uno = zzero
            do igk1 = 1, ngk
              ig = Gkqset%igkig(igk1,1,ik)
              igf = Gset%igfft(ig)
              fft_occ(igf) = green_gg_occ(igk, igk1)
              fft_uno(igf) = green_gg_uno(igk, igk1)
              ! ! green_gr_occ(igk, igf) = green_gg_occ(igk, igk1)
            enddo

            call zfftifc(3, Gset%ngrid, -1, fft_occ)
            call zfftifc(3, Gset%ngrid, -1, fft_uno)
            ! ! call zfftifc(3, Gset%ngrid, -1, green_gr_occ(igk, :))
            !! TODO: do I need Gset%ngrid or Gkset%ngrid???

            green_gr_occ(igk, :) = fft_occ(:)
            green_gr_uno(igk, :) = fft_uno(:)
            ! green_gr_occ(igk, :) = fft_occ(:) * Gset%ngrtot
            ! green_gr_uno(igk, :) = fft_uno(:) * Gset%ngrtot
            !! TODO:(Manoar) multiplicatiojn by "Gset%ngrtot" is not present in the eqn but
            !! to match the pola I did it but don't know why !!!

          enddo

          write(31,*) 'green_gr_occ', minval(real(green_gr_occ)), maxval(real(green_gr_occ))
          write(31,*) 'green_gr_occ', minval(imag(green_gr_occ)), maxval(imag(green_gr_occ))
          write(32,*) 'green_gr_uno', minval(real(green_gr_uno)), maxval(real(green_gr_uno))
          write(32,*) 'green_gr_uno', minval(imag(green_gr_uno)), maxval(imag(green_gr_uno))

          print*, 'test green_k_occ', minval(real(green_gr_occ)), maxval(real(green_gr_occ))
          print*, 'test green_k_uno', minval(real(green_gr_uno)), maxval(real(green_gr_uno))
          print*, 'sum gr occ I-I 1', ik, sum(real(green_gr_occ)), sum(imag(green_gr_occ))
          print*, 'sum gr uno I-I 1', ik, sum(real(green_gr_uno)), sum(imag(green_gr_uno))

          call CPU_TIME(t_i)

          !! FFT from G to r-space of Green's function
          do ir_grid1 = 1, Gset%ngrtot
            fft_occ = zzero
            fft_uno = zzero
            do igk = 1, ngk
              ig = Gkqset%igkig(igk,1,ik)
              igf = Gset%igfft(ig)
              fft_occ(igf) = green_gr_occ(igk, ir_grid1)
              fft_uno(igf) = green_gr_uno(igk, ir_grid1)
              ! ! green_R_occ(igf, ir_grid1, ik) = green_gr_occ(igk, ir_grid1)
            enddo

            ! call zfftifc(3, Gset%ngrid, -1, fft_occ)   !!! Do we need "-1" or "1"; see copy for explanation
            ! call zfftifc(3, Gset%ngrid, -1, fft_uno)
            call zfftifc(3, Gset%ngrid, 1, fft_occ) 
            call zfftifc(3, Gset%ngrid, 1, fft_uno)
            !! TODO: do I need Gset%ngrid or Gkset%ngrid???
            ! ! call zfftifc(3, Gset%ngrid, -1, green_R_occ(:, ir_grid1, ik))

            ! green_R_occ(:, ir_grid1, ik) = fft_occ(:)
            ! green_R_uno(:, ir_grid1, ik) = fft_uno(:)
            green_R_occ(:, ir_grid1, ik) = fft_occ(:) * Gset%ngrtot
            green_R_uno(:, ir_grid1, ik) = fft_uno(:) * Gset%ngrtot
            !! TODO:(Manoar) multiplicatiojn by "Gset%ngrtot" is not present in the eqn but
            !! to match the pola I did it but don't know why !!!

          enddo

          write(31,*) 'green_rr_occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
          write(31,*) 'green_rr_occ', minval(imag(green_R_occ)), maxval(imag(green_R_occ))
          write(32,*) 'green_rr_uno', minval(real(green_R_uno)), maxval(real(green_R_uno))
          write(32,*) 'green_rr_uno', minval(imag(green_R_uno)), maxval(imag(green_R_uno))

          call CPU_TIME(t_f)
          print*, 'test green_k_occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
          print*, 'test green_k_uno', minval(real(green_R_uno)), maxval(real(green_R_uno))
          print*, 'time FFT green I-I 2', t_f - t_i, sum(real(green_R_occ)), sum(imag(green_R_occ))
          print*, 'time FFT green I-I 2', t_f - t_i, sum(real(green_R_uno)), sum(imag(green_R_uno))


          deallocate(green_n_g)
          deallocate(green_gg_occ)
          deallocate(green_gg_uno)
          deallocate(green_gr_occ)
          deallocate(green_gr_uno)
          
        enddo   ! k-loop ends 

        call CPU_TIME(t_i)
        !! FFT of the Green's function from k to R-space
        do ir_grid = 1, Gset%ngrtot
          do ir_grid1 = 1, Gset%ngrtot
            call zfftifc(3, input%gw%ngridq, 1, green_R_occ(ir_grid, ir_grid1, :))
            call zfftifc(3, input%gw%ngridq, 1, green_R_uno(ir_grid, ir_grid1, :))
          enddo
        enddo

        write(31,*) 'green_R_occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
        write(31,*) 'green_R_occ', minval(imag(green_R_occ)), maxval(imag(green_R_occ))
        write(32,*) 'green_R_uno', minval(real(green_R_uno)), maxval(real(green_R_uno))
        write(32,*) 'green_R_uno', minval(imag(green_R_uno)), maxval(imag(green_R_uno))

        call CPU_TIME(t_f)
        print*, 'time, sum Green_R_occ', t_f - t_i, sum(green_R_occ)
        print*, 'time, sum Green_R_uno', t_f - t_i, sum(green_R_uno)
        print*, 'test green_R_occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
        print*, 'test green_R_uno', minval(real(green_R_uno)), maxval(real(green_R_uno))


        deallocate(fft_occ)
        deallocate(fft_uno)
        deallocate(eveck)
        deallocate(green_n)

        
    end subroutine green_R_ii
    !
end module mod_green_R_ii
