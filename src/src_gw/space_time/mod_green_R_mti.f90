module mod_green_R_mti
    use precision, only: wp
    use constants, only: zzero, zone, real_one
    use modinput, only: input   !! This is for "input%gw%ngridq"
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_lattice, only: omega
    use mod_bands, only : nomax, numin, nstdf, eveck, eveckalm
    use mod_APW_LO, only: apwordmax, apword, nlorb, lorbl, nlomax, lolmmax, maxapword
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_eigensystem, only: nmatmax, idxlo
    use mod_muffin_tin, only: lmmaxapw, idxlm
    use mod_kpoint, only: ivk, ivknr, nkpt
    use modgw, only: kset, kqset, Gkqset, Gkset, Gset
    use mod_bravais_lattice, only: bravais_lattice
    use mod_green_band, only: green_band

    !! TEST and delete
    use modmain, only : ngkmax
    use mod_coulomb_potential, only : barc
    use mod_overlap_ii, only: overlap_ii

    implicit none
    private
    public :: green_R_mti

    contains

    !!--------------------------------------------------------------------------------
    !> Green' function (GF) in k-space from the band-representation GF
    !> \begin{equation}
		!>  \begin{aligned}
		!>   G^{\mathbf k}_{\alpha \xi l m , \mathbf{G}'}(\tau) = \frac{1}{\sqrt{\Omega_0}}
    !>     \sum_{n} \mathcal{A}_{\alpha \xi l m}^{n\mathbf k} G^{\mathbf k}_{n}(\tau) 
    !>     C_{\mathbf{G}'}^{* n \mathbf k}
		!>   \nonumber
		!>  \end{aligned}
	  !> \end{equation}
    !>
    !> Let's say    L = \mathcal{A}_{\alpha \xi l m}^{n\mathbf k}
    !>              M = G^{\mathbf k}_{n}(\tau)
    !>              R = C_{\mathbf{G}'}^{* n \mathbf k}
    !> Therefore, G = L * M * R
    !> L can be APW as well as LO. Depending on this G will have 2 combinations 
    !> (1) APW-contribution  (2) LO-contribution
    !>
    !> Then Fourier trnasform of GF from G-space to r-space
    !> \begin{equation}
    !>   G^{\mathbf k}_{\alpha \xi l m , \mathbf{r}'}(\tau) = 
    !>   \sum_{\mathbf{G}'} e^{-i(\mathbf{k}+ \mathbf{G}')\mathbf{r}'} 
    !>   G^{\mathbf k}_{\alpha \xi l m , \mathbf{G}'}(\tau)
    !> \end{equation}
    !>
    !> And finally GF in R-space from k-space by FFT 
    !> \begin{equation}
    !>   G^{\mathbf R}_{\alpha \xi l m , \mathbf{r}'}(\tau) = 
    !>   \frac{1}{N_{\mathbf k}} \sum_{\mathbf{k}} e^{i \mathbf{k} \mathbf{R}} 
    !>   G^{\mathbf k}_{\alpha \xi l m , \mathbf{r}'}(\tau)
    !> \end{equation}
    !!--------------------------------------------------------------------------------

    subroutine green_R_mti(tau, green_R_occ, green_R_uno)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Green function in imaginary time and R representation occupied part
        complex(wp), allocatable, intent(out) :: green_R_occ(:,:,:,:,:,:)
        !> Green function in imaginary time and R representation unoccupied part
        complex(wp), allocatable, intent(out) :: green_R_uno(:,:,:,:,:,:)
        
        !> Local variables
        integer :: ik, ik_irreduced
        integer :: ig, igk, igf, ngk
        integer :: io, ilo, is, ia, ias, l, m, lm 
        !> APW and LO combination variables
        integer :: i
        !> Index for lattice vector R
        integer :: ir
        !> index for r-mesh 
        integer :: ir_grid
        !> Max of lm including APW and LO
        integer :: lm_lo_max
        !> Combined dimension 
        integer :: n_dim
        !> APW and LO combinations
        integer :: n_apw_lo
        !> One by square root of volume (omega)
        real(wp) :: sqrt_omega_inv
        real(wp) :: bigR(3)
        real(wp), allocatable :: green_n(:)
        complex(wp), allocatable :: green_n_g(:,:)
        complex(wp), allocatable :: evecklo(:,:,:,:)
        complex(wp), allocatable :: green_pg_occ(:,:,:,:,:)
        complex(wp), allocatable :: green_pg_uno(:,:,:,:,:)
        complex(wp), allocatable :: fft_occ(:), fft_uno(:)
        !! Test and delete--------------------------
        ! Test and delete
        real :: t_i, t_f
        !! Test and delete--------------------------
        
        !> External routines for matrix matrix multiplication
        external :: zgemm

        n_apw_lo = 2
        ! lm_max = max(lmmaxapw, lolmmax)
        lm_lo_max = max(lmmaxapw, nlomax)

        !> nstdf = # states used to calculate the dielectric function
        allocate(green_n(nstdf))
        if(allocated(eveck)) deallocate(eveck)
        !> nmatmax = maximum nmat over all k-points, nstfv = # first-variational states
        allocate(eveck(nmatmax, nstfv))
        if(allocated(eveckalm)) deallocate(eveckalm)
        !> apwordmax = maximum of apword over all angular momenta and species
        !> lmmaxapw = maximum angular momentum for augmented plane waves
        !> natmtot = total number of atoms 
        allocate(eveckalm(nstfv, apwordmax, lmmaxapw, natmtot))
        allocate(evecklo(nstfv, apwordmax, lm_lo_max, natmtot))
        allocate(fft_occ(Gset%ngrtot))
        allocate(fft_uno(Gset%ngrtot))
        allocate(green_R_occ(Gset%ngrtot, apwordmax, lm_lo_max, natmtot, nkpt, n_apw_lo))
        allocate(green_R_uno(Gset%ngrtot, apwordmax, lm_lo_max, natmtot, nkpt, n_apw_lo))

        green_R_occ = zzero
        green_R_uno = zzero
        
        !> Total q/k-points
        !! kqset%nkpt
        !> Irreducible k-points
        !! kset%nkpt
        
        ! !> Calculates the Bravais lattice (or translational vector) R
        ! call bravais_lattice(ir, bigR)
        ! ! write(36, '(4i3, 3f16.11)') ir, ivk(:,ir), bigR
        
        print*, 'nstfv, nstdf (MT-I)', nstfv, nstdf
        write(13,*) 'Gkqset%ngk test from modgw->Gkqset', Gkqset%ngk
        !> ngk from mod_Gkvector and Gkqset%ngk from modgw are exactly same

        

        !> k loop starts 
        ! TODO (Manoar): if we need kqset or kset or nkpt or nkptnr (reduced or irreduced) ?
        do ik = 1, nkpt !kqset%nkpt  !! k-loop starts

          ngk = Gkqset%ngk(1,ik)                           !! Number of G-vectors 
          ik_irreduced = kset%ik2ikp(ik)                   !! If symmetry is considered in ground-state
          if(input%groundstate%nosym) ik_irreduced = ik    !! If symmetry is not considered in ground-state
          !! TODO:(Manoar) nscf needs to be included here too

          print*, 'Gk ==> Gkqset%ngk(1,ik)', ik, Gkqset%ngk(1,ik), ngk
          
          allocate(green_n_g(nstfv, ngk))
          allocate(green_pg_occ(apwordmax, lm_lo_max, natmtot, ngk, n_apw_lo))
          allocate(green_pg_uno(apwordmax, lm_lo_max, natmtot, ngk, n_apw_lo))


          
          call green_band(tau, ik_irreduced, green_n)

          !> Collecting KS eigenvectors to perform expand_evec
          call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)

          print*, ' ++> shape eveck', shape(eveck), 'size', size(eveck), 'ngkmax', ngkmax
          do igk = 1, nstfv
            do ig = 1, nmatmax
              write(112,'(2i4,2es25.15e3)') ig, igk, real(eveck(ig, igk)), imag(eveck(ig, igk))
            enddo
          enddo

          !> Calculate the product of an eigenvector with the corresponding 
          !> matching coefficients ---> Compute products ==> \[ \sum_G C_{k}n * A_{lm} \]
          call expand_evec(ik, 'T')
          !> Result is written in eveckalm(ist,io,lm,ias) array

          evecklo = zzero
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)
              do l = 0, input%groundstate%lmaxapw
                do m = -l, l
                  lm = idxlm(l,m)
                  io = 1
                  do ilo = 1, nlorb(is)
                    if(lorbl(ilo, is)==l) then
                      !! Collecting eigen-vectors of LO in array "evecklo"
                      ! evecklo(1:nstdf, io, lm, ias) = eveck(Gkqset%ngk(1,ik)+idxlo(lm, ilo, ias), 1:nstdf)
                      ! evecklo(1:nstdf, io, lm, ias) = eveck(ngk+idxlo(lm, ilo, ias), 1:nstdf)
                      evecklo(1:nstdf, io, ilo, ias) = eveck(ngk+idxlo(lm, ilo, ias), 1:nstdf)
                    endif 
                  enddo 
                enddo 
              enddo 
            enddo 
          enddo

          !! M * R = MR = green_n_g
          green_n_g = zzero
          do ig = 1, ngk
            green_n_g(1:nstdf, ig) = green_n(1:nstdf) * conjg(eveck(ig, 1:nstdf))
          enddo

          print*, 'sum tmp green', ik, sum(real(green_n_g)), sum(imag(green_n_g))
          print*, 'sum tmp green', ik, sum(real(green_n_g(1:numin,:))), sum(imag(green_n_g(1:numin,:)))
          print*, 'sum tmp green', ik, sum(real(green_n_g(numin:,:))), sum(imag(green_n_g(numin:,:)))
          
          write(65,*) '-------------------------------------',  ik, tau
          write(65,*) 'min, max real tg', minval(real(green_n_g)), maxval(real(green_n_g))
          write(65,*) 'min, max imag tg', minval(imag(green_n_g)), maxval(imag(green_n_g))
          write(65,*) 'min, max real eveckalm', minval(real(eveckalm)), maxval(real(eveckalm))
          write(65,*) 'min, max imag eveckalm', minval(imag(eveckalm)), maxval(imag(eveckalm))
          print*, 'shape eveckalm', shape(eveckalm), 'size', size(eveckalm)
          print*, 'sum eveckalm', sum(real(eveckalm)), sum(imag(eveckalm))

          
          green_pg_occ = zzero
          green_pg_uno = zzero

          sqrt_omega_inv = real_one / dsqrt(omega)

          n_dim = apwordmax * lm_lo_max * natmtot
          !!------------------------------------------------------------------------------------------
          !! L_{APW} * MR = G_{APW} 
          !> For occupied part
          call zgemm('t', 'n', n_dim, ngk, nomax,  sqrt_omega_inv, eveckalm(1,1,1,1), nstdf, &
          ! call zgemm('t', 'n', n_dim, ngk, nomax,  zone, eveckalm(1,1,1,1), nstdf, &
                     green_n_g(1,1), nstdf,  zzero, green_pg_occ(1,1,1,1,1), n_dim)

          !> For unoccupied part
          call zgemm('t', 'n', n_dim, ngk, nstdf-nomax,  sqrt_omega_inv, eveckalm(numin,1,1,1), nstdf, &
          ! call zgemm('t', 'n', n_dim, ngk, nstdf-nomax,  zone, eveckalm(numin,1,1,1), nstdf, &
                     green_n_g(numin,1), nstdf,  zzero, green_pg_uno(1,1,1,1,1), n_dim)
          !!------------------------------------------------------------------------------------------
          !! L_{LO} * MR = G_{LO} 
          !> For occupied part
          call zgemm('t', 'n', n_dim, ngk, nomax,  sqrt_omega_inv, evecklo(1,1,1,1), nstdf, &
          ! call zgemm('t', 'n', n_dim, ngk, nomax,  zone, evecklo(1,1,1,1), nstdf, &
                     green_n_g(1,1), nstdf,  zzero, green_pg_occ(1,1,1,1,2), n_dim)

          !> For unoccupied part
          call zgemm('t', 'n', n_dim, ngk, nstdf-nomax,  sqrt_omega_inv, evecklo(numin,1,1,1), nstdf, &
          ! call zgemm('t', 'n', n_dim, ngk, nstdf-nomax,  zone, evecklo(numin,1,1,1), nstdf, &
                     green_n_g(numin,1), nstdf,  zzero, green_pg_uno(1,1,1,1,2), n_dim)
          !!------------------------------------------------------------------------------------------

          ! green_pg_occ = green_pg_occ / dsqrt(omega)
          ! green_pg_uno = green_pg_uno / dsqrt(omega)
          ! ! green_pg_occ = green_pg_occ 
          ! ! green_pg_uno = green_pg_uno 

          print*, 'omega and sqrt(omega)', omega, dsqrt(omega)

          write(330,*) 'real pgocc', minval(real(green_pg_occ)), maxval(real(green_pg_occ))
          write(330,*) 'imag pgocc', minval(imag(green_pg_occ)), maxval(imag(green_pg_occ))
          write(340,*) 'real pguno', minval(real(green_pg_uno)), maxval(real(green_pg_uno))
          write(340,*) 'imag pguno', minval(imag(green_pg_uno)), maxval(imag(green_pg_uno))
          ! write(330,*) 'real pgocc', minval(real(green_pg_occ(:,:,:,:,1))), maxval(real(green_pg_occ(:,:,:,:,1)))
          ! write(330,*) 'imag pgocc', minval(imag(green_pg_occ(:,:,:,:,1))), maxval(imag(green_pg_occ(:,:,:,:,1)))
          ! write(340,*) 'real pguno', minval(real(green_pg_uno(:,:,:,:,1))), maxval(real(green_pg_uno(:,:,:,:,1)))
          ! write(340,*) 'imag pguno', minval(imag(green_pg_uno(:,:,:,:,1))), maxval(imag(green_pg_uno(:,:,:,:,1)))
          ! write(3301,*) 'real pgocc', minval(real(green_pg_occ(:,:,:,:,2))), maxval(real(green_pg_occ(:,:,:,:,2)))
          ! write(3301,*) 'imag pgocc', minval(imag(green_pg_occ(:,:,:,:,2))), maxval(imag(green_pg_occ(:,:,:,:,2)))
          ! write(3401,*) 'real pguno', minval(real(green_pg_uno(:,:,:,:,2))), maxval(real(green_pg_uno(:,:,:,:,2)))
          ! write(3401,*) 'imag pguno', minval(imag(green_pg_uno(:,:,:,:,2))), maxval(imag(green_pg_uno(:,:,:,:,2)))


          !! FFT from G' to r'-space of Green's function
          do i = 1, n_apw_lo   !! APW/LO
            do ias = 1, natmtot
              do lm = 1, lm_lo_max
                do io = 1, apwordmax
                  fft_occ = zzero
                  fft_uno = zzero
                  do igk = 1, ngk
                    ig = Gkqset%igkig(igk,1,ik)
                    igf = Gset%igfft(ig)
                    fft_occ(igf) = green_pg_occ(io, lm, ias, igk, i)
                    fft_uno(igf) = green_pg_uno(io, lm, ias, igk, i)
                  enddo

                  call zfftifc(3, Gset%ngrid, -1, fft_occ)
                  call zfftifc(3, Gset%ngrid, -1, fft_uno)

                  ! green_R_occ(:, io, lm ,ias, ik, i) = fft_occ(:) 
                  ! green_R_uno(:, io, lm ,ias, ik, i) = fft_uno(:)
                  green_R_occ(:, io, lm ,ias, ik, i) = fft_occ(:) * Gset%ngrtot
                  green_R_uno(:, io, lm ,ias, ik, i) = fft_uno(:) * Gset%ngrtot
                  !! TODO:(Manoar) multiplicatiojn by "Gset%ngrtot" is not present in the eqn but
                  !! to match the pola I did it but don't know why !!!

                enddo 
              enddo
            enddo
          enddo 

          ! write(71010,*) 'real procc', minval(real(green_R_occ(:,:,:,:,ik,3))), maxval(real(green_R_occ(:,:,:,:,ik,3)))
          ! write(71010,*) 'imag procc', minval(imag(green_R_occ(:,:,:,:,ik,3))), maxval(imag(green_R_occ(:,:,:,:,ik,3)))
          ! write(71020,*) 'real pruno', minval(real(green_R_uno(:,:,:,:,ik,3))), maxval(real(green_R_uno(:,:,:,:,ik,3)))
          ! write(71020,*) 'imag pruno', minval(imag(green_R_uno(:,:,:,:,ik,3))), maxval(imag(green_R_uno(:,:,:,:,ik,3)))

          deallocate(green_n_g)
          deallocate(green_pg_occ)
          deallocate(green_pg_uno)

          
        enddo   ! k-loop ends 

        write(331,*) 'real occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
        write(331,*) 'imag occ', minval(imag(green_R_occ)), maxval(imag(green_R_occ))
        write(341,*) 'real uno', minval(real(green_R_uno)), maxval(real(green_R_uno))
        write(341,*) 'imag uno', minval(imag(green_R_uno)), maxval(imag(green_R_uno))


        ! print*, 'omega', omega
        print*, 'r grid', Gset%ngrid, Gset%ngrtot


        ! deallocate(fft_occ)
        ! deallocate(fft_uno)
        ! allocate(fft_occ(nkpt))
        ! allocate(fft_uno(nkpt))

        call CPU_TIME(t_i)

        !! FFT of the Green's function from k to R-space 
        do i = 1, n_apw_lo    !! APW/LO
          do ias = 1, natmtot
            do lm = 1, lm_lo_max
              do io = 1, apwordmax
                do ir_grid = 1, Gset%ngrtot

                  ! fft_occ(:) = green_R_occ(ir_grid, io, lm, ias, :, i)
                  ! fft_uno(:) = green_R_uno(ir_grid, io, lm, ias, :, i)

                  call zfftifc(3, input%gw%ngridq, 1, green_R_occ(ir_grid, io, lm, ias, :, i))
                  call zfftifc(3, input%gw%ngridq, 1, green_R_uno(ir_grid, io, lm, ias, :, i))

                  ! call zfftifc(3, input%gw%ngridq, 1, fft_occ)
                  ! call zfftifc(3, input%gw%ngridq, 1, fft_uno)

                  ! green_R_occ(ir_grid, io, lm, ias, :, i) = fft_occ(:)
                  ! green_R_uno(ir_grid, io, lm, ias, :, i) = fft_uno(:)

                enddo 
              enddo 
            enddo
          enddo
        enddo
        call CPU_TIME(t_f)
        write(81010,*) 'min max real occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
        write(81010,*) 'min max imag occ', minval(imag(green_R_occ)), maxval(imag(green_R_occ))
        write(81011,*) 'min max real uno', minval(real(green_R_uno)), maxval(real(green_R_uno))
        write(81011,*) 'min max imag uno', minval(imag(green_R_uno)), maxval(imag(green_R_uno))
        print*, 'time FFT I-I greenR', t_f-t_i
        print*, 'sum greenR occ', sum(real(green_R_occ(:,:,:,:,:,1))), sum(imag(green_R_occ(:,:,:,:,:,1)))
        print*, 'sum greenR uno', sum(real(green_R_uno(:,:,:,:,:,1))), sum(imag(green_R_uno(:,:,:,:,:,1)))


        green_R_occ = green_R_occ / nkpt
        green_R_uno = green_R_uno / nkpt


        deallocate(fft_occ)
        deallocate(fft_uno)
        deallocate(eveck)
        deallocate(evecklo)
        deallocate(eveckalm)
        deallocate(green_n)

    end subroutine green_R_mti
end module mod_green_R_mti
