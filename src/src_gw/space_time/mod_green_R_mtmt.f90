module mod_green_R_mtmt
    use precision, only: wp
    use modinput, only:  input  !! This is for "input%gw%ngridq"
    use constants, only: zzero, zone
    use mod_atoms, only: natmtot, nspecies, natoms, idxas 
    use mod_bands, only : nomax, numin, nstdf, eveck, eveckalm
    use mod_eigensystem, only: nmatmax, idxlo
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_muffin_tin, only: lmmaxapw, idxlm
    use mod_APW_LO, only: apwordmax, apword, nlorb, lorbl, nlomax, lolmmax, maxapword
    use modgw, only: kset, kqset, Gkqset, Gkset
    use mod_kpoint, only: ivk, ivknr, nkpt, nkptnr
    use mod_bravais_lattice, only: bravais_lattice
    use mod_green_band, only: green_band

    implicit none
    private
    public :: green_R_mtmt

    contains

    !!--------------------------------------------------------------------------------
    !> Green' function (GF) in k-space from the band-representation GF
    !> \begin{equation}
		!>  \begin{aligned}
		!>   G^{\mathbf k}_{\alpha \xi l m , \alpha\prime \xi\prime l' m'}(\tau) = 
    !>     \sum_{n} \mathcal{A}_{\alpha \xi l m}^{n\mathbf k} G^{\mathbf k}_{n}(\tau) 
    !>     \mathcal{A}_{\alpha\prime \xi\prime l' m'}^{* n \mathbf k}
		!>   \nonumber
		!>  \end{aligned}
	  !> \end{equation}
    !>
    !> Let's say    L = \mathcal{A}_{\alpha \xi l m}^{n\mathbf k}
    !>              M = G^{\mathbf k}_{n}(\tau)
    !>              R = \mathcal{A}_{\alpha\prime \xi\prime l' m'}^{* n \mathbf k}
    !> Therefore, G = L * M * R
    !> L and R can be APW and LO. Depending on this G will have 4 combinations 
    !> (1) APW-APW  (2) APW-LO  (3) LO-APW  (4) LO-LO
    !>
    !> Then GF in R-space from k-space (FFT)
    !> \begin{equation}
    !>   G^{\mathbf R}_{\alpha \xi l m , \alpha\prime \xi\prime l' m'}(\tau) = 
    !>   \frac{1}{N_{\mathbf k}} \sum_{\mathbf k} e^{i \mathbf k \mathbf R} 
    !>   G^{\mathbf k}_{\alpha \xi l m , \alpha\prime \xi\prime l' m'}(\tau) 
    !> \end{equation}
    !!--------------------------------------------------------------------------------

    subroutine green_R_mtmt(tau, green_R_occ, green_R_uno)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Green function in imaginary time and R representation occupied part
        complex(wp), allocatable, intent(out) :: green_R_occ(:,:,:,:,:,:,:,:)
        !> Green function in imaginary time and R representation unoccupied part
        complex(wp), allocatable, intent(out) :: green_R_uno(:,:,:,:,:,:,:,:)
        
        !> Local variables
        !> APW and LO combination variables
        integer :: i, j 
        integer :: ik, ik_irreduced
        integer :: io, io1, ilo, is, ia, ias, ias1, l, m, lm, l1, m1, lm1
        !> Combined dimension 
        integer :: n_dim
        !> Max of lm and LO
        integer :: lm_lo_max
        !> APW and LO combinations
        integer :: n_apw_lo
        real(wp), allocatable :: green_n(:)
        complex(wp), allocatable :: evecklo(:,:,:,:) 
        complex(wp), allocatable :: green_apw(:,:,:,:)
        complex(wp), allocatable :: green_lo(:,:,:,:)
        complex(wp), allocatable :: green_pp_occ(:,:,:,:,:,:,:)
        complex(wp), allocatable :: green_pp_uno(:,:,:,:,:,:,:)
        complex(wp), allocatable :: fft_occ(:), fft_uno(:)
        ! Test and delete
        real :: t_i, t_f
        
        !> External routines for matrix matrix multiplication
        external :: zgemm

        n_apw_lo = 4 
        ! nxi_max = max(input%groundstate%lmaxapw, nlomax)
        ! lm_max = max(lmmaxapw, lolmmax)
        lm_lo_max = max(lmmaxapw, nlomax)
        print*, 'lm_lo_max, lmmaxapw, lolmmax, nlomax', lm_lo_max, lmmaxapw, lolmmax, nlomax
        print*, 'apwordmax, maxapword', apwordmax, maxapword
        !! what is the diff b/t apwordmax and maxapword???


        !> nstdf = # states used to calculate the dielectric function
        allocate(green_n(nstdf))
        if(allocated(eveck)) deallocate(eveck)
        !> nmatmax = maximum nmat over all k-points, nstfv = # first-variational states
        allocate(eveck(nmatmax, nstfv))
        if(allocated(eveckalm)) deallocate(eveckalm)
        !> apwordmax = maximum of apword over all angular momenta and species
        !> lmmaxapw = maximum angular momentum for augmented plane waves
        !> lolmmax = maximum angular momentum for local orbitals
        !> lm_lo_max = maximum angular momentum out of APW & LO
        !> natmtot = totaal number of atoms 
        allocate(eveckalm(nstfv, apwordmax, lmmaxapw, natmtot))
        allocate(evecklo(nstfv, apwordmax, lm_lo_max, natmtot))
        allocate(green_apw(nstfv, apwordmax, lm_lo_max, natmtot))
        allocate(green_lo(nstfv, apwordmax, lm_lo_max, natmtot))
        allocate(green_pp_occ(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, natmtot, n_apw_lo))
        allocate(green_pp_uno(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, natmtot, n_apw_lo))
        allocate(green_R_occ(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, natmtot, nkpt, n_apw_lo))
        allocate(green_R_uno(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, natmtot, nkpt, n_apw_lo))
        allocate(fft_occ(nkpt), fft_uno(nkpt))

        green_R_occ = zzero
        green_R_uno = zzero
        
        
        !> Total q/k-points
        !! kqset%nkpt
        !> Irreducible k-points
        !! kset%nkpt
        ! print*, 'k points from G_R_mtmt routine', kset%nkpt, kqset%nkpt, nkpt, nkptnr
        ! print*, 'input%groundstate%lmaxapw, lmmaxapw', input%groundstate%lmaxapw, lmmaxapw
        
        !> k loop starts 
        ! TODO (Manoar): Do we need kqset or kset or nkpt or nkptnr (reduced or irreduced) ?
        print*, 'kqset or kset or nkpt or nkptnr', kqset%nkpt, kset%nkpt, nkpt, nkptnr


        do ik = 1, kqset%nkpt !nkpt !kqset%nkpt, kset%nkpt

          
          ik_irreduced = kset%ik2ikp(ik)                   !! If symmetry is considered in ground-state
          if(input%groundstate%nosym) ik_irreduced = ik    !! If symmetry is not considered in ground-state
          !! TODO:(Manoar) nscf needs to be included here too

          !! TEST and DELETE
          is = kset%ik2ikp(ik) 
          ia = kset%ikp2ik(ik)
          write(36,*) is, ia, ik, ik_irreduced

          ! call green_band(tau, ik, green_n)
          call green_band(tau, ik_irreduced, green_n)

          !> Collecting KS eigenvectors to perform expand_evec
          call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)
          ! TODO (Manoar): Do we need kqset or kset?

          !> Calculate the product of an eigenvector with the corresponding 
          !> matching coefficients ---> Compute products ==> \[ \sum_G C_{k}n * A_{lm} \]
          !> Result is written in eveckalm(ist,io,lm,ias) array
          call expand_evec(ik, 'T')


          print*, 'species atoms ', nspecies, maxval(natoms), maxval(nlorb), nstdf, nstfv
          print*, 'real minval(eveck), maxval(eveck)', minval(real(eveck)), maxval(real(eveck))
          print*, 'imag minval(eveck), maxval(eveck)', minval(imag(eveck)), maxval(imag(eveck))
          write(35,*) minval(real(eveck)), maxval(real(eveck))
          write(35,*) minval(imag(eveck)), maxval(imag(eveck))


          green_apw = zzero
          green_lo = zzero
          evecklo = zzero

          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)
              do l = 0, input%groundstate%lmaxapw
                do m = -l, l
                  lm = idxlm(l,m)
                  !! M * R_{APW} = MR_{APW} = green_apw
                  do io = 1, apword(l, is)
                    green_apw(1:nstdf, io, lm, ias) = green_n(1:nstdf) * conjg(eveckalm(1:nstdf, io, lm, ias))
                  enddo 
                  !! M * R_{LO} = MR_{LO} = green_lo
                  io = 1
                  do ilo = 1, nlorb(is)
                    if(lorbl(ilo, is)==l) then
                      ! evecklo(1:nstdf, io, lm, ias) = eveck(Gkqset%ngk(1,ik)+idxlo(lm, ilo, ias), 1:nstdf)
                      ! green_lo(1:nstdf, io, lm, ias) = green_n(1:nstdf) * conjg(evecklo(1:nstdf, io, lm, ias))
                      evecklo(1:nstdf, io, ilo, ias) = eveck(Gkqset%ngk(1,ik)+idxlo(lm, ilo, ias), 1:nstdf)
                      green_lo(1:nstdf, io, ilo, ias) = green_n(1:nstdf) * conjg(evecklo(1:nstdf, io, ilo, ias))
                      !! TEST and DELETE------------------------------------------------------------------------
                      do i = 1, nstdf
                        if(abs(green_lo(i, io, ilo, ias))>5.0_wp) then
                          write(33,'(5i5, 2es25.15e3)') ias, ilo, i, io, lm, green_lo(i, io, ilo, ias)
                          write(34,'(2es25.15e3)') evecklo(i, io, ilo, ias)
                        endif 
                      enddo 
                      !! TEST and DELETE----
                    endif
                  enddo 
                enddo 
              enddo 
            enddo 
          enddo

          print*, 'min max GF apw', minval(real(green_apw)), maxval(real(green_apw))
          print*, 'min max GF lo ', minval(real(green_lo)), maxval(real(green_lo))
          print*, 'sum green_apw, green_lo', sum(green_apw), sum(green_lo)
                

          green_pp_occ = zzero
          green_pp_uno = zzero


          !!------------------------------------------------------------------------------------------
          n_dim = apwordmax * lm_lo_max * natmtot 
          !! L_{APW} * MR_{APW} = G_{APW-APW} 
          !> For occupied part
          call zgemm('t', 'n', n_dim, n_dim, nomax,  zone, eveckalm(1,1,1,1), nstdf, &
                     green_apw(1,1,1,1), nstdf,  zzero, green_pp_occ(1,1,1,1,1,1,1), n_dim)

          !> For unoccupied part
          call zgemm('t', 'n', n_dim, n_dim, nstdf-nomax,  zone, eveckalm(numin,1,1,1), nstdf, &
                     green_apw(numin,1,1,1), nstdf,  zzero, green_pp_uno(1,1,1,1,1,1,1), n_dim)
          !!------------------------------------------------------------------------------------------
          !! L_{APW} * MR_{LO} = G_{APW-LO} 
          !> For occupied part
          call zgemm('t', 'n', n_dim, n_dim, nomax,  zone, eveckalm(1,1,1,1), nstdf, &
                     green_lo(1,1,1,1), nstdf,  zzero, green_pp_occ(1,1,1,1,1,1,2), n_dim)

          !> For unoccupied part
          call zgemm('t', 'n', n_dim, n_dim, nstdf-nomax,  zone, eveckalm(numin,1,1,1), nstdf, &
                     green_lo(numin,1,1,1), nstdf,  zzero, green_pp_uno(1,1,1,1,1,1,2), n_dim)
          !!------------------------------------------------------------------------------------------
          !! L_{LO} * MR_{APW} = G_{LO-APW} 
          !> For occupied part
          call zgemm('t', 'n', n_dim, n_dim, nomax,  zone, evecklo(1,1,1,1), nstdf, &
                     green_apw(1,1,1,1), nstdf,  zzero, green_pp_occ(1,1,1,1,1,1,3), n_dim)

          !> For unoccupied part
          call zgemm('t', 'n', n_dim, n_dim, nstdf-nomax,  zone, evecklo(numin,1,1,1), nstdf, &
                     green_apw(numin,1,1,1), nstdf,  zzero, green_pp_uno(1,1,1,1,1,1,3), n_dim)
          !!------------------------------------------------------------------------------------------
          !! L_{LO} * MR_{LO} = G_{LO-LO} 
          !> For occupied part
          call zgemm('t', 'n', n_dim, n_dim, nomax,  zone, evecklo(1,1,1,1), nstdf, &
                     green_lo(1,1,1,1), nstdf,  zzero, green_pp_occ(1,1,1,1,1,1,4), n_dim)

          !> For unoccupied part
          call zgemm('t', 'n', n_dim, n_dim, nstdf-nomax,  zone, evecklo(numin,1,1,1), nstdf, &
                     green_lo(numin,1,1,1), nstdf,  zzero, green_pp_uno(1,1,1,1,1,1,4), n_dim)
          !!------------------------------------------------------------------------------------------


          !! TEST and DELETE -----------------------------------------------------------------------------
          write(371,*) 'min, max real occ', minval(real(green_pp_occ)), maxval(real(green_pp_occ))
          write(371,*) 'min, max imag occ', minval(imag(green_pp_occ)), maxval(imag(green_pp_occ))
          write(372,*) 'min, max real uno', minval(real(green_pp_uno)), maxval(real(green_pp_uno))
          write(372,*) 'min, max imag uno', minval(imag(green_pp_uno)), maxval(imag(green_pp_uno))
          print*, 'min, max real gf mtmt occ', minval(real(green_pp_occ)), maxval(real(green_pp_occ))
          print*, 'min, max imag gf mtmt occ', minval(imag(green_pp_occ)), maxval(imag(green_pp_occ))
          print*, 'min, max real gf mtmt uno', minval(real(green_pp_uno)), maxval(real(green_pp_uno))
          print*, 'min, max imag gf mtmt uno', minval(imag(green_pp_uno)), maxval(imag(green_pp_uno))
          print*, 'shape(green_pp_uno)', shape(green_pp_uno)
          ! allocate(green_pp_uno(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, natmtot, n_apw_lo))
          do i = 1, n_apw_lo
            do ias = 1, natmtot
              do l = 1, lm_lo_max
                do io = 1, apwordmax
                  do ias1 = 1, natmtot
                    do m = 1, lm_lo_max
                      do io1 = 1, apwordmax
                        write(37,*) green_pp_uno(io1, m, ias1, io, l, ias, i)
                        if(abs(green_pp_uno(io1, m, ias1, io, l, ias, i)) > 50.0_wp) write(32, '(7i4, 2es25.15e3)') &
                        i, ias, l, io, ias1, m, io1, green_pp_uno(io1, m, ias1, io, l, ias, i)
                      enddo 
                    enddo 
                  enddo 
                enddo 
              enddo 
            enddo 
          enddo 
          !! TEST and DELETE -----------------------------------------------------------------------------
          
          green_R_occ(:,:,:,:,:,:,ik,:) = green_pp_occ(:,:,:,:,:,:,:)
          green_R_uno(:,:,:,:,:,:,ik,:) = green_pp_uno(:,:,:,:,:,:,:)

        enddo !! ik (i.e., k) loop ends


        write(373,*) 'min, max real occ', minval(real(green_pp_occ)), maxval(real(green_pp_occ))
        write(373,*) 'min, max imag occ', minval(imag(green_pp_occ)), maxval(imag(green_pp_occ))
        write(374,*) 'min, max real uno', minval(real(green_pp_uno)), maxval(real(green_pp_uno))
        write(374,*) 'min, max imag uno', minval(imag(green_pp_uno)), maxval(imag(green_pp_uno))
        ! STOP   !! TEST and DELETE

        call CPU_TIME(t_i)

        !! FFT from k to R space of Green's function
        j = 0 
        do i = 1, n_apw_lo   !! APW/LO 
          do ias = 1, natmtot
            do lm = 1, lm_lo_max
              do io = 1, apwordmax
                do ias1 = 1, natmtot
                  do lm1 = 1, lm_lo_max
                    do io1 = 1, apwordmax
                      j = j + 1
                      ! call zfftifc(3, input%gw%ngridq, 1, green_R_occ(io1, lm1, ias1, io, lm, ias, :, i))
                      ! call zfftifc(3, input%gw%ngridq, 1, green_R_uno(io1, lm1, ias1, io, lm, ias, :, i))
                      fft_occ(:) = green_R_occ(io1, lm1, ias1, io, lm, ias, :, i)
                      fft_uno(:) = green_R_uno(io1, lm1, ias1, io, lm, ias, :, i)
                      if(j==1) write(61,*) fft_occ
                      if(j==1) write(62,*) fft_uno
                      if(j==1) write(61,*) '------'
                      if(j==1) write(62,*) '------'
                      call zfftifc(3, input%gw%ngridq, 1, fft_occ)
                      call zfftifc(3, input%gw%ngridq, 1, fft_uno)
                      if(j==1) write(61,*) fft_occ
                      if(j==1) write(62,*) fft_uno
                      green_R_occ(io1, lm1, ias1, io, lm, ias, :, i) = fft_occ(:)
                      green_R_uno(io1, lm1, ias1, io, lm, ias, :, i) = fft_uno(:)
                    enddo 
                  enddo 
                enddo 
              enddo 
            enddo
          enddo
        enddo
        green_R_occ = green_R_occ / nkpt
        green_R_uno = green_R_uno / nkpt

        call CPU_TIME(t_f)
        print*, 'time of FFT GF_R and sum', t_f - t_i, sum(green_R_occ), sum(green_R_uno)
        print*, 'test green_R_occ', minval(real(green_R_occ)), maxval(real(green_R_occ))
        print*, 'test green_R_uno', minval(real(green_R_uno)), maxval(real(green_R_uno))


        deallocate(fft_occ)
        deallocate(fft_uno)
        deallocate(eveck)
        deallocate(eveckalm)
        deallocate(evecklo)
        deallocate(green_n)
        deallocate(green_apw)
        deallocate(green_lo)
        deallocate(green_pp_occ)
        deallocate(green_pp_uno)
        
      end subroutine green_R_mtmt

end module mod_green_R_mtmt
