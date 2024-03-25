module mod_polarizability_R_mti
    use precision, only: wp
    use modinput, only: input
    use mod_atoms, only: natmtot, idxas, natoms
    use mod_APW_LO, only: apwordmax, nlomax, lolmmax, nlorb, lorbl
    use modgw, only: kset, kqset, Gset
    use mod_muffin_tin, only: idxlm, lmmaxapw
    use mod_gaunt_coefficients, only: getgauntcoef
    use mod_product_basis, only: locmatsiz, mbindex, nmix, bigl, bradketa, bradketlo
    use index_contracting, only: nlm_index_contracting
    use mod_green_R_mti, only: green_R_mti
    use mod_overlap_mtmt, only: overlap_mtmt
    !! TEST and DELETE
    use mod_Gkvector, only: ngk, gkmax
    use modxs, only: fftmap_type


    implicit none
    private
    public :: polarizability_R_mti
    !! TEST and DELETE
    type(fftmap_type) :: fftmap
    

    contains

    !!---------------------------------------------------------------------------------------
    !! Polarizability in MT-I region in R-space
    !> \begin{equation}
    !>   P_{\alpha I, \mathbf{r}^{\prime}}^{{\mathbf{R}}}\left(\tau\right) = - 
    !>   \sum_{\xi \ell m} \sum_{\xi' \ell' m'} S^{\alpha I}_{\xi \ell m ; \xi' \ell' m'} 
    !>   G_{\alpha \xi \ell m ; \mathbf{r}'}^{\mathbf{R}}\left(\tau\right) 
    !>   G_{\alpha \xi' \ell' m' ; \mathbf{r}'}^{\mathbf{R}}\left(-\tau\right)
    !> \end{equation}
    !!---------------------------------------------------------------------------------------

    subroutine polarizability_R_mti(tau, pola_R_mti)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Polarizability in imaginary time and R representation in MT-I
        complex(wp), allocatable, intent(out) :: pola_R_mti(:,:,:,:)
        
        ! Local variables
        !> APW and LO combination variables
        integer :: i, j 
        integer :: imix
        integer :: ia, is, ias, ias1
        integer :: mbn, mbl, mbm
        integer :: l, m, lm, l1, m1, lm1, l2, m2, lm2
        integer :: l1min, l1max
        integer :: io, io1, io2
        integer :: ilo, ilo1 
        integer :: nlm, nlm1, nlm_tot
        !> Max of lm including APW and LO
        integer :: lm_lo_max
        !> Index for r-mesh or G-mesh
        integer :: ir_grid
        !> Running index for the Bravias vector R 
        integer :: ir
        !> Number of Bravais lattice R
        integer :: nbigR
        !> Gaunt coefficient
        real(wp) :: gaunt_coefficient
        real(wp), parameter :: tolerance = 0.00000001_wp
        !> Green function in imaginary time and R representation occupied part
        complex(wp), allocatable :: green_R_occ(:,:,:,:,:,:)
        !> Green function in imaginary time and R representation unoccupied part
        complex(wp), allocatable :: green_R_uno(:,:,:,:,:,:)
        !> Overlap matrix in MT-MT region
        real(wp), allocatable :: s_mtmt(:,:,:,:,:,:,:,:) 
        ! Combined NLM indices 
        integer, allocatable :: nlm_indices(:, :, :)
        
        ! Test and delete
        real :: t_i, t_f

        lm_lo_max = max(lmmaxapw, lolmmax)

        write(13,*) 'ngk test from mod_Gkvector', ngk
        !> ngk from mod_Gkvector and Gkqset%ngk from modgw are exactly same 
        call genfftmap(fftmap,2.01d0*gkmax)
        print*, '"gkmax" from mod_Gkvector', gkmax
        write(14,*) fftmap%ngrid, fftmap%ngrtot
        write(14,*) Gset%ngrid, Gset%ngrtot, Gset%ngvec
        
        nbigR = kqset%nkpt
        !!print*, '*******', kqset%nkpt, nbigR
        call init_kqpoint_set()  !!! TEST and delete  !! this print the "gkmax" and "gqmax" and "gqmaxbarc"
        
        call nlm_index_contracting(nlm_tot, nlm_indices)
        
        allocate(pola_R_mti(Gset%ngrtot, nlm_tot, natmtot, nbigR))
        !> Initialization of "pola_R_mtmt" array
        pola_R_mti = cmplx(0.0_wp, 0.0_wp, kind=wp)

        print*, 'nlm dimension test', lmmaxapw, nlm_tot, input%groundstate%lmaxapw
        
        call overlap_mtmt(nlm_indices, s_mtmt)
        
        call green_R_mti(tau, green_R_occ, green_R_uno)
        ! STOP   !! TEST NAD DELETE

        print*, 'sum greenRocc call', sum(real(green_R_occ)), sum(imag(green_R_occ))
        print*, 'sum greenRocc call', sum(real(green_R_uno)), sum(imag(green_R_uno))
        print*, 'real space r grid', Gset%ngrid, product(Gset%ngrid), Gset%ngrtot
        
        do ir = 1, nbigR
          
          call CPU_TIME(t_i)
          
          do ias = 1, natmtot
            do nlm = 1, nlm_tot
              is  = nlm_indices(nlm, ias, 1)
              mbn = nlm_indices(nlm, ias, 3)
              mbl = nlm_indices(nlm, ias, 4)
              mbm = nlm_indices(nlm, ias, 5)
              do l = 0, input%groundstate%lmaxapw
                do m = - l, l
                  lm = idxlm(l, m)
                  m1 = - mbm + m
                  l1min = abs(m1)
                  l1max = min(mbl+l, input%groundstate%lmaxapw)
                  do l1 = l1min, l1max 
                    gaunt_coefficient = getgauntcoef(l1, mbl, l, m1, mbm)
                    if(abs(gaunt_coefficient)>tolerance) then
                      lm1 = idxlm(l1, m1)

                      do io = 1, apwordmax
                        do io1 = 1, apwordmax
                          !! {APW-APW} (1st term)
                          do ir_grid = 1, Gset%ngrtot
                            pola_R_mti(ir_grid, nlm, ias, ir) = pola_R_mti(ir_grid, nlm, ias, ir) - &
                            bradketa(2, mbn, l, io, l1, io1, ias) * gaunt_coefficient * &
                            green_R_uno(ir_grid, io, lm, ias, ir, 1) * conjg(green_R_occ(ir_grid, io1, lm1, ias, ir, 1))
                          enddo
                        enddo 
                        io1 = 1
                        do ilo1 = 1, nlorb(is)
                          if(lorbl(ilo1, is)==l1) then
                            !! {APW-LO} (2nd term)
                            do ir_grid = 1, Gset%ngrtot
                              pola_R_mti(ir_grid, nlm, ias, ir) = pola_R_mti(ir_grid, nlm, ias, ir) - &
                              bradketa(3, mbn, l, io, ilo1, io1, ias) * gaunt_coefficient * &
                              green_R_uno(ir_grid, io, lm, ias, ir, 1) * conjg(green_R_occ(ir_grid, io1, lm1, ias, ir, 2))
                              !! lm1 --> ilo1
                            enddo 
                          endif 
                        enddo 
                      enddo 

                      io = 1
                      do ilo = 1, nlorb(is)
                        if(lorbl(ilo, is)==l) then
                          do io1 = 1, apwordmax
                            !! {LO-APW} (3rd term)
                            do ir_grid = 1, Gset%ngrtot
                              pola_R_mti(ir_grid, nlm, ias, ir) = pola_R_mti(ir_grid, nlm, ias, ir) - & 
                              bradketlo(2, mbn, ilo, l1, io1, ias) * gaunt_coefficient * & 
                              green_R_uno(ir_grid, io, lm, ias, ir, 2) * conjg(green_R_occ(ir_grid, io1, lm1, ias, ir, 1))
                              !! lm --> ilo
                            enddo 
                          enddo 
                          io1 = 1
                          do ilo1 = 1, nlorb(is)
                            if(lorbl(ilo1, is)==l1) then
                              !! {LO-LO} (4th term)
                              do ir_grid = 1, Gset%ngrtot
                                pola_R_mti(ir_grid, nlm, ias, ir) = pola_R_mti(ir_grid, nlm, ias, ir) - & 
                                bradketlo(3, mbn, ilo, ilo1, io1, ias) * gaunt_coefficient * & 
                                green_R_uno(ir_grid, io, lm, ias, ir, 2) * conjg(green_R_occ(ir_grid, io1, lm1, ias, ir, 2))
                                !! lm --> ilo & lm1 --> ilo1
                              enddo 
                            endif 
                          enddo 
                        endif 
                      enddo 

                    endif 
                  enddo 
                enddo 
              enddo 
            enddo 
          enddo 

          !! Spin
          pola_R_mti = 2.0_wp * pola_R_mti   !! spin included here but (TODO:) do it properly 

          call CPU_TIME(t_i)
          print*, 'sum pola_R_mti', sum(pola_R_mti(:,:,:,ir)), 'time', t_f - t_i
          write(80111,*) sum(pola_R_mti(:,:,:,ir)), ir

          ! !! old 
          ! call CPU_TIME(t_i)
          ! do i = 1, 2   !! APW/LO 
          !   do j = 1, 2   !! APW/LO 
          
          !     do ias = 1, natmtot
          !       do nlm = 1, nlm_tot
          !         do lm1 = 1, lm_lo_max
          !           do io1 = 1, apwordmax !apword(l1,is)   !! needs to fix these 
          !             do lm2 = 1, lm_lo_max
          !               do io2 = 1, apwordmax !apword(l2,is) !! needs to fix these
          !                 do ir_grid = 1, Gset%ngrtot

          !                   pola_R_mti(ir_grid, nlm, ias, ir) = pola_R_mti(ir_grid, nlm, ias, ir) - &
          !                   s_mtmt(io2, lm2, io1, lm1, nlm, ias, j, i) * &  
          !                   green_R_uno(ir_grid, io1, lm1, ias, ir, j) * &
          !                   conjg(green_R_occ(ir_grid, io2, lm2, ias, ir, j))
          !                   !! check if this needs complex conjugate of "s_mtmt" matrix 

          !                 enddo 
          !               enddo 
          !             enddo 
          !           enddo 
          !         enddo 
          !       enddo 
          !     enddo 
          !   enddo 
          ! enddo 
          ! call CPU_TIME(t_f)
          ! print*, 'sum pola_R_mti', sum(pola_R_mti(:,:,:,ir)), 'time', t_f - t_i
          ! write(80111,*) sum(pola_R_mti(:,:,:,ir)) 


        
        enddo   ! ir loop ends 

        ! pola_R_mti = cmplx(0.0_wp, 0.0_wp, kind=wp)  !! TEST AND DELETE
        ! do ir = 1, nbigR

        !   ! OLD --------------------------------------------
        !   call CPU_TIME(t_i)
        !   do ias = 1, natmtot
        !     do nlm = 1, nlm_tot
        !       do lm1 = 1, lmmaxapw
        !         do lm2 = 1, lmmaxapw
        !           do io1 = 1, 1 !apword(l1,is)   !! needs to fix these 
        !             do io2 = 1, 1 !apword(l2,is) !! needs to fix these
        !               do ir_grid = 1, Gset%ngrtot

        !                 pola_R_mti(ir_grid, nlm, ias, ir) = pola_R_mti(ir_grid, nlm, ias, ir) + &
        !                                                     s_mtmt(lm2, lm1, nlm, ias,1,1,1,1) * &  !! TODO: change this according origanl overlap matrix
        !                                                     green_R_uno(ir_grid, lm1, ias, ir,1,1) * &
        !                                                     conjg(green_R_occ(ir_grid, lm2, ias, ir,1,1))
        !                   !! check if this needs complex conjugate of "s_mtmt" matrix 
        !               enddo   !! ir_grid
        !             enddo   !! io2
        !           enddo   !! io1
        !         enddo   !! lm2
        !       enddo   !! lm1
        !     enddo   !! nlm
        !   enddo   !! ias
        !   call CPU_TIME(t_f)
        !   print*, 'sum pola_R_mti', sum(pola_R_mti(:,:,:,ir)), 'time', t_f - t_i
        !   write(91,*) sum(pola_R_mti(:,:,:,ir)) 
          
          
        ! enddo   ! ir loop ends
        print*, 'sum pola_R_mti', sum(pola_R_mti)

        ! do ir = 1, nbigR
        !  write(402,*) sum(real(pola_R_mtmt(:,:,:,:,ir))), sum(imag(pola_R_mtmt(:,:,:,:,ir))), ir
        ! enddo
        ! write(41,*) sum(real(pola_R_mtmt)), sum(imag(pola_R_mtmt))

        deallocate(green_R_occ)
        deallocate(green_R_uno)
        deallocate(s_mtmt)
        deallocate(nlm_indices)

    end subroutine polarizability_R_mti

end module mod_polarizability_R_mti
