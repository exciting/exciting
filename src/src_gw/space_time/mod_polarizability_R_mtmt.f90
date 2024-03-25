module mod_polarizability_R_mtmt
    use precision, only: wp
    use modinput, only: input
    use mod_atoms, only: nspecies, natmtot, idxas, natoms
    use mod_APW_LO, only: apwordmax, nlomax, lolmmax, nlorb, lorbl
    use modgw, only: kset, kqset
    use mod_muffin_tin, only: idxlm, lmmaxapw
    use mod_gaunt_coefficients, only: getgauntcoef
    use mod_product_basis, only: locmatsiz, mbindex, nmix, bigl, bradketa, bradketlo
    use index_contracting, only: nlm_index_contracting
    use mod_green_R_mtmt, only: green_R_mtmt
    use mod_overlap_mtmt, only: overlap_mtmt

    implicit none
    private
    public :: polarizability_R_mtmt

    contains

    !!-----------------------------------------------------------------------------------------------------
    !! Polarizability in MT-MT region in R-space
    !> \begin{equation}
    !>   P_{\alpha I, \alpha^{\prime} I^{\prime}}^{{\mathbf{R}}}\left(\tau\right) = -
    !>   \sum_{\xi \ell m} \sum_{\xi'' \ell'' m''} S^{\alpha I}_{\xi \ell m ; \xi'' \ell'' m''} 
    !>   \sum_{\xi' \ell' m'} G_{\alpha \xi \ell m ; \alpha' \xi' \ell' m'}^{\mathbf{R}}\left(\tau\right) 
    !>   \sum_{\xi''' \ell''' m'''} G_{\alpha \xi'' \ell'' m'' ; \alpha' \xi''' \ell''' m'''}^
    !>   {\mathbf{R}}\left(-\tau\right) S^{\alpha' I'^*}_{\xi' \ell' m' ; \xi''' \ell''' m'''} 
    !> \end{equation}
    !!-----------------------------------------------------------------------------------------------------

    subroutine polarizability_R_mtmt(tau, pola_R_mtmt)
        !> Imaginary time 
        real(wp), intent(in) :: tau
        !> Polarizability in imaginary time and R representation in MT-MT
        complex(wp), allocatable, intent(out) :: pola_R_mtmt(:,:,:,:,:)
        
        ! Local variables
        !> APW and LO combination variables
        integer :: i, j 
        integer :: is, ia, ias, is1, ias1
        integer :: l, l1, l2, l3
        integer :: m, m1, m2, m3
        integer :: mbn, mbl, mbm
        integer :: mbn1, mbl1, mbm1
        integer :: lm, lm1, lm2, lm3 
        integer :: io, io1, io2, io3
        integer :: ilo, ilo1, ilo2, ilo3
        integer :: nlm, nlm1, nlm_tot
        integer :: l3min, l3max
        integer :: l2min, l2max
        !> Max of lm including APW and LO
        integer :: lm_lo_max
        !> Running index for the Bravias vector R 
        integer :: ir
        !> Number of Bravais lattice R
        integer :: nbigR
        !> Gaunt coefficient 
        real(wp) :: gaunt_coefficient
        real(wp), parameter :: tolerance = 0.00000001_wp
        !> Green function in imaginary time and R representation occupied part
        complex(wp), allocatable :: green_R_occ(:,:,:,:,:,:,:,:)
        !> Green function in imaginary time and R representation unoccupied part
        complex(wp), allocatable :: green_R_uno(:,:,:,:,:,:,:,:)
        !> Overlap matrix 
        real(wp), allocatable :: s_mtmt(:,:,:,:,:,:,:,:) 
        !> Temporary polarizability arrays
        ! complex(wp), allocatable :: tmp_p1(:,:,:), tmp_p2(:,:,:,:,:)  
        ! complex(wp), allocatable :: tmp_p1(:,:,:,:,:,:,:), tmp_p2(:,:,:,:,:,:,:)
        complex(wp), allocatable :: tmp_p1(:,:,:,:,:,:,:,:), tmp_p2(:,:,:,:,:,:,:,:)
        !> Combined NLM indices 
        integer, allocatable :: nlm_indices(:, :, :)
        
        ! Test and delete
        real :: t_i, t_f
        complex(wp) :: tmp1, tmp2
        
        nbigR = kqset%nkpt
        !!print*, '*******', kqset%nkpt, nbigR
        ! lm_lo_max = max(lmmaxapw, lolmmax)
        lm_lo_max = max(lmmaxapw, nlomax)
        
        call nlm_index_contracting(nlm_tot, nlm_indices)

        ! allocate(tmp_p1(apwordmax, lm_lo_max, natmtot))
        ! allocate(tmp_p2(apwordmax, lm_lo_max, apwordmax, lm_lo_max, natmtot))
        ! allocate(tmp_p1(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, nlm_tot, natmtot))
        allocate(tmp_p1(apwordmax, lm_lo_max, natmtot, apwordmax, lm_lo_max, nlm_tot, natmtot, 4))
        ! allocate(tmp_p2(apwordmax, lm_lo_max, apwordmax, lm_lo_max, natmtot, nlm_tot, natmtot))
        ! allocate(tmp_p2(apwordmax, lm_lo_max, apwordmax, lm_lo_max, natmtot, nlm_tot, natmtot))
        allocate(tmp_p2(apwordmax, lm_lo_max, apwordmax, lm_lo_max, natmtot, nlm_tot, natmtot, 4))
        allocate(pola_R_mtmt(nlm_tot, natmtot, nlm_tot, natmtot, nbigR))

        !> Initialization of "pola_R_mtmt" array
        pola_R_mtmt = cmplx(0.0_wp, 0.0_wp, kind=wp)

        print*, 'lmaxapw, lmmaxapw, and nlm dimension ', input%groundstate%lmaxapw, lmmaxapw, nlm_tot
        ! print*, 'gaunt_coefficient tolerance', tolerance
        
        call overlap_mtmt(nlm_indices, s_mtmt)
        
        call green_R_mtmt(tau, green_R_occ, green_R_uno) 


        !! Separate-separate multiplication --- modifying according to Andris's suggestion ----------
        do ir = 1, nbigR
          call CPU_TIME(t_i)

          tmp_p1 = cmplx(0.0_wp, 0.0_wp, kind=wp)
          io2 = 1 
          do ias1 = 1, natmtot
            do nlm1 = 1, nlm_tot
              is1  = nlm_indices(nlm1, ias1, 1)
              mbn1 = nlm_indices(nlm1, ias1, 3)
              mbl1 = nlm_indices(nlm1, ias1, 4)
              mbm1 = nlm_indices(nlm1, ias1, 5)
              do l1 = 0, input%groundstate%lmaxapw
                do m1 = - l1, l1
                  lm1 = idxlm(l1, m1)
                  m3 = - mbm1 + m1 
                  l3min = abs(m3)
                  l3max = min(mbl1+l1, input%groundstate%lmaxapw)
                  do l3 = l3min, l3max 
                    gaunt_coefficient = getgauntcoef(l3, mbl1, l1, m3, mbm1)
                    if(abs(gaunt_coefficient)>tolerance) then
                      lm3 = idxlm(l3, m3)

                      do io1 = 1, apwordmax !apword(l1,is)
                        do io3 = 1, apwordmax !apword(l3,is)
                          !! p1_{APW-APW} (1st term)
                          tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 1) = & 
                          tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 1) + &
                          conjg(green_R_occ(:, :, :, io3, lm3, ias1, ir, 1)) * &
                          bradketa(2, mbn1, l1, io1, l3, io3, ias1) * gaunt_coefficient

                          !! p1_{LO-APW} (1st term)
                          ! io2 = 1  ! move outside the loop
                          tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 3) = & 
                          tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 3) + & 
                          conjg(green_R_occ(io2, :, :, io3, lm3, ias1, ir, 3)) * & 
                          bradketa(2, mbn1, l1, io1, l3, io3, ias1) * gaunt_coefficient
                        enddo 
                        io3 = 1
                        do ilo3 = 1, nlorb(is1)
                          if(lorbl(ilo3, is1)==l3) then
                            !! p1_{APW-APW} (2nd term)
                            tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 1) = & 
                            tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 1) + &
                            conjg(green_R_occ(:, :, :, io3, lm3, ias1, ir, 2)) * &  !! lm3 --> ilo3
                            bradketa(3, mbn1, l1, io1, ilo3, io3, ias1) * gaunt_coefficient

                            !! p1_{LO-APW} (2nd term)
                            ! io2 = 1  ! move outside the loop
                            tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 3) = & 
                            tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 3) + & 
                            conjg(green_R_occ(io2, :, :, io3, lm3, ias1, ir, 4)) * &  !! lm3 --> ilo3
                            bradketa(3, mbn1, l1, io1, ilo3, io3, ias1) * gaunt_coefficient
                          endif 
                        enddo 
                      enddo 

                      io1 = 1
                      do ilo1 = 1, nlorb(is1)
                        if(lorbl(ilo1, is1)==l1) then
                          do io3 = 1, apwordmax !apword(l3,is)
                            !! p1_{APW-LO} (1st term)
                            tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 2) = &   !! lm1 --> ilo1
                            tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 2) + &   !! lm1 --> ilo1
                            conjg(green_R_occ(:, :, :, io3, lm3, ias1, ir, 1)) * &
                            bradketlo(2, mbn1, ilo1, l3, io3, ias1) * gaunt_coefficient

                            !! p1_{LO-LO} (1st term) 
                            ! io2 = 1  ! move outside the loop
                            tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 4) = &    !! lm1 --> ilo1
                            tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 4) + &    !! lm1 --> ilo1
                            conjg(green_R_occ(io2, :, :, io3, lm3, ias1, ir, 3)) * &
                            bradketlo(2, mbn1, ilo1, l3, io3, ias1) * gaunt_coefficient
                          enddo 
                          io3 = 1
                          do ilo3 = 1, nlorb(is1)
                            if(lorbl(ilo3, is1)==l3) then
                              !! p1_{APW-LO} (2nd term)
                              tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 2) = &   !! lm1 --> ilo1
                              tmp_p1(:, :, :, io1, lm1, nlm1, ias1, 2) + &   !! lm1 --> ilo1
                              conjg(green_R_occ(:, :, :, io3, lm3, ias1, ir, 2)) * &  !! lm3 --> ilo3
                              bradketlo(3, mbn1, ilo1, ilo3, io3, ias1) * gaunt_coefficient

                              !! p1_{LO-LO} (2nd term) 
                              ! io2 = 1  ! move outside the loop
                              tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 4) = &    !! lm1 --> ilo1
                              tmp_p1(io2, :, :, io1, lm1, nlm1, ias1, 4) + &    !! lm1 --> ilo1
                              conjg(green_R_occ(io2, :, :, io3, lm3, ias1, ir, 4)) * &  !! lm3 --> ilo3
                              bradketlo(3, mbn1, ilo1, ilo3, io3, ias1) * gaunt_coefficient
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

          call CPU_TIME(t_f)
          ! print*, 'min, max value of real tmp_p1', minval(real(tmp_p1)), maxval(real(tmp_p1))
          ! print*, 'min, max value of imag tmp_p1', minval(imag(tmp_p1)), maxval(imag(tmp_p1))
          ! print*, 'sum tmp_p1', sum(tmp_p1), 'time', t_f - t_i
          write(9010,*) sum(tmp_p1)
          call CPU_TIME(t_i)

          tmp_p2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
          io = 1 
          do ias1 = 1, natmtot
            do nlm1 = 1, nlm_tot
              is1  = nlm_indices(nlm1, ias1, 1)
              do l1 = 0, input%groundstate%lmaxapw
                do m1 = - l1, l1
                  lm1 = idxlm(l1, m1)
                  do is = 1, nspecies
                    do ia = 1, natoms(is)
                      ias = idxas(ia,is)
                      do l2 = 0, input%groundstate%lmaxapw
                        do m2 = - l2, l2 
                          lm2 = idxlm(l2, m2)

                          do io2 = 1, apwordmax !apword(l2,is)
                            do io1 = 1, apwordmax !apword(l1,is)
                              !! p2_{APW-APW} (1st term)
                              tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 1) = &
                              tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 1) + &
                              green_R_uno(:, :, ias, io1, lm1, ias1, ir, 1) * &
                              tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 1)

                              !! p2_{LO-APW} (1st term)
                              tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 3) = &
                              tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 3) + &
                              green_R_uno(io, :, ias, io1, lm1, ias1, ir, 3) * &
                              tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 1)
                            enddo 
                            io1 = 1
                            do ilo1 = 1, nlorb(is1)
                              if(lorbl(ilo1, is1)==l1) then
                                !! p2_{APW-APW} (2nd term)
                                tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 1) = &
                                tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 1) + &
                                green_R_uno(:, :, ias, io1, lm1, ias1, ir, 2) * &  !! lm1 --> ilo1
                                tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 2)     !! lm1 --> ilo1

                                !! p2_{LO-APW} (2nd term)
                                tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 3) = &
                                tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 3) + &
                                green_R_uno(io, :, ias, io1, lm1, ias1, ir, 4) * &  !! lm1 --> ilo1
                                tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 2)  !! lm1 --> ilo1
                              endif 
                            enddo 
                          enddo

                          io2 = 1
                          do ilo2 = 1, nlorb(is)
                            if(lorbl(ilo2, is)==l2) then
                              do io1 = 1, apwordmax
                                !! p2_{APW-LO} (1st term)
                                tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 2) = &  !! lm2 --> ilo2
                                tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 2) + &  !! lm2 --> ilo2
                                green_R_uno(:, :, ias, io1, lm1, ias1, ir, 1) * &
                                tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 3)  !! lm2 --> ilo2

                                !! p2_{LO-LO} (1st term)
                                tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 4) = &  !! lm2 --> ilo2
                                tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 4) + &  !! lm2 --> ilo2
                                green_R_uno(io, :, ias, io1, lm1, ias1, ir, 3) * &
                                tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 3)  !! lm2 --> ilo2
                              enddo 
                              io1 = 1
                              do ilo1 = 1, nlorb(is1)
                                if(lorbl(ilo1, is1)==l1) then
                                  !! p2_{APW-LO} (2nd term)
                                  tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 2) = &  !! lm2 --> ilo2
                                  tmp_p2(:, :, io2, lm2, ias, nlm1, ias1, 2) + &  !! lm2 --> ilo2
                                  green_R_uno(:, :, ias, io1, lm1, ias1, ir, 2) * &  !! lm1 --> ilo1
                                  tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 4)  !! lm2 --> ilo2 & lm1 --> ilo1

                                  !! p2_{LO-LO} (2nd term)
                                  tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 4) = &  !! lm2 --> ilo2
                                  tmp_p2(io, :, io2, lm2, ias, nlm1, ias1, 4) + &  !! lm2 --> ilo2
                                  green_R_uno(io, :, ias, io1, lm1, ias1, ir, 4) * &  !! lm1 --> ilo1
                                  tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1, 4)  !! lm2 --> ilo2 & lm1 --> ilo1
                                endif
                              enddo 
                            endif 
                          enddo 

                        enddo 
                      enddo 
                    enddo 
                  enddo
                enddo
              enddo
            enddo
          enddo
          call CPU_TIME(t_f)
          ! print*, 'min, max value of real tmp_p2', minval(real(tmp_p2)), maxval(real(tmp_p2))
          ! print*, 'min, max value of imag tmp_p2', minval(imag(tmp_p2)), maxval(imag(tmp_p2))
          ! print*, 'sum tmp_p2', sum(tmp_p2), 'time', t_f - t_i
          ! write(3031,*) sum(real(tmp_p2)), sum(imag(tmp_p2)), ir
          write(100020,*) sum(tmp_p2)

          call CPU_TIME(t_i)
          
          do ias1 = 1, natmtot
            do nlm1 = 1, nlm_tot
              do ias = 1, natmtot
                do nlm = 1, nlm_tot
                  is  = nlm_indices(nlm, ias, 1)
                  mbn = nlm_indices(nlm, ias, 3)
                  mbl = nlm_indices(nlm, ias, 4)
                  mbm = nlm_indices(nlm, ias, 5)
                  do l = 0, input%groundstate%lmaxapw
                    do m = - l, l
                      lm = idxlm(l, m)
                      m2 = - mbm + m
                      l2min = abs(m2)
                      l2max = min(mbl+l, input%groundstate%lmaxapw)
                      do l2 = l2min, l2max 
                        gaunt_coefficient = getgauntcoef(l2, mbl, l, m2, mbm)
                        if(abs(gaunt_coefficient)>tolerance) then
                          lm2 = idxlm(l2, m2)

                          do io = 1, apwordmax !apword(l,is)
                            do io2 = 1, apwordmax !apword(l2,is)
                              !! {APW-APW} (1st term)
                              pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
                              pola_R_mtmt(nlm, ias, nlm1, ias1, ir) - &
                              bradketa(2, mbn, l, io, l2, io2, ias) * gaunt_coefficient * &
                              tmp_p2(io, lm, io2, lm2, ias, nlm1, ias1, 1)
                            enddo 
                            io2 = 1
                            do ilo2 = 1, nlorb(is)
                              if(lorbl(ilo2, is)==l2) then
                                !! {APW-LO} (2nd term)
                                pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
                                pola_R_mtmt(nlm, ias, nlm1, ias1, ir) - &
                                bradketa(3, mbn, l, io, ilo2, io2, ias) * gaunt_coefficient * & 
                                tmp_p2(io, lm, io2, lm2, ias, nlm1, ias1, 2)  !! lm2 --> ilo2
                              endif 
                            enddo 
                          enddo

                          io = 1
                          do ilo = 1, nlorb(is)
                            if(lorbl(ilo, is)==l) then
                              do io2 = 1, apwordmax !apword(l2,is)
                                !! {LO-APW} (3rd term)
                                pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
                                pola_R_mtmt(nlm, ias, nlm1, ias1, ir) - &
                                bradketlo(2, mbn, ilo, l2, io2, ias) * gaunt_coefficient * & 
                                tmp_p2(io, lm, io2, lm2, ias, nlm1, ias1, 3)  !! lm --> ilo
                              enddo 
                              io2 = 1
                              do ilo2 = 1, nlorb(is)
                                if(lorbl(ilo2, is)==l2) then
                                  !! {LO-LO} (4th term)
                                  pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
                                  pola_R_mtmt(nlm, ias, nlm1, ias1, ir) - &
                                  bradketlo(3, mbn, ilo, ilo2, io2, ias) * gaunt_coefficient * & 
                                  tmp_p2(io, lm, io2, lm2, ias, nlm1, ias1, 4)  !! lm --> ilo & lm2 --> ilo2
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
            enddo 
          enddo 

          call CPU_TIME(t_f)
          print*, 'time 1', t_f - t_i
          call CPU_TIME(t_i)
          write(400,*) sum(real(pola_R_mtmt(:,:,:,:,ir))), sum(imag(pola_R_mtmt(:,:,:,:,ir))), ir
          write(4010,*) sum(real(pola_R_mtmt)), sum(imag(pola_R_mtmt)), ir
          ! print*, 'min, max value of real pR', minval(real(pola_R_mtmt)), maxval(real(pola_R_mtmt))
          ! print*, 'min, max value of imag pR', minval(imag(pola_R_mtmt)), maxval(imag(pola_R_mtmt))
          print*, 'sum pola_R mt-mt', sum(pola_R_mtmt)
          print*, 'test pola_R_mtmt', minval(real(pola_R_mtmt)), maxval(imag(pola_R_mtmt))
          call CPU_TIME(t_f)
          print*, 'time 2', t_f - t_i

        enddo   ! ir loop ends

        !! Spin 
        pola_R_mtmt = 2.0_wp * pola_R_mtmt   !! spin included here but (TODO:) do it properly 

        ! !! Polarizability in R and tau -------------------------------------
        ! i = nint(tau*2.0_wp)
        ! do ir = 1, nbigR
        !   do ias1 = 1, natmtot
        !     do nlm1 = 1, nlm_tot
        !       do ias = 1, natmtot
        !         do nlm = 1, nlm_tot
        !           write(1900+i,*) pola_R_mtmt(nlm, ias, nlm1, ias1, ir)
        !         enddo 
        !       enddo 
        !     enddo 
        !   enddo 
        ! enddo 
        !!----------------------------------------------------------------------


        !! Semi-separate multiplication
        ! do ir = 1, nbigR
          
        !   do i = 2, 3   !! APW/LO 
        !     do j = 2, 3   !! APW/LO 
          
        !       call CPU_TIME(t_i)
        !       do ias1 = 1, natmtot
        !         do nlm1 = 1, nlm_tot

        !           tmp_p2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        !           do lm1 = 1, lm_lo_max
        !             do io1 = 1, apwordmax !apword(l1,is)

        !               tmp_p1 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        !               do lm2 = 1, lm_lo_max
        !                 do io2 = 1, apwordmax !apword(l2,is)
        !                   do ias = 1, natmtot
        !                     do lm = 1, lm_lo_max
        !                       do io = 1, apwordmax !apword(l,is)
        !                         tmp_p1(io, lm, ias) = tmp_p1(io, lm, ias) + &
        !                         ! conjg(green_R_occ(io, lm, ias, io2, lm2, ias1, ir, j, i)) * &
        !                         !! Do we need conjugate of "green_R_occ" ???
        !                         green_R_occ(io, lm, ias, io2, lm2, ias1, ir, j, i) * &
        !                         s_mtmt(io2, lm2, io1, lm1, nlm1, ias1, j, i)
        !                         !! check if this needs complex conjugate of "s_mtmt" matrix
        !                       enddo 
        !                     enddo 
        !                   enddo 
        !                 enddo 
        !               enddo 

        !               do ias = 1, natmtot
        !                 do lm = 1, lm_lo_max
        !                   do io = 1, apwordmax !apword(l,is)
        !                     do lm2 = 1, lm_lo_max
        !                       do io2 = 1, apwordmax !apword(l2,is)
        !                         tmp_p2(io2, lm2, io, lm, ias) = tmp_p2(io2, lm2, io, lm, ias) + &
        !                         green_R_uno(io, lm, ias, io1, lm1, ias1, ir, j, i) * tmp_p1(io2, lm2, ias)
        !                       enddo 
        !                     enddo 
        !                   enddo 
        !                 enddo
        !               enddo

        !             enddo   !! io1
        !           enddo   !! lm1

        !           do ias = 1, natmtot
        !             do nlm = 1, nlm_tot
        !               do lm = 1, lm_lo_max
        !                 do io = 1, apwordmax !apword(l,is)
        !                   do lm2 = 1, lm_lo_max
        !                     do io2 = 1, apwordmax !apword(l2,is)
        !                       pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
        !                       pola_R_mtmt(nlm, ias, nlm1, ias1, ir) - &  
        !                       ! pola_R_mtmt(nlm, ias, nlm1, ias1, ir) + &  !! it shoud b e "-" sign 
        !                       s_mtmt(io2, lm2, io, lm, nlm, ias, j, i) * &
        !                       tmp_p2(io2, lm2, io, lm, ias)
        !                     enddo 
        !                   enddo 
        !                 enddo
        !               enddo
        !             enddo
        !           enddo

        !         enddo   !! nlm1
        !       enddo   !! ias1
        !       call CPU_TIME(t_f)

        !       print*, 'time 1', t_f - t_i
        !       call CPU_TIME(t_i)
        !       write(400,*) j, i, sum(real(pola_R_mtmt(:,:,:,:,ir))), sum(imag(pola_R_mtmt(:,:,:,:,ir))), ir
        !       write(4010,*) sum(real(pola_R_mtmt)), sum(imag(pola_R_mtmt)), ir
        !       ! print*, 'min, max value of real pR', minval(real(pola_R_mtmt)), maxval(real(pola_R_mtmt))
        !       ! print*, 'min, max value of imag pR', minval(imag(pola_R_mtmt)), maxval(imag(pola_R_mtmt))
        !       print*, 'sum pola_R mt-mt', sum(pola_R_mtmt)
        !       call CPU_TIME(t_f)
        !       print*, 'time 2', t_f - t_i

        !     enddo   !! j loop ends
        !   enddo   !! i loop ends
          
        ! enddo   ! ir loop ends
        !-----------------------------------------------------------------------------------------

        ! !! Separate-separate multiplication
        ! do ir = 1, nbigR
          
        !   do i = 3, 3   !! APW/LO 
        !     do j = 3, 3   !! APW/LO 
          
        !       call CPU_TIME(t_i)
        !       tmp_p1 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        !       do ias1 = 1, natmtot
        !         do nlm1 = 1, nlm_tot
        !           do lm1 = 1, lm_lo_max
        !             do io1 = 1, apwordmax !apword(l1,is)
        !               do lm2 = 1, lm_lo_max
        !                 do io2 = 1, apwordmax !apword(l2,is)
        !                   do ias = 1, natmtot
        !                     do lm = 1, lm_lo_max
        !                       do io = 1, apwordmax !apword(l,is)
        !                         ! tmp_p1(lm, ias, lm1, nlm, ias1) = tmp_p1(lm, ias, lm1, nlm, ias1) + &
        !                         tmp_p1(io, lm, ias, io1, lm1, nlm1, ias1) = & 
        !                         tmp_p1(io, lm, ias, io1, lm1, nlm1, ias1) + &
        !                         conjg(green_R_occ(io, lm, ias, io2, lm2, ias1, ir, j, i)) * &
        !                         s_mtmt(io2, lm2, io1, lm1, nlm1, ias1, j, i)
        !                         !! check if this needs complex conjugate of "s_mtmt" matrix
        !                       enddo   !! io
        !                     enddo   !! lm
        !                   enddo   !! ias
        !                 enddo   !! io2
        !               enddo   !! lm2
        !             enddo   !! io1
        !           enddo   !! lm1
        !         enddo   !! nlm1
        !       enddo   !! ias1
        !       call CPU_TIME(t_f)
        !       ! print*, 'min, max value of real tmp_p1', minval(real(tmp_p1)), maxval(real(tmp_p1))
        !       ! print*, 'min, max value of imag tmp_p1', minval(imag(tmp_p1)), maxval(imag(tmp_p1))
        !       print*, 'sum tmp_p1', sum(tmp_p1), 'time', t_f - t_i
        !       write(100010,*) sum(tmp_p1)
        !       call CPU_TIME(t_i)

        !       tmp_p2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        !       do ias1 = 1, natmtot
        !         do nlm1 = 1, nlm_tot
        !           do lm1 = 1, lm_lo_max
        !             do io1 = 1, apwordmax !apword(l1,is)
        !               do ias = 1, natmtot
        !                 do lm = 1, lm_lo_max
        !                   do io = 1, apwordmax !apword(l,is)
        !                     do lm2 = 1, lm_lo_max
        !                       do io2 = 1, apwordmax !apword(l2,is)
        !                         ! tmp_p2(lm2, lm, ias, nlm1, ias1) = &
        !                         tmp_p2(io2, lm2, io, lm, ias, nlm1, ias1) = &
        !                         tmp_p2(io2, lm2, io, lm, ias, nlm1, ias1) + &
        !                         green_R_uno(io, lm, ias, io1, lm1, ias1, ir, j, i) * &
        !                         tmp_p1(io2, lm2, ias, io1, lm1, nlm1, ias1)
        !                       enddo 
        !                     enddo 
        !                   enddo 
        !                 enddo
        !               enddo
        !             enddo
        !           enddo
        !         enddo
        !       enddo
        !       call CPU_TIME(t_f)
        !       ! print*, 'min, max value of real tmp_p2', minval(real(tmp_p2)), maxval(real(tmp_p2))
        !       ! print*, 'min, max value of imag tmp_p2', minval(imag(tmp_p2)), maxval(imag(tmp_p2))
        !       print*, 'sum tmp_p2', sum(tmp_p2), 'time', t_f - t_i
        !       ! write(3031,*) sum(real(tmp_p2)), sum(imag(tmp_p2)), ir
        !       write(100020,*) sum(tmp_p2)

        !       call CPU_TIME(t_i)
        !       do ias1 = 1, natmtot
        !         do nlm1 = 1, nlm_tot
        !           do ias = 1, natmtot
        !             do nlm = 1, nlm_tot
        !               do lm = 1, lm_lo_max
        !                 do io = 1, apwordmax !apword(l,is)
        !                   do lm2 = 1, lm_lo_max
        !                     do io2 = 1, apwordmax !apword(l2,is)
        !                       pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
        !                       pola_R_mtmt(nlm, ias, nlm1, ias1, ir) - &
        !                       s_mtmt(io2, lm2, io, lm, nlm, ias, j, i) * &
        !                       tmp_p2(io2, lm2, io, lm, ias, nlm1, ias1)
        !                     enddo 
        !                   enddo 
        !                 enddo
        !               enddo
        !             enddo
        !           enddo
        !         enddo
        !       enddo

        !       call CPU_TIME(t_f)
        !       print*, 'time 1', t_f - t_i
        !       call CPU_TIME(t_i)
        !       write(400,*) sum(real(pola_R_mtmt(:,:,:,:,ir))), sum(imag(pola_R_mtmt(:,:,:,:,ir))), ir
        !       write(4010,*) sum(real(pola_R_mtmt)), sum(imag(pola_R_mtmt)), ir
        !       ! print*, 'min, max value of real pR', minval(real(pola_R_mtmt)), maxval(real(pola_R_mtmt))
        !       ! print*, 'min, max value of imag pR', minval(imag(pola_R_mtmt)), maxval(imag(pola_R_mtmt))
        !       print*, 'sum pola_R mt-mt', sum(pola_R_mtmt)
        !       call CPU_TIME(t_f)
        !       print*, 'time 2', t_f - t_i

        !     enddo   !! j loop ends
        !   enddo   !! i loop ends
          
        ! enddo   ! ir loop ends
        ! !!----------------------------------------------------------------------

        ! pola_R_mtmt = 2.0_wp * pola_R_mtmt   !! spin included here but (TODO:) do it properly 

        ! !! Polarizability in R and tau -------------------------------------
        ! i = nint(tau*2.0_wp)
        ! do ir = 1, nbigR
        !   do ias1 = 1, natmtot
        !     do nlm1 = 1, nlm_tot
        !       do ias = 1, natmtot
        !         do nlm = 1, nlm_tot
        !           write(1900+i,*) pola_R_mtmt(nlm, ias, nlm1, ias1, ir)
        !         enddo 
        !       enddo 
        !     enddo 
        !   enddo 
        ! enddo 
        !! Polarizability in R and tau -------------------------------------


        ! This part is not working
        ! ! ! do ir = 1, nbigR
          
        ! ! !   do i = 3, 3   !! APW/LO 
        ! ! !     do j = 3, 3   !! APW/LO 
          
        ! ! !       call CPU_TIME(t_i)
        ! ! !       do ias1 = 1, natmtot
        ! ! !         do nlm1 = 1, nlm_tot
        ! ! !           do ias = 1, natmtot
        ! ! !             do nlm = 1, nlm_tot

        ! ! !               do lm = 1, lm_lo_max
        ! ! !                 do io = 1, apwordmax !apword(l,is)
        ! ! !                   do lm2 = 1, lm_lo_max
        ! ! !                     do io2 = 1, apwordmax !apword(l2,is)

        ! ! !                       tmp2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        ! ! !                       do lm1 = 1, lm_lo_max
        ! ! !                         do io1 = 1, apwordmax !apword(l1,is)

        ! ! !                           tmp1 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        ! ! !                           do lm3 = 1, lm_lo_max
        ! ! !                             do io3 = 1, apwordmax !apword(l3,is)
        ! ! !                              tmp1 = tmp1 + conjg(green_R_occ(io2, lm2, ias, io3, lm3, ias1, ir, j, i)) * &
        ! ! !                              s_mtmt(io1, lm1, io3, lm3, nlm1, ias1, j, i)
        ! ! !                              !! check if this needs complex conjugate of "s_mtmt" matrix
        ! ! !                             enddo 
        ! ! !                           enddo 

        ! ! !                           tmp2 = tmp2 + green_R_uno(io, lm, ias, io1, lm1, ias1, ir, j, i) * tmp1

        ! ! !                         enddo 
        ! ! !                       enddo  

        ! ! !                       pola_R_mtmt(nlm, ias, nlm1, ias1, ir) = &
        ! ! !                       pola_R_mtmt(nlm, ias, nlm1, ias1, ir) + &
        ! ! !                       s_mtmt(io2, lm2, io, lm, nlm, ias, j, i) * tmp2

        ! ! !                     enddo   !! io2
        ! ! !                   enddo   !! lm2
        ! ! !                 enddo   !! io
        ! ! !               enddo   !! lm

        ! ! !             enddo   !! nlm
        ! ! !           enddo   !! ias
        ! ! !         enddo   !! nlm1
        ! ! !       enddo   !! ias1
        ! ! !       call CPU_TIME(t_f)

        ! ! !       print*, 'time 1', t_f - t_i
        ! ! !       call CPU_TIME(t_i)
        ! ! !       write(400,*) sum(real(pola_R_mtmt(:,:,:,:,ir))), sum(imag(pola_R_mtmt(:,:,:,:,ir))), ir
        ! ! !       write(4010,*) sum(real(pola_R_mtmt)), sum(imag(pola_R_mtmt)), ir
        ! ! !       ! print*, 'min, max value of real pR', minval(real(pola_R_mtmt)), maxval(real(pola_R_mtmt))
        ! ! !       ! print*, 'min, max value of imag pR', minval(imag(pola_R_mtmt)), maxval(imag(pola_R_mtmt))
        ! ! !       print*, 'sum pola_R mt-mt', sum(pola_R_mtmt)
        ! ! !       call CPU_TIME(t_f)
        ! ! !       print*, 'time 2', t_f - t_i

        ! ! !     enddo   !! j loop ends
        ! ! !   enddo   !! i loop ends
          
        ! ! ! enddo   ! ir loop ends

        do ir = 1, nbigR
         write(402,*) sum(real(pola_R_mtmt(:,:,:,:,ir))), sum(imag(pola_R_mtmt(:,:,:,:,ir))), ir
        enddo
        write(41,*) sum(real(pola_R_mtmt)), sum(imag(pola_R_mtmt))
        !
        deallocate(green_R_occ)
        deallocate(green_R_uno)
        deallocate(s_mtmt)
        deallocate(nlm_indices)
        deallocate(tmp_p1)
        deallocate(tmp_p2)
        !
    end subroutine polarizability_R_mtmt

end module mod_polarizability_R_mtmt
