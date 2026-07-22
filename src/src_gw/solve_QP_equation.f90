!> Solve the quasiparticle equation for the current GW setup.
subroutine solve_QP_equation()
    use modinput
    use modmain,        only: efermi
    use modgw,          only: ibgw, nbgw, nvelgw, nbandsgw, evalqp, eferqp
    use mod_vxc,        only: vxcnn
    use mod_selfenergy, only: selfex, selfec, sigc, znorm, freq_selfc, deltaE
    use mod_bands,      only: nomax, ikvbm, evalfv
    use modmpi,         only: terminate
    use mod_pade
    use precision, only: i32, dp
    use self_consistent_eigenvalue_gw0, only: use_evgw0_input_qp, get_evalqp_evgw0_pointer
    use to_char_conversion, only: to_char
    implicit none
    integer(i32), parameter :: nitermax = 1000_i32
    real(dp), parameter :: etol = 1.0e-4_dp
    integer(i32) :: iter, ik, ib, nz, n_kpoints
    real(dp) :: enk, eqp, eqp_prev, diff, dzf2
    complex(dp) :: sx, sc, de
    complex(dp) :: dsigma, znk
    logical    :: converged
    real(dp) :: eqp_evgw0
    real(dp), pointer :: evalqp_evgw0(:, :)

    if ( use_evgw0_input_qp() ) then
       call get_evalqp_evgw0_pointer(evalqp_evgw0)
    end if

    !-----------------------------------------
    ! Alignment of the chemical potential:
    !   ef + de = ef + Sigma(kf, ef + de)
    !------------------------------------------
    n_kpoints = size( selfec, 3 )
    select case (input%gw%selfenergy%eshift)
        case(0)
            ! no shift
            de = 0.0_dp
        case(1)
            ! following Lucia Reining's book
            enk = evalfv(nomax,ikvbm)-efermi
            sx = enk + selfex(nomax,ikvbm) - vxcnn%diag_elements(nomax,ikvbm)
            eqp = enk
            eqp_prev = eqp
            converged = .false.
            do iter = 1, nitermax
                call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(nomax,:,ikvbm), &
                                eqp, sc, dsigma )
                eqp = sx + sc
                diff = eqp - eqp_prev
                if ( abs(diff) < etol ) then
                    converged = .true.
                    exit
                else
                    ! Next iteration
                    eqp_prev = eqp
                end if
            end do
            if (.not.converged) write(*,*) 'Problem with convergence!'
            de = eqp - enk
        case(2)
            ! following Bruneval&Gatti's article
            enk = evalfv(nomax,ikvbm)-efermi
            call get_selfc(freq_selfc%nomeg, freq_selfc%freqs, selfec(nomax,:,ikvbm), &
                        enk, sc, dsigma)
            de = selfex(nomax,ikvbm) + sc - vxcnn%diag_elements(nomax,ikvbm)
        case default
            call terminate('Non supported values of eshift=' // trim(to_char(input%gw%selfenergy%eshift)))
    end select
    deltaE = de%re
    ! print*, 'QP energy shift delta_e = ', deltaE

    !--------------------------------------------
    ! Solve QP equation
    !--------------------------------------------
    do ik = 1, n_kpoints
        do ib = ibgw, nbgw

            enk = evalfv(ib,ik)-efermi
            eqp = enk
            eqp_prev = eqp

            converged = .false.
            do iter = 1, nitermax

                if (use_evgw0_input_qp()) then
                    eqp_evgw0 = evalqp_evgw0(ib,ik)
                    call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                         eqp_evgw0, sigc(ib,ik), dsigma )
                    znk = zone / (zone-dsigma)
                    znorm(ib,ik) = znk%re
                    eqp = eqp_evgw0 + znorm(ib,ik) * (selfex(ib,ik)%re + sigc(ib,ik)%re - &
                         vxcnn%diag_elements(ib,ik)%re + enk - eqp_evgw0)
                    converged = .true.
                    exit
                end if

                select case (input%gw%selfenergy%eqpsolver)
                    case(0)
                        ! Perturbative solution (single iteration)
                        call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                                        enk, sigc(ib,ik), dsigma )
                        znk = zone / (zone-dsigma)
                        znorm(ib,ik) = znk%re
                        eqp = enk + znorm(ib,ik) * (selfex(ib,ik)%re + sigc(ib,ik)%re - &
                              vxcnn%diag_elements(ib,ik)%re) + (1.0_dp - znorm(ib,ik)) * deltaE
                        converged = .true.
                        exit
                    case(1)
                        ! Perturbative solution without renormalization
                        call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                                        enk, sigc(ib,ik), dsigma )
                        eqp = enk + selfex(ib,ik)%re + sigc(ib,ik)%re - vxcnn%diag_elements(ib,ik)%re
                        znorm(ib,ik) = 1.0_dp
                        converged = .true.
                        exit
                    case(2)
                        ! Iterative solution
                        call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                                        eqp-deltaE, sigc(ib,ik), dsigma )
                        eqp = enk + selfex(ib,ik)%re + sigc(ib,ik)%re - vxcnn%diag_elements(ib,ik)%re
                        znorm(ib,ik) = 1.0_dp
                    case default
                        call terminate('Error(solve_QP_equation) Non supported value: eqpsolver =' // &
                          trim(to_char(input%gw%selfenergy%eqpsolver)))
                end select

                ! Error function
                diff = eqp - eqp_prev

                if ( abs(diff) < etol ) then
                    converged = .true.
                    print*, '# iterations ', iter
                    exit
                else
                    ! Next iteration
                    eqp_prev = eqp
                end if

            end do ! iter

            if (.not.converged) then
                write(*,*)
                write(*,'(a,i0,a)') 'Warning(solve_QP_equation) Solution of the quasiparticle equation is not converged after ', nitermax, ' iterations!'
                write(*,'(a,f8.4,a,f8.4)') 'Absolute error = ', abs(diff), ' > ', etol
                write(*,'(a,2i4)') 'ib, ik = ', ib, ik
                write(*,'(a,2f8.4)') 'diff = ', diff
            end if

            ! Quasiparticle energy
            evalqp(ib,ik) = eqp

        end do ! ib
    end do ! ik

end subroutine
