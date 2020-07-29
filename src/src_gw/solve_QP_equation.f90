
subroutine solve_QP_equation()
    use modinput
    use modmain,        only: efermi
    use modgw,          only: kset, ibgw, nbgw, nvelgw, nbandsgw, evalqp, eferqp, evalfv
    use mod_vxc,        only: vxcnn
    use mod_selfenergy, only: selfex, selfec, sigc, znorm, freq_selfc, deltaE
    use mod_bands,      only: nomax, ikvbm
    use mod_pade
    use m_getunit
    implicit none
    integer(4), parameter :: nitermax = 1000
    real(8),    parameter :: etol = 1.d-4
    integer(4) :: iter, ik, ib, nz
    real(8)    :: enk, eqp, eqp_prev, diff, dzf2
    complex(8) :: sx, sc, de
    complex(8) :: dsigma, znk
    logical    :: converged

    !-----------------------------------------
    ! Alignment of the chemical potential:
    !   ef + de = ef + Sigma(kf, ef + de)
    !------------------------------------------
    select case (input%gw%selfenergy%eshift)
        case(0)
            ! no shift
            de = 0.d0
        case(1)
            ! following Lucia Reining's book
            enk = evalfv(nomax,ikvbm)-efermi
            sx = enk + selfex(nomax,ikvbm) - vxcnn(nomax,ikvbm)
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
            de = selfex(nomax,ikvbm) + sc - vxcnn(nomax,ikvbm)
        case default
            write(*,*) 'Non supported values of eshift=', input%gw%selfenergy%eshift
            stop
    end select
    deltaE = dble(de)
    ! print*, 'QP energy shift delta_e = ', deltaE

    !--------------------------------------------
    ! Solve QP equation
    !--------------------------------------------
    do ik = 1, kset%nkpt
        do ib = ibgw, nbgw

            enk = evalfv(ib,ik)-efermi
            eqp = enk
            eqp_prev = eqp

            converged = .false.
            do iter = 1, nitermax

                select case (input%gw%selfenergy%eqpsolver)
                    case(0)
                        ! Perturbative solution (single iteration)
                        call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                                        enk, sigc(ib,ik), dsigma )
                        znk = zone / (zone-dsigma)
                        znorm(ib,ik)  = dble(znk)
                        eqp = enk + znorm(ib,ik)*dble(selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik)) + &
                              (1.d0-znorm(ib,ik))*deltaE
                        converged = .true.
                        exit
                    case(1)
                        ! Perturbative solution without renormalization
                        call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                                        enk, sigc(ib,ik), dsigma )
                        eqp = enk + dble(selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik))
                        znorm(ib,ik)  = 1.d0
                        converged = .true.
                        exit
                    case(2)
                        ! Iterative solution
                        call get_selfc( freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                                        eqp-deltaE, sigc(ib,ik), dsigma )
                        eqp = enk + selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik)
                        znorm(ib,ik)  = 1.d0
                    case default
                        write(*,*) 'Error(solve_QP_equation) Non supported value: eqpsolver =', input%gw%selfenergy%eqpsolver
                        stop
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