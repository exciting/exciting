
subroutine solve_QP_equation()
    
    use modinput
    use modmain, only: evalsv, efermi
    use modgw,   only: kset, ibgw, nbgw, nvelgw, nbandsgw, evalqp, eferqp
    use mod_vxc, only: vxcnn
    use mod_selfenergy, only: selfex, selfec, sigc, znorm, freq_selfc 
    use mod_pade
    use m_getunit

    implicit none
    integer(4), parameter :: nitermax = 100
    real(8),    parameter :: etol = 1.d-4
    integer(4) :: iter, ik, ib, nz
    real(8)    :: enk, dzf2
    complex(8) :: ein, ein_prev, dsigma, znk, zf, dzf
    logical    :: converged
    
    do ik = 1, kset%nkpt
        do ib = ibgw, nbgw
            
            enk = evalsv(ib,ik)-efermi
            ein = cmplx(enk, 0.d0, 8)
            ein_prev = ein

            converged = .false.
            do iter = 1, nitermax
                
                ! get the value of Sigma_c at the energy e_nk
                call get_selfc(freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                dble(ein), sigc(ib,ik), dsigma)
                
                ! print*, ''
                ! print*, 'iter=', iter
                ! print*, 'ik, ib, ein=', ik, ib, ein
                ! print*, 'sigc=', sigc(ib,ik)
                ! print*, 'dsigma=', dsigma

                if (input%gw%selfenergy%iopes == 0) then
                    ! Perturbative solution (single iteration)
                    znk = zone / (zone-dsigma)
                    znorm(ib,ik)  = dble(znk)
                    ein = ein + znk * (selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik))
                    converged = .true.
                    exit
                else if (input%gw%selfenergy%iopes == 1) then
                    ! Perturbative solution without renormalization
                    znorm(ib,ik)  = 1.d0
                    ein = ein + selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik)
                end if

                ! Error function
                zf = ein - ein_prev

                if (abs(zf) < etol) then
                    converged = .true.
                    exit
                else
                    ! Prepare for next iteration
                    ein_prev = ein
                    ! apply next Newton-Raphson step
                    dzf  = dsigma - zone
                    dzf2 = abs(dzf)**2
                    ein  = ein - (zf*conjg(dzf))/dzf2
                end if

            end do ! iter
            
            if (.not.converged) then
                write(*,*)
                write(*,'(a,i0,a)') 'Warning(solve_QP_equation) Newton-Raphson solution of the quasiparticle equation is not converged after ', nitermax, ' iterations!'
                write(*,'(a,f8.4,a,f8.4)') 'Absolute error = ', abs(zf), ' > ', etol
            end if

            ! Quasiparticle energy
            evalqp(ib,ik) = dble(ein)

        end do ! ib
    end do ! ik

end subroutine