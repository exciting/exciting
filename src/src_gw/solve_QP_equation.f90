
subroutine solve_QP_equation()
    
    use modinput
    use modmain, only: evalsv, efermi
    use modgw,   only: kset, ibgw, nbgw, nvelgw, nbandsgw, evalqp, eferqp, &
                       freq, selfex, selfec, sigc, znorm, freq_selfc
    use mod_vxc, only: vxcnn
    use mod_pade
    use m_getunit

    implicit none
    integer(4), parameter :: nitermax = 100
    real(8),    parameter :: etol = 1.d-4
    integer(4) :: iter, ik, ib, nz
    real(8)    :: enk, egap, efdos, dzf2, scissor
    complex(8) :: ein, ein_prev, dsigma, znk, zf, dzf
    complex(8), allocatable :: z(:), f(:)
    logical    :: converged
    integer    :: fid
    
    nz = freq_selfc%nomeg
    allocate(z(nz), f(nz))
    z(:) = cmplx( 0.d0, freq_selfc%freqs(:), 8)

    call getunit(fid)
    open(fid,file='EVALQP-CMPLX.DAT',action='WRITE',form='FORMATTED')
    
    do ik = 1, kset%nkpt

        write(fid,'(a,i4,4f12.6)') 'k-point #', ik, kset%vkl(:,ik), kset%wkpt(ik)
        write(fid,'(a,T11,a,T35,a,T49,a,T73,a)') 'state', 'E_QP', 'Self_x', 'Self_c', 'V_xc'

        do ib = ibgw, nbgw
            
            enk = evalsv(ib,ik)-efermi
            ein = cmplx(enk, 0.d0, 8)
            ein_prev = ein

            converged = .false.
            do iter = 1, nitermax

                ! print*, ''
                ! print*, 'iter=', iter
                ! print*, 'ik, ib, ein=', ik, ib, ein

                ! Analytical continuation
                call eval_ac(nz, z, selfec(ib,:,ik), ein, sigc(ib,ik), dsigma)
                ! print*, 'sigc=', sigc(ib,ik)
                ! print*, 'dsigma=', dsigma

                if (input%gw%selfenergy%iopes == 0) then
                    ! Perturbative solution (single iteration)
                    znk = zone / (zone-dsigma)
                    znorm(ib,ik)  = dble(znk)
                    ein = enk + znk * (selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik))
                    converged = .true.
                    exit
                else if (input%gw%selfenergy%iopes == 1) then
                    ! Perturbative solution without renormalization
                    znorm(ib,ik)  = 1.d0
                    ein = enk + selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik)
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

            write(fid,'(i4,4x,2f10.6,4x,f10.6,4x,2f10.6,4x,f10.6)') ib, ein, dble(selfex(ib,ik)), sigc(ib,ik), dble(vxcnn(ib,ik))

        end do ! ib

        write(fid,*) ; write(fid,*)

    end do ! ik

    close(fid)
    deallocate(z, f)
        
    ! Calculate Fermi energy
    call fermi_exciting(input%groundstate%tevecsv, &
    &                   nvelgw, &
    &                   nbandsgw, kset%nkpt, evalqp(ibgw:nbgw,:), &
    &                   kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
    &                   eferqp, egap, efdos)

    ! Shift QP energies to align KS and QP Fermi energies
    do ik = 1, kset%nkpt
        do ib = ibgw, nbgw
            evalqp(ib,ik) = evalqp(ib,ik) - eferqp + efermi
        end do
    end do    

end subroutine