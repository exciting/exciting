
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
    real(8),    parameter :: etol = 1.d-6
    integer(4) :: iter, ik, ib, nz
    real(8)    :: enk, egap, efdos, dzf2, scissor
    complex(8) :: ein, dsigma, znk, zf, dzf
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
            
            ! I think, the shift wrt efermi is important due to the requirement to
            ! have a self-consistency at efermi where Self_c = 0
            enk = evalsv(ib,ik) - efermi
            ein = cmplx(enk,0.d0,8)

            converged = .false.
            do iter = 1, nitermax

                ! print*, ik, ib, iter, ein

                if ( enk > 0.d0) then
                    f(:) = selfec(ib,:,ik)
                    call pade_approximant(nz, z, f, ein, sigc(ib,ik), dsigma)
                else
                    f(:) = conjg(selfec(ib,:,ik))
                    call pade_approximant(nz, conjg(z), f, ein, sigc(ib,ik), dsigma)
                end if

                if (input%gw%selfenergy%iopes == 0 ) then
                    ! Perturbative solution (single iteration)
                    znk = 1.0d0 / (1.0d0-dsigma)
                    znorm(ib,ik)  = dble(znk)
                    ein = enk + znk * ( selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik) )
                    converged = .true.
                    exit
                end if

                ! f(z) = E^{KS} + \Sigma_x + \Sigma_c(z) - Vxc - z
                zf = enk + selfex(ib,ik) + sigc(ib,ik) - vxcnn(ib,ik) - ein

                if ( abs(zf) < etol ) then
                    converged = .true.
                    exit
                else
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

            ! Quasiparticle energy (shifted back to have it in absolute units)
            evalqp(ib,ik) = dble(ein) + efermi

            write(fid,'(i4,4x,2f10.6,4x,f10.6,4x,2f10.6,4x,f10.6)') ib, ein+efermi, dble(selfex(ib,ik)), sigc(ib,ik), dble(vxcnn(ib,ik))

        end do ! ib

        write(fid,*) ; write(fid,*)

    end do ! ik
        
    ! Calculate Fermi energy
    call fermi_exciting(input%groundstate%tevecsv, &
    &                   nvelgw, &
    &                   nbandsgw, kset%nkpt, evalqp(ibgw:nbgw,:), &
    &                   kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
    &                   eferqp, egap, efdos)

    deallocate(z, f)
    close(fid)

end subroutine