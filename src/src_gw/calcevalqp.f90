!BOP
!
!!ROUTINE: calcevalqp
!
!!INTERFACE:
!
subroutine calcevalqp
!
!!DESCRIPTION:
!
! Given the matrix elements $\langle
! \Psi_{n\vec{k}}|\Sigma(\vec{k},\omega)|\Psi_{n\vec{k}}\rangle$, $\langle
! \Psi_{n\vec{k}}|V^{xc}|\Psi_{n\vec{k}}\rangle$ and
! $\varepsilon^{DFT}_{n\vec{k}}$. this subroutine calculates the
! quasi-particle energies $\varepsilon^{qp}_{n\vec{k}}$
!
!!USES:
    use modinput
    use modmain, only: efermi, zzero
    use modgw,   only: ibgw, nbgw, kset, evalfv, evalqp, eferqp, &
                       sigc, znorm, selfex, selfec, &
                       nbandsgw, nvelgw, &
                       sigsx, sigch, fgw
    use mod_vxc, only: vxcnn
    implicit none
    integer :: nb, ie, ik
    real(8) :: egap, df

!!REVISION HISTORY:
!
! Created: 16.08.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC
    select case (input%gw%taskname)

      case('g0w0','evalqp')
        call solve_QP_equation()

      case('g0w0-x')

        do ik = 1, kset%nkpt
          do ie = ibgw, nbgw
            evalqp(ie,ik) = evalfv(ie,ik) + &
                            dble(selfex(ie,ik) - vxcnn(ie,ik))
          end do ! ie
        end do ! ik

      case('cohsex')

        do ik = 1, kset%nkpt
          do ie = ibgw, nbgw
            sigsx(ie,ik) = sigsx(ie,ik)+selfex(ie,ik)
            evalqp(ie,ik) = evalfv(ie,ik) + &
                            dble(sigsx(ie,ik) + sigch(ie,ik) - vxcnn(ie,ik))
          end do ! ie
        end do ! ik

    end select

    ! Calculate Fermi energy
    call fermi_exciting(.false., &
    &                   nvelgw, &
    &                   nbandsgw, kset%nkpt, evalqp(ibgw:nbgw,:), &
    &                   kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
    &                   eferqp, egap, df)

end subroutine
!EOC
