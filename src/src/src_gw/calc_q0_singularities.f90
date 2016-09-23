!--------------------------------------------------------------
! "General treatment of the singularities in Hartree-Fock and 
!  exact-exchange Kohn-Sham methods for solids",
!  by Carrier et al., PRB 75, 205126 (2007)
!--------------------------------------------------------------
subroutine calc_q0_singularities
    use modinput
    use modmain
    use modgw
    use mod_mpi_gw, only : myrank
    
    ! local variables
    integer :: iq, nz
    real(8) :: f, sum1, sumf2, intf1, intf2
    
    if (myrank==0) then
      write(*,*)
      write(*,*) 'Use the auxiliary function by Carrier et al. (2007)'
      write(*,*)
    end if
    
    sumf1 = 0.d0
    sumf2 = 0.d0
    do iq = 1, kqset%nkpt
      if (.not.gammapoint(kqset%vqc(1:3,iq))) then
        f = faux(kqset%vql(:,iq))
        sumf1 = sumf1 + dsqrt(f)
        sumf2 = sumf2 + f
      end if
    end do ! iq
    sumf1 = sumf1/dble(kqset%nkpt)
    sumf2 = sumf2/dble(kqset%nkpt)

    ! Gauss-Legendre quadrature
    call integrate_BZ_GL(256,256,256,intf1,intf2)
    
    ! adaptive grid
    !intf2 = integral_BZ(60,3,1.d-4)
    
    singc1 = intf1-sumf1
    singc2 = intf2-sumf2
    
    if (myrank==0) then
      write(*,*) 'Info(calc_q0_singularities): Integrals of the auxiliary function'
      write(*,1) 
      write(*,2) intf1, intf2
      write(*,3) sumf1, sumf2
      write(*,4) singc1, singc2
      1 format(30x,'q^(-1)',12x,'q^(-2)')
      2 format('Analitic integration: ',2f18.12)
      3 format('Numerical integration: ',2f18.12)
      4 format('Correction factor: ',2f18.12)
    end if
    
contains

    !---------------------------------------------------------
    ! Auxiliary function
    !---------------------------------------------------------
    real(8) function faux(vql)
        use modinput
        use modmain,     only: pi, twopi, bvec
        use mod_misc_gw, only: avec
        ! input parameters
        real(8), intent(in) :: vql(3)
        ! local variables
        real(8) :: vqc(3)
        
        faux = 4.d0*( dot_product(bvec(:,1),bvec(:,1))*dsin(pi*vql(1))*dsin(pi*vql(1)) &
        &            +dot_product(bvec(:,2),bvec(:,2))*dsin(pi*vql(2))*dsin(pi*vql(2)) &
        &            +dot_product(bvec(:,3),bvec(:,3))*dsin(pi*vql(3))*dsin(pi*vql(3)) &
        &           ) &
        &     +2.d0*( dot_product(bvec(:,1),bvec(:,2))*dsin(twopi*vql(1))*dsin(twopi*vql(2)) &
        &            +dot_product(bvec(:,2),bvec(:,3))*dsin(twopi*vql(2))*dsin(twopi*vql(3)) &
        &            +dot_product(bvec(:,3),bvec(:,1))*dsin(twopi*vql(3))*dsin(twopi*vql(1)) &
        &           )
        
        if (dabs(faux)<1.0d-16) then
          write(*,*) 'WARNING(calc_q0_singularities::faux)'
          write(*,*) '  Division by zero!'
          write(*,'(a,3f8.4)') '  vql=', vql
          faux = 0.d0
          !stop
        else
          faux = (2.d0*pi)**2 / faux
        end if
        
        return
    end function
    
    !----------------------------------------------------------------
    ! BZ integration: Gauss-Legendre quadrature
    !----------------------------------------------------------------
    subroutine integrate_BZ_GL(n1,n2,n3,intf1,intf2)
        use modmain, only : bvec
        integer, intent(in) :: n1, n2, n3
        real(8), intent(out):: intf1, intf2
        ! local
        integer :: i, j, k
        real(8) :: q1(n1), w1(n1)
        real(8) :: q2(n2), w2(n2)
        real(8) :: q3(n3), w3(n3)
        real(8) :: f, vql(3)
        
        ! generate grid
        call gauleg(-0.5d0,0.5d0,q1,w1,n1)
        call gauleg(-0.5d0,0.5d0,q2,w2,n2)
        call gauleg(-0.5d0,0.5d0,q3,w3,n3)
        
        ! integration
        intf1 = 0.d0
        intf2 = 0.d0
        do i = 1, n1
          vql(1) = q1(i)
          do j = 1, n2
            vql(2) = q2(j)
            do k = 1, n3
              vql(3) = q3(k)
              f = faux(vql)
              intf1 = intf1 + w1(i)*w2(j)*w3(k)*dsqrt(f)
              intf2 = intf2 + w1(i)*w2(j)*w3(k)*f
            end do
          end do
        end do
        
        return
    end subroutine
    
    !----------------------------------------------------------------
    ! BZ integration: Iterative algorithm using an adaptive scheme
    !----------------------------------------------------------------
    real(8) function integral_BZ(n0,div,eps)
        integer, intent(in) :: n0
        integer, intent(in) :: div
        real(8), intent(in) :: eps
        ! local variables
        integer :: iter, n_out, n, i, j, k, nn
        integer, allocatable :: i_grid(:,:)
        real(8) :: vql(3)
        real(8) :: sum, sum_prev, sum_out, delta
        logical :: tconv
        
        write(*,*)
        write(*,*) 'Info(calc_q0_singularities::integral_BZ):'
        write(*,*) ' Number of points: ', n0
        write(*,*) ' Accuracy: ', eps
        write(*,*)
        
        ! outer region grid only (grid is symmetric)
        ! Note also that BZ is set in [-b/2,b/2]
        n_out = (2*n0+1)**3 - (2*n0/div+1)**3
        allocate(i_grid(n_out,3))
        
        n = 0
        do i = -n0, n0
        do j = -n0, n0
        do k = -n0, n0 
          if ((abs(i) > n0/div).or. &
          &   (abs(j) > n0/div).or. &
          &   (abs(k) > n0/div)) then
            n = n+1
            i_grid(n,1) = i
            i_grid(n,2) = j
            i_grid(n,3) = k
            !write(*,*) 'i_grid=', i_grid(n,:)
          end if
        end do
        end do
        end do
        !write(*,*) n, n_out
        !stop
        
        n = n0
        sum = 0.d0
        sum_prev = 0.d0
        tconv = .false.
        
        do iter = 0, 10
        
          ! integral over the outer region
          sum_out = 0.d0
          do i = 1, n_out
            vql(1) = dble(i_grid(i,1))/dble(2*n)
            vql(2) = dble(i_grid(i,2))/dble(2*n)
            vql(3) = dble(i_grid(i,3))/dble(2*n)
            !write(*,*) 'vql=', vql
            sum_out = sum_out + faux(vql)
          end do ! i
          
          sum_out = sum_out / dble(2*n)**3
          
          ! update the BZ integral value 
          sum = sum + sum_out

          ! check for convergence
          delta = sum - sum_prev
          
          !write(*,*) 'iteration=', iter
          !write(*,*) 'sum_out=', sum_out
          !write(*,*) 'sum=', sum
          !write(*,*) 'delta=', delta
          !write(*,*)
          
          if (dabs(delta)<eps) then
            tconv = .true.
            exit
          else
            sum_prev = sum
            !------------------------------------------
            ! Next iteration step: triple the sampling 
            !------------------------------------------
            n = n*div
          end if
          
        end do ! iter
        
        deallocate(i_grid)
        
        if (.not.tconv) then
          write(*,*)
          write(*,*) 'WARNING(calc_q0_singularities::integral_BZ)'
          write(*,*) '  Integral over BZ is not converged!'
          write(*,*)
        end if
        
        integral_BZ = sum
        
    end function    
    
end subroutine
