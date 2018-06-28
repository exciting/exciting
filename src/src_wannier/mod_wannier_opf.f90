module mod_wannier_opf
  use mod_wannier_variables
  use mod_wannier_maxloc, only: wannier_loc, wannier_phases
  use m_linalg

  use mod_constants,      only: zone, zzero, zi, twopi
  use modinput

  implicit none

! module variables
    
contains
    !BOP
    ! !ROUTINE: wannier_gen_opf
    ! !INTERFACE:
    !
    subroutine wannier_gen_opf( maxit, nowrite)
      ! !USES:
      use m_plotmat
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      ! !DESCRIPTION:
      !   Does the same thing as {\tt genwf} does but the transformation
      !   matrices are used for the generation of maximally localized Wannier
      !   functions. The matrices are computed in a self consistent loop. 
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      integer, optional, intent( in) :: maxit
      logical, optional, intent( in) :: nowrite
      ! local variables
      integer :: iknr, is, n, nx, iproj, i, j, minit, idxn, ist, a
      real(8) :: lambda, sum_n_wgt, phi, phimean, uncertainty, theta, evalmin, omegamin, x1
      real(8) :: opf_min_q(3,3), opf_min_p(3,1), opf_min_c, opf_min_sol(3), eval(3), revec(3,3), opf_min_big(6,6), opf_min_aux(3,3), r, tp(2)
      real(8) :: t0, t1, t2
      complex(8) :: opf_min_z(3,1), evec_big(6,6), eval_big(6), r1, r2, tmp
      logical :: opf_nowrite

      ! allocatable arrays
      complex(8), allocatable :: projection(:,:), opf_transform(:,:,:), opf_projm(:,:,:), lsvec(:,:), rsvec(:,:), opf_mixing(:,:), opf_x(:,:,:), opf_x0(:,:,:), opf_m(:,:)
      complex(8), allocatable :: auxmat(:,:), projm(:,:), lsvec2(:,:), rsvec2(:,:)
      real(8), allocatable :: sval(:), opf_t(:), sval2(:), phihist(:)

      opf_nowrite = .false.
      if( present( nowrite)) opf_nowrite = nowrite

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(" calculate optimized projection functions (OPF)...")')
      call timesec( t0)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      !********************************************************************
      ! initialize M matrices and Wannier centers
      !********************************************************************
      call wannier_loc

      !********************************************************************
      ! build rectangular transformation matrices
      !********************************************************************
      allocate( opf_transform( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, wf_kset%nkpt))
      allocate( opf_projm( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, wf_kset%nkpt))
      allocate( lsvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf), &
                rsvec( wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nprojused), &
                sval( wf_groups( wf_group)%nwf))
      allocate( projection( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%nprojused))

      do iknr = 1, wf_kset%nkpt
        projection = zzero
        n = 0
        do iproj = 1, wf_nprojtot
          if( wf_groups( wf_group)%projused( iproj) .eq. 1) then
            n = n + 1
            projection( :, n) = wf_projection( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, iproj, iknr)
          end if
        end do
        call zgemm( 'c', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nst, zone, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, iknr), wf_nst, &
               projection, wf_groups( wf_group)%nst, zzero, &
               opf_projm( :, :, iknr), wf_groups( wf_group)%nwf)
        call zsvd( opf_projm( :, :, iknr), sval, lsvec, rsvec)
        ! for numerical stability
        !do ist = 1, wf_nst
        !  if( sval( ist) .lt. 1.d-12) lsvec( :, ist) = zzero
        !end do
        
        call zgemm( 'n', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nwf, zone, &
               lsvec, wf_groups( wf_group)%nwf, &
               rsvec, wf_groups( wf_group)%nprojused, zzero, &
               opf_transform( :, :, iknr), wf_groups( wf_group)%nwf)
      end do
      deallocate( projection)

      !call writematlab( opf_projm(:,:,12), 'pmopf')
      !call writematlab( opf_transform(:,:,12), 'tmopf')

      sum_n_wgt = 2.d0*sum( wf_n_wgt)
      lambda = sum_n_wgt*input%properties%wannier%grouparray( wf_group)%group%lambdaopf
      minit = 2

      !********************************************************************
      ! build enlarged overlap and constraint matrices
      !********************************************************************
      allocate( opf_x( wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nprojused, wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_x0( wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nprojused, wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_t( wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_m( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))

      ! build enlarged overlap matrices
      allocate( auxmat( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused))
      i = 0
      do iknr = 1, wf_kset%nkpt
        do idxn = 1, wf_n_ntot
          i = i + 1
          !do ist = 1, opf_nst
          !  do jst = 1, opf_nst
          !    opf_m( ist, jst) = wf_m0( opf_bands( ist), opf_bands( jst), iknr, idxn)
          !  end do
          !end do
          call zgemm( 'n', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nwf, zone, &
                 wf_m( wf_groups( wf_group)%fwf, wf_groups( wf_group)%fwf, iknr, idxn), wf_nwf, &
                 opf_transform( :, :, wf_n_ik( idxn, iknr)), wf_groups( wf_group)%nwf, zzero, &
                 auxmat, wf_groups( wf_group)%nwf)
          call zgemm( 'c', 'n', wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nwf, zone, &
                 opf_transform( :, :, iknr), wf_groups( wf_group)%nwf, &
                 auxmat, wf_groups( wf_group)%nwf, zzero, &
                 opf_x0( :, :, i), wf_groups( wf_group)%nprojused)
          opf_t( i) = -2.d0*wf_n_wgt( idxn)/wf_kset%nkpt
        end do
      end do
      ! build constraint matrices
      do iknr = 1, wf_kset%nkpt
        i = i + 1
        call zgemm( 'c', 'n', wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nwf, zone, &
               opf_projm( :, :, iknr), wf_groups( wf_group)%nwf, &
               opf_projm( :, :, iknr), wf_groups( wf_group)%nwf, zzero, &
               opf_x0( :, :, i), wf_groups( wf_group)%nprojused)
        do is = 1, wf_groups( wf_group)%nprojused
          opf_x0( is, is, i) = opf_x0( is, is, i) - zone
        end do
        opf_t( i) = lambda/wf_kset%nkpt
      end do
      nx = i

      !********************************************************************
      ! minimize Lagrangian
      !********************************************************************
      if( allocated( wf_opf)) deallocate( wf_opf)
      allocate( wf_opf( wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nwf))

      allocate( opf_mixing( wf_groups( wf_group)%nprojused, wf_groups( wf_group)%nprojused))
      allocate( phihist( minit))
      deallocate( lsvec, rsvec)
      allocate( lsvec(3,3), rsvec(3,3))
      omegamin = 1.d13
      opf_mixing = zzero
      opf_x = opf_x0
      do i = 1, wf_groups( wf_group)%nprojused
        opf_mixing( i, i) = zone
      end do
      n = 0
      phihist = 0.d0
      uncertainty = 1.d0
      ! start minimization
      call timesec( t2)
      !write(*,*) "opf: start minimization"
      !write(*,'(" nx = ",I6)') nx
      !write(*,'(" omega = ",4F13.6)') sum( wf_omega   ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                                sum( wf_omega_i ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                                sum( wf_omega_d ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                                sum( wf_omega_od( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      do while( (n .lt. minit) .or. (uncertainty .gt. input%properties%wannier%grouparray( wf_group)%group%epsopf))
        n = n + 1
        if( present( maxit)) then
          if( n .gt. maxit) exit
        end if
        do i = 1, wf_groups( wf_group)%nwf
          do j = i+1, wf_groups( wf_group)%nprojused
            opf_min_q = 0.d0
            opf_min_p = 0.d0
            opf_min_c = 0.d0
            do is = 1, nx
              opf_min_z(1,1) = 0.5d0*(opf_x( i, i, is) - opf_x( j, j, is))
              opf_min_z(2,1) = -0.5d0*(opf_x( i, j, is) + opf_x( j, i, is))
              opf_min_z(3,1) = 0.5d0*zi*(opf_x( i, j, is) - opf_x( j, i, is))
              opf_min_q = opf_min_q + opf_t( is)*dble( matmul( opf_min_z, transpose( conjg( opf_min_z))))
              opf_min_p = opf_min_p + opf_t( is)*dble( conjg( opf_x( i, i, is) + opf_x( j, j, is))*opf_min_z)
              !opf_min_c = 0.25d0*opf_t( is)*abs( opf_x( i, i, is) + opf_x( j, j, is))**2
            end do

            call rsydiag( opf_min_q, eval, revec)
            if( j .le. wf_groups( wf_group)%nwf) then
              opf_min_sol(:) = revec( :, 1)
              if( opf_min_sol(1) .lt. 0.d0) opf_min_sol = -opf_min_sol
            else
              opf_min_big = 0.d0
              opf_min_big( 1:3, 1:3) = 2*opf_min_q
              opf_min_big( 1:3, 4:6) = 0.25d0*matmul( opf_min_p, transpose( opf_min_p)) - matmul( opf_min_q, opf_min_q)
              opf_min_big( 4, 1) = 1.d0
              opf_min_big( 5, 2) = 1.d0
              opf_min_big( 6, 3) = 1.d0
              call zgediag( cmplx( opf_min_big, 0, 8), eval_big, evec_big)
              evalmin = maxval( abs( eval_big))
              do is = 1, 6
                if( abs( aimag( eval_big( is))) .lt. input%structure%epslat) evalmin = min( evalmin, dble( eval_big( is)))
              end do
              do is = 1, 3
                opf_min_q( is, is) = opf_min_q( is, is) - evalmin
              end do
              if( minval( eval(:) - evalmin) .lt. 1.d-6) then
                call rpinv( opf_min_q, opf_min_aux)
                call r3mv( opf_min_aux, -0.5d0*opf_min_p(:,1), opf_min_sol)
                call r3mv( opf_min_q, opf_min_sol, opf_min_aux(:,1))
                if( (norm2( opf_min_aux(:,1) + 0.5d0*opf_min_p(:,1)) .ge. 1.d-6) .or. (norm2( opf_min_sol) .gt. 1.d0)) then
                  opf_min_sol = 10.d0
                else
                  opf_min_sol = opf_min_sol + revec(:,1)/norm2( revec(:,1))*sqrt( 1.d0 - norm2( opf_min_sol))
                end if
              else
                call r3minv( opf_min_q, opf_min_aux)
                call r3mv( opf_min_aux, -0.5d0*opf_min_p(:,1), opf_min_sol)
              end if
            end if
            if( abs( 1.d0 - norm2( opf_min_sol)) .le. input%structure%epslat) then
              call sphcrd( (/opf_min_sol(2), opf_min_sol(3), opf_min_sol(1)/), r, tp)
              x1 = tp(2)
              theta = tp(1)/2
              r1 = cmplx( cos( theta), 0, 8)
              r2 = exp( zi*x1)*sin( theta)
              ! update mixing matrix
              do a = 1, wf_groups( wf_group)%nprojused
                tmp = opf_mixing( a, i)
                opf_mixing( a, i) = tmp*r1 - opf_mixing( a, j)*conjg( r2)
                opf_mixing( a, j) = opf_mixing( a, j)*r1 + tmp*r2
              end do
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, tmp, a)
!$OMP DO  
#endif
              ! update X-matrices
              do is = 1, nx
                do a = 1, wf_groups( wf_group)%nprojused
                  tmp = opf_x( a, i, is)
                  opf_x( a, i, is) = tmp*r1 - opf_x( a, j, is)*conjg( r2)
                  opf_x( a, j, is) = opf_x( a, j, is)*r1 + tmp*r2
                end do
                do a = 1, wf_groups( wf_group)%nprojused
                  tmp = opf_x( i, a, is)
                  opf_x( i, a, is) = tmp*r1 - opf_x( j, a, is)*r2
                  opf_x( j, a, is) = opf_x( j, a, is)*r1 + tmp*conjg( r2)
                end do
              end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            else
              !write(*,*) "Doh!"
            end if
          end do
        end do
        ! calculate spread
        wf_opf(:,:) = opf_mixing( :, 1:wf_groups( wf_group)%nwf)

        phi = 0.d0
        do is = 1, nx
          do ist = 1, wf_groups( wf_group)%nwf
            phi = phi + opf_t( is)*opf_x( ist, ist, is)*conjg( opf_x( ist, ist, is))
          end do
        end do
        phi = phi + wf_groups( wf_group)%nwf*sum_n_wgt

        ! convergence analysis
        if( n .eq. 1) phihist(:) = phi
        phihist = cshift( phihist, -1)
        phihist(1) = phi
        phimean = sum( phihist(:))/minit
        if( n .gt. 1) then
          !uncertainty = sqrt( sum( (phihist(:)-phimean)**2)/(min( n, minit)-1))/abs( phi)
          uncertainty = abs( phihist(2)/phihist(1) - 1.d0)
        else
          uncertainty = 1.d0
        end if

        call timesec( t1)
        !write(*,'(I7,3x,40F23.16)') n, t1-t2, phi, sum( wf_omega), sum( wf_omega_i), sum( wf_omega_d), sum( wf_omega_od), uncertainty
        !write(*,'(i7,3x,40g23.16)') n, t1-t2, phi, uncertainty

      end do

      allocate( projm( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))
      allocate( sval2( wf_groups( wf_group)%nwf))
      allocate( lsvec2( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))
      allocate( rsvec2( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))
      deallocate( auxmat)
      allocate( auxmat( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%nwf))
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, projm, sval2, lsvec2, rsvec2, auxmat)
!$OMP DO  
#endif
      do iknr = 1, wf_kset%nkpt
        call zgemm( 'n', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, zone, &
               opf_projm( :, :, iknr), wf_groups( wf_group)%nwf, &
               wf_opf, wf_groups( wf_group)%nprojused, zzero, &
               projm, wf_groups( wf_group)%nwf)
        call zsvd( projm, sval2, lsvec2, rsvec2)
        ! for numerical stability
        !do ist = 1, wf_nst
        !  if( sval2( ist) .lt. 1.d-12) lsvec2( :, ist) = zzero
        !end do
        !write(*,*) iknr
        call zgemm( 'n', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               lsvec2, wf_groups( wf_group)%nwf, &
               rsvec2, wf_groups( wf_group)%nwf, zzero, &
               projm, wf_groups( wf_group)%nwf)
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, iknr), wf_nst, &
               projm, wf_groups( wf_group)%nwf, zzero, &
               auxmat, wf_groups( wf_group)%nst)
        wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, iknr) = auxmat
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      ! update M matrices
      call wannier_loc
      ! determine phases for logarithms
      call wannier_phases
      ! calculate spread
      call wannier_loc

      !write(*,'(I7,3x,40F23.16)') n, t1-t2, phi, &
      !                            sum( wf_omega ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                            sum( wf_omega_i ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                            sum( wf_omega_d ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                            sum( wf_omega_od( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)), &
      !                            uncertainty

      call timesec( t1)
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"minimum iterations: ",T40,7x,I6)') minit
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') n
      write( wf_info, '(5x,"cutoff uncertainty: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%epsopf
      write( wf_info, '(5x,"achieved uncertainty: ",T40,E13.6)') uncertainty
      write( wf_info, '(5x,"Lagrangian: ",T40,F13.6)') phi
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      write( wf_info, *)
      call flushifc( wf_info)

!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      deallocate( lsvec2, rsvec2, projm, opf_x, auxmat, phihist, opf_mixing, opf_x0, opf_t, opf_transform, lsvec, rsvec, sval, sval2)
      return
      !EOC
    end subroutine wannier_gen_opf
    !EOP

end module mod_wannier_opf
