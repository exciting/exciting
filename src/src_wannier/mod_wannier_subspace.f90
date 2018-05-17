module mod_wannier_subspace
  use mod_wannier_variables
  use mod_wannier_maxloc
  use mod_wannier_opf
  use mod_wannier_helper
  use mod_eigenvalue_occupancy, only: efermi

  implicit none

! module variables
    
contains

    subroutine wannier_subspace
      use m_plotmat

      integer :: ik, ikb, i, nik, nikb, nok, nokb, n, is, ni, idxn, jst, fst, lst, it, maxit, iproj
      real(8) :: mixing, maxdiff, omegai0, t0, t1, fac, win(2)
      real(8) :: omegai, omegaik, omegaiinner
      complex(8) :: tmpz, zdotc

      integer, allocatable :: idxo(:) 
      real(8), allocatable :: sval(:), eval(:), evalmem(:,:), score(:,:), scoresum(:)
      complex(8), allocatable :: projm(:,:), z(:,:,:), z0(:,:,:), evec(:,:), lsvec(:,:), rsvec(:,:), auxmat(:,:), auxmat2(:,:), m(:,:,:,:), subspace(:,:,:), subspace_mem(:,:,:)

      ! find optimal set of nwf bands to generate initial guess from it via OPF method
      allocate( idxo( wf_groups( wf_group)%nwf))
      allocate( scoresum( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst))
      allocate( score( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_kset%nkpt))
      fac = 4.d0
      score = 0.d0
      idxo = 0
      win = wf_groups( wf_group)%win_i
      if( sum( abs( win)) .lt. 1.d-23) then
        win = wf_groups( wf_group)%win_o
      end if
      call wannier_geteval( evalmem, fst, lst)
      if( wf_fermizero) then
        call readfermi
        evalmem = evalmem - efermi
      end if
      do ik = 1, wf_kset%nkpt
        do i = 1, wf_groups( wf_group)%win_ni( ik)
          n = wf_groups( wf_group)%win_ii( i, ik)
          score( n, ik) = score( n, ik) + &
            1.d0/( 1.d0 + exp( -fac*( evalmem( n, ik) - win(1))/( win(2) - win(1)))) * &
            1.d0/( 1.d0 + exp(  fac*( evalmem( n, ik) - win(2))/( win(2) - win(1))))
        end do
        do i = 1, wf_groups( wf_group)%win_no( ik)
          n = wf_groups( wf_group)%win_io( i, ik)
          score( n, ik) = score( n, ik) + &
            1.d0/( 1.d0 + exp( -fac*( evalmem( n, ik) - win(1))/( win(2) - win(1)))) * &
            1.d0/( 1.d0 + exp(  fac*( evalmem( n, ik) - win(2))/( win(2) - win(1))))
        end do
      end do
      scoresum = sum( score, 2)/wf_kset%nkpt
      i = 0
      do while( i .lt. wf_groups( wf_group)%nwf)
        n = wf_groups( wf_group)%fst + maxloc( scoresum, 1) - 1
        i = i + 1
        idxo( i) = n
        scoresum( n) = 0.d0
      end do
      wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :) = zzero
      do i = 1, wf_groups( wf_group)%nwf
        wf_transform( idxo( i), wf_groups( wf_group)%fwf+i-1, :) = zone
      end do
      deallocate( evalmem)
        
      maxit = 10000

      ! search for nwf OPFs
      call wannier_loc
      omegai0 = sum( wf_omega_i( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      call wannier_gen_opf

      allocate( projm( wf_groups( wf_group)%nst, wf_groups( wf_group)%nprojused), &
                m( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst, wf_kset%nkpt, 2*wf_n_ntot))
      allocate( lsvec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst), &
                rsvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf), &
                sval( wf_groups( wf_group)%nst))
      allocate( subspace( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( subspace_mem( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( auxmat( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      allocate( auxmat2( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst))
      m = zzero
      
      ! get initial subspace
      do ik = 1, wf_kset%nkpt
        !write(*,*) ik
        nok = wf_groups( wf_group)%win_no( ik)
        nik = wf_groups( wf_group)%win_ni( ik)

        projm = zzero
        subspace( :, :, ik) = zzero

        do n = 1, wf_groups( wf_group)%win_ni( ik)
          ni = 0
          do iproj = 1, wf_nprojtot
            if( wf_groups( wf_group)%projused( iproj) .eq. 1) then
              ni = ni + 1
              projm( n, ni) = wf_projection( wf_groups( wf_group)%win_ii( n, ik), iproj, ik)
            end if
          end do
        end do
        do n = 1, wf_groups( wf_group)%win_no( ik)
          ni = 0
          do iproj = 1, wf_nprojtot
            if( wf_groups( wf_group)%projused( iproj) .eq. 1) then
              ni = ni + 1
              projm( wf_groups( wf_group)%win_ni( ik)+n, ni) = wf_projection( wf_groups( wf_group)%win_io( n, ik), iproj, ik)
            end if
          end do
        end do
        
        auxmat = zzero
        call zgemm( 'n', 'n', nok+nik, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nprojused, zone, &
               projm, wf_groups( wf_group)%nst, &
               wf_opf, wf_groups( wf_group)%nprojused, zzero, & 
               auxmat, wf_groups( wf_group)%nst)
        call zgesdd_wrapper( auxmat( 1:(nok+nik), :), &
               nok+nik, wf_groups( wf_group)%nwf, &
               sval( 1:wf_groups( wf_group)%nwf), &
               lsvec( 1:(nok+nik), 1:(nok+nik)), &
               rsvec)
        call zgemm( 'n', 'n', nok+nik, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               lsvec, wf_groups( wf_group)%nst, &
               rsvec, wf_groups( wf_group)%nwf, zzero, &
               auxmat, wf_groups( wf_group)%nst)
        call zgemm( 'n', 'c', nok+nik, nok+nik, wf_groups( wf_group)%nwf, zone, &
               auxmat, wf_groups( wf_group)%nst, &
               auxmat, wf_groups( wf_group)%nst, zzero, &
               auxmat2, wf_groups( wf_group)%nst)

        call zgemm( 'n', 'c', wf_groups( wf_group)%nst, wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, zone, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst, zzero, &
               lsvec, wf_groups( wf_group)%nst)
        !call plotmat( lsvec)
        !write(*,*)
        !auxmat2 = zzero
        !do n = 1, nk
        !  do i = 1, nk
        !    auxmat2( mk+n, mk+i) = lsvec( wf_groups( wf_group)%win_io( n, ik)-wf_groups( wf_group)%fst+1, wf_groups( wf_group)%win_io( i, ik)-wf_groups( wf_group)%fst+1)
        !  end do
        !end do

        lsvec = zzero
        call diaghermat( nok, auxmat2( (nik+1):(nok+nik), (nik+1):(nok+nik)), sval( 1:nok), lsvec( 1:nok, 1:nok))
        !write(*,*) nk, mk, wf_groups( wf_group)%nwf - mk
        !call plotmat( auxmat2( (mk+1):(nk+mk), (mk+1):(nk+mk)))
        !write(*,*)
        !write(*,'(1000f13.6)') sval( 1:nk)
        !call plotmat( lsvec(1:nk,1:nk))
        !write(*,*)
        !stop
        do n = 1, wf_groups( wf_group)%nwf - nik
          subspace( 1:nok, n, ik) = lsvec( 1:nok, nok-n+1)
        end do
        !call plotmat( subspace( :, :, ik))
        !write(*,*)

        !subspace(:,:,ik) = zzero
        !do n = 1, wf_groups( wf_group)%nwf
        !  i = wf_groups( wf_group)%fst + maxloc( score( :, ik), 1) - 1
        !  score( i, ik) = 0.d0
        !  if( n .gt. mk) then
        !    do ni = 1, nk
        !      if( wf_groups( wf_group)%win_io( ni, ik) .eq. i) exit
        !    end do
        !    subspace( ni, n-mk, ik) = zone
        !  end if
        !end do
      end do

      !subspace = zzero
      !do n = 1, wf_groups( wf_group)%nwf
      !  subspace( n, n, :) = zone
      !end do

      ! set up plane-wave matrix-elements
      m = zzero
      do ik = 1, wf_kset%nkpt
        nik = wf_groups( wf_group)%win_ni( ik)
        nok = wf_groups( wf_group)%win_no( ik)
        do idxn = 1, wf_n_ntot
          ! positive neighbors
          ikb = wf_n_ik( idxn, ik)
          nikb = wf_groups( wf_group)%win_ni( ikb)
          nokb = wf_groups( wf_group)%win_no( ikb)
          do i = 1, nik
            do n = 1, nikb
              m( i, n, ik, idxn) = wf_m0( wf_groups( wf_group)%win_ii( i, ik), wf_groups( wf_group)%win_ii( n, ikb), ik, idxn)
            end do
            do n = 1, nokb
              m( i, nikb+n, ik, idxn) = wf_m0( wf_groups( wf_group)%win_ii( i, ik), wf_groups( wf_group)%win_io( n, ikb), ik, idxn)
            end do
          end do
          do i = 1, nok
            do n = 1, nikb
              m( nik+i, n, ik, idxn) = wf_m0( wf_groups( wf_group)%win_io( i, ik), wf_groups( wf_group)%win_ii( n, ikb), ik, idxn)
            end do
            do n = 1, nokb
              m( nik+i, nikb+n, ik, idxn) = wf_m0( wf_groups( wf_group)%win_io( i, ik), wf_groups( wf_group)%win_io( n, ikb), ik, idxn)
            end do
          end do
          ! negative neighbors
          ikb = wf_n_ik2( idxn, ik)
          nikb = wf_groups( wf_group)%win_ni( ikb)
          nokb = wf_groups( wf_group)%win_no( ikb)
          do i = 1, nik
            do n = 1, nikb
              m( i, n, ik, wf_n_ntot+idxn) = conjg( wf_m0( wf_groups( wf_group)%win_ii( n, ikb), wf_groups( wf_group)%win_ii( i, ik), ikb, idxn))
            end do
            do n = 1, nokb
              m( i, nikb+n, ik, wf_n_ntot+idxn) = conjg( wf_m0( wf_groups( wf_group)%win_io( n, ikb), wf_groups( wf_group)%win_ii( i, ik), ikb, idxn))
            end do
          end do
          do i = 1, nok
            do n = 1, nikb
              m( nik+i, n, ik, wf_n_ntot+idxn) = conjg( wf_m0( wf_groups( wf_group)%win_ii( n, ikb), wf_groups( wf_group)%win_io( i, ik), ikb, idxn))
            end do
            do n = 1, nokb
              m( nik+i, nikb+n, ik, wf_n_ntot+idxn) = conjg( wf_m0( wf_groups( wf_group)%win_io( n, ikb), wf_groups( wf_group)%win_io( i, ik), ikb, idxn))
            end do
          end do
        end do
      end do

      ! calculate Z_0
      allocate( z0( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst, wf_kset%nkpt))
      z0 = zzero
      do ik = 1, wf_kset%nkpt
        nik = wf_groups( wf_group)%win_ni( ik)
        nok = wf_groups( wf_group)%win_no( ik)
        do idxn = 1, wf_n_ntot
          ! positive neighbors
          ikb = wf_n_ik( idxn, ik)
          nikb = wf_groups( wf_group)%win_ni( ikb)
          if( nikb .gt. 0) then
            call zgemm( 'n', 'c', nok, nok, nikb, cmplx( wf_n_wgt( idxn), 0.d0, 8), &
                   m( nik+1, 1, ik, idxn), wf_groups( wf_group)%nst, &
                   m( nik+1, 1, ik, idxn), wf_groups( wf_group)%nst, zone, &
                   z0( :, :, ik), wf_groups( wf_group)%nst)
          end if
          ! negative neighbors
          ikb = wf_n_ik2( idxn, ik)
          nikb = wf_groups( wf_group)%win_ni( ikb)
          if( nikb .gt. 0) then
            call zgemm( 'n', 'c', nok, nok, nikb, cmplx( wf_n_wgt( idxn), 0.d0, 8), &
                   m( nik+1, 1, ik, wf_n_ntot+idxn), wf_groups( wf_group)%nst, &
                   m( nik+1, 1, ik, wf_n_ntot+idxn), wf_groups( wf_group)%nst, zone, &
                   z0( :, :, ik), wf_groups( wf_group)%nst)
          end if
        end do
      end do

      ! idependent contribution to Omega_I from inner-window states
      omegaiinner = 0.d0
      do ik = 1, wf_kset%nkpt
        nik = wf_groups( wf_group)%win_ni( ik)
        do idxn = 1, wf_n_ntot
          ikb = wf_n_ik( idxn, ik)
          nikb = wf_groups( wf_group)%win_ni( ikb)
          do i = 1, nik
            omegaiinner = omegaiinner - wf_n_wgt( idxn)*dble( zdotc( nikb, m( i, 1:nikb, ik, idxn), 1, m( i, 1:nikb, ik, idxn), 1))
          end do
          ikb = wf_n_ik2( idxn, ik)
          nikb = wf_groups( wf_group)%win_ni( ikb)
          do i = 1, nik
            omegaiinner = omegaiinner - wf_n_wgt( idxn)*dble( zdotc( nikb, m( i, 1:nikb, ik, wf_n_ntot+idxn), 1, m( i, 1:nikb, ik, wf_n_ntot+idxn), 1))
          end do
        end do
      end do
 
      ! perfom disentanglement
      allocate( evec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst), &
                eval( wf_groups( wf_group)%nst), &
                evalmem( wf_groups( wf_group)%nst, wf_kset%nkpt))
      allocate( z( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst, wf_kset%nkpt))
      z = zzero
      subspace_mem = subspace
      maxdiff = 1.d0
      evalmem = 0.d0
      it = 0
      write( wf_info, '(" disentangle optimal subspace...")')
      call timesec( t0)
      do while( (maxdiff .gt. input%properties%wannier%grouparray( wf_group)%group%epsdis) .and. (it .lt. maxit))

        !if( mod( it, 100) .eq. 0) then
        !  !call plotmat( z( :, :, 12))
        !  wf_transform( :, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :) = zzero
        !  do ik = 1, wf_kset%nkpt
        !    do n = 1, wf_groups( wf_group)%win_ni( ik)
        !      wf_transform( wf_groups( wf_group)%win_ii( n, ik), wf_groups( wf_group)%fwf+n-1, ik) = zone
        !    end do
        !    do n = 1, wf_groups( wf_group)%win_no( ik)
        !      do i = 1, wf_groups( wf_group)%nwf - wf_groups( wf_group)%win_ni( ik)
        !        wf_transform( wf_groups( wf_group)%win_io( n, ik), wf_groups( wf_group)%fwf+wf_groups( wf_group)%win_ni( ik)+i-1, ik) = subspace( n, i, ik)
        !      end do
        !    end do
        !  end do
        !  call wannier_loc

        !  write(*,'(i,g13.6,2f13.6,g13.6)') it, maxdiff, sum( wf_omega_i( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))!, omegai, omegai/sum( wf_omega_i( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))-1.d0
        !end if
        
        it = it + 1
        maxdiff = 0.d0
        mixing = 0.9d0
        omegai = 0.d0
        if( it .eq. 1) mixing = 1.d0
#ifdef USEOMP
!$omp parallel default( shared) private( ik, ikb, i, n, nik, nok, nikb, nokb, idxn, auxmat, eval, evec, omegaik, tmpz) reduction(+:omegai)
!$omp do
#endif
        do ik = 1, wf_kset%nkpt
          omegaik = 2.d0*sum( wf_n_wgt( 1:wf_n_ntot))*wf_groups( wf_group)%nwf
          nik = wf_groups( wf_group)%win_ni( ik)
          nok = wf_groups( wf_group)%win_no( ik)
          
          if( nok .ge. 1) then
            z( :, :, ik) = (1.d0 - mixing)*z( :, :, ik) + mixing*z0( :, :, ik)
            do idxn = 1, wf_n_ntot
              ikb = wf_n_ik( idxn, ik)
              nikb = wf_groups( wf_group)%win_ni( ikb)
              nokb = wf_groups( wf_group)%win_no( ikb)
              n = wf_groups( wf_group)%nwf-nikb
              if( nokb .ge. 1) then
                call zgemm( 'n', 'n', nik+nok, n, nokb, zone, &
                       m( 1, nikb+1, ik, idxn), wf_groups( wf_group)%nst, &
                       subspace( 1, 1, ikb), wf_groups( wf_group)%nst, zzero, &
                       auxmat, wf_groups( wf_group)%nst)
                do i = 1, nik
                  omegaik = omegaik - wf_n_wgt( idxn)*dble( zdotc( n, auxmat( i, :), 1, auxmat( i, :), 1))
                end do
                call zgemm( 'n', 'c', nok, nok, n, cmplx( mixing*wf_n_wgt( idxn), 0, 8), &
                       auxmat( nik+1, 1), wf_groups( wf_group)%nst, &
                       auxmat( nik+1, 1), wf_groups( wf_group)%nst, zone, &
                       z( 1, 1, ik), wf_groups( wf_group)%nst)
              end if

              ikb = wf_n_ik2( idxn, ik)
              nikb = wf_groups( wf_group)%win_ni( ikb)
              nokb = wf_groups( wf_group)%win_no( ikb)
              n = wf_groups( wf_group)%nwf-nikb
              if( nokb .ge. 1) then
                call zgemm( 'n', 'n', nik+nok, n, nokb, zone, &
                       m( 1, nikb+1, ik, wf_n_ntot+idxn), wf_groups( wf_group)%nst, &
                       subspace( 1, 1, ikb), wf_groups( wf_group)%nst, zzero, &
                       auxmat, wf_groups( wf_group)%nst)
                do i = 1, nik
                  omegaik = omegaik - wf_n_wgt( idxn)*dble( zdotc( n, auxmat( i, :), 1, auxmat( i, :), 1))
                end do
                call zgemm( 'n', 'c', nok, nok, n, cmplx( mixing*wf_n_wgt( idxn), 0, 8), &
                       auxmat( nik+1, 1), wf_groups( wf_group)%nst, &
                       auxmat( nik+1, 1), wf_groups( wf_group)%nst, zone, &
                       z( 1, 1, ik), wf_groups( wf_group)%nst)
              end if

            end do
            n = wf_groups( wf_group)%nwf - nik
            call diaghermat( nok, z( 1:nok, 1:nok, ik), eval( 1:nok), evec( 1:nok, 1:nok))
            !write(*,'(1000f13.6)') sum( wf_n_wgt( 1:wf_n_ntot)), eval( 1:nk)
            do i = 1, n
              subspace_mem( 1:nok, i, ik) = evec( 1:nok, nok-i+1)
              omegaik = omegaik - eval( nok-i+1)
            end do
#ifdef USEOMP
!$omp atomic update
#endif
            maxdiff = max( maxdiff, dsqrt( sum( (eval( (nok-n+1):nok) - evalmem( (nok-n+1):nok, ik))**2)/n))
#ifdef USEOMP
!$omp end atomic
#endif
            evalmem( 1:nok, ik) = eval( 1:nok)
          end if

          omegai = omegai + omegaik

        end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
        omegai = (omegai + omegaiinner)/wf_kset%nkpt
        write(*,'(1a1,i8,f23.16,$)') char(13), it, omegai

        subspace = subspace_mem
      end do
      write(*,*)

      wf_transform( :, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :) = zzero
      do ik = 1, wf_kset%nkpt
        do n = 1, wf_groups( wf_group)%win_ni( ik)
          wf_transform( wf_groups( wf_group)%win_ii( n, ik), wf_groups( wf_group)%fwf+n-1, ik) = zone
        end do
        do n = 1, wf_groups( wf_group)%win_no( ik)
          do i = 1, wf_groups( wf_group)%nwf - wf_groups( wf_group)%win_ni( ik)
            wf_transform( wf_groups( wf_group)%win_io( n, ik), wf_groups( wf_group)%fwf+wf_groups( wf_group)%win_ni( ik)+i-1, ik) = subspace( n, i, ik)
          end do
        end do
      end do

      !ik = 12
      !write(*,*) wf_groups( wf_group)%win_ii( :, ik)
      !write(*,*) wf_groups( wf_group)%win_io( :, ik)
      !write(*,'(100f13.6)') score( :, ik)
      !call plotmat( wf_transform( :, :, ik), matlab=.true.)
      !write(*,*)

      deallocate( auxmat, m, subspace, evec, evalmem, eval, lsvec, rsvec, sval, idxo, projm, z, z0, subspace_mem)
      call wannier_loc
      call timesec( t1)

      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') it
      write( wf_info, '(5x,"convergence cutoff: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%epsdis
      write( wf_info, '(5x,"Omega_I before: ",T40,F13.6)') omegai0
      write( wf_info, '(5x,"Omega_I after: ",T40,F13.6)') sum( wf_omega_i( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      write( wf_info, '(5x,"reduction: ",T40,7x,I5,"%")') nint( 100.d0*(omegai0-sum( wf_omega_i( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)))/omegai0)
      write( wf_info, *)
      call flushifc( wf_info)

    end subroutine wannier_subspace

end module mod_wannier_subspace
