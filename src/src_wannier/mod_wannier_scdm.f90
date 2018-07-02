module mod_wannier_scdm
  use mod_wannier_variables
  use mod_wannier_helper,       only: wannier_geteval, wannier_getevec
  use mod_wannier_maxloc,       only: wannier_loc
  use m_linalg

  use mod_APW_LO,               only: apwordmax
  use mod_atoms,                only: natmtot, nspecies, natoms, idxas, spr 
  use mod_constants,            only: twopi, zzero, zone
  use mod_muffin_tin,           only: lmmaxapw, nrmtmax, nrmt
  use mod_Gvector,              only: ngrtot
  use mod_Gkvector,             only: ngkmax_ptr
  use mod_eigensystem,          only: nmatmax_ptr
  use mod_spin,                 only: nspinor
  use mod_eigenvalue_occupancy, only: nstfv
  use m_getunit
  use m_plotmat

  implicit none

! methods
  contains

    subroutine wannier_scdm
      use mod_rgrid
      use mod_lattice, only: omega

      type( rgrid) :: grid
      integer :: ngrid(3)
      real(8) :: boxl(4,3)

      integer :: ik, ir, irmt, ist, is, ia, ias, igk, la1, la2, la3, n1, n2, n3, i, j, np2, ir0, lm, fst, lst
      real(8) :: rho, r, t1, t2,  vl(3), vc(3)
      real(8) :: sigma(2), fe
      complex(8) :: zsum

      real(8), allocatable :: xa(:), ya(:), c(:), sval(:), eval(:)
      real(8), allocatable :: evalfv(:,:)
      complex(8), allocatable :: wfongrid(:,:), zeta(:,:), x(:,:), y(:,:), ur(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:), wfmt(:,:,:), wfir(:), zylm(:), lsvec(:,:), rsvec(:,:)

      integer :: lp_lwork, lp_info
      integer, allocatable :: lp_jpvt(:)
      real(8), allocatable :: lp_rwork(:)
      complex(8), allocatable :: lp_tau(:), lp_work(:)

      real(8), external :: polynom

      ik = 1
      rho = 500.d0

      la1 = norm2( input%structure%crystal%basevect( :, 1))
      la2 = norm2( input%structure%crystal%basevect( :, 2))
      la3 = norm2( input%structure%crystal%basevect( :, 3))
      n1 = max( 1, nint( (rho*la1*la1*omega/la2/la3)**(1.d0/3)))
      n2 = max( 1, nint( (rho*la2*la2*omega/la3/la1)**(1.d0/3)))
      n3 = max( 1, nint( (rho*la3*la3*omega/la1/la2)**(1.d0/3)))
      ngrid = (/n1, n2, n3/)
      boxl( 1, :) = (/0.d0, 0.d0, 0.d0/)
      boxl( 2, :) = (/1.d0, 0.d0, 0.d0/)
      boxl( 3, :) = (/0.d0, 1.d0, 0.d0/)
      boxl( 4, :) = (/0.d0, 0.d0, 1.d0/)

      grid = gen_3d( ngrid, boxl, 1)

      rho = dsqrt( omega/grid%npt)
      rho = 1.d0
      write(*,*) ngrid
      write(*,*) wf_groups( wf_group)%win_i
      write(*,*) wf_groups( wf_group)%win_o
      write(*,*) wf_groups( wf_group)%fst, wf_groups( wf_group)%lst, wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf
      sigma = 1.d-10 + 0.1d0*abs( wf_groups( wf_group)%win_o - wf_groups( wf_group)%win_i)
      sigma = 1.d0/sigma
      write(*,*) sigma

      allocate( evecfv( nmatmax_ptr, nstfv, nspinor))
      allocate( apwalm( ngkmax_ptr, apwordmax, lmmaxapw, natmtot))
      allocate( wfmt( lmmaxapw, nrmtmax, natmtot))
      allocate( wfir( ngrtot))
      allocate( wfongrid( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, grid%npt))

      apwalm = zzero
      call wannier_getevec( ik, evecfv)
      call wannier_geteval( evalfv, fst, lst)
      call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)

      do ist = wf_groups( wf_group)%fst, wf_groups( wf_group)%lst
        if( norm2( wf_groups( wf_group)%win_o) .gt. 1.d-20) then
          if( (evalfv( ist, ik) .lt. wf_groups( wf_group)%win_o(1)) .or. (evalfv( ist, ik) .gt. wf_groups( wf_group)%win_o(2))) then
            fe = 0.d0
          else
            fe = ( exp( sigma(1)*( wf_groups( wf_group)%win_i(1) - evalfv( ist, ik))) + 1.d0)*( exp( sigma(2)*( evalfv( ist, ik) - wf_groups( wf_group)%win_i(2))) + 1.d0)
            fe = 1.d0/fe
          end if
        else
          fe = 1.d0
        end if
        if( fe .gt. 1.d-6) then
          wfmt = zzero
          wfir = zzero
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              call wavefmt( 1, input%groundstate%lmaxapw, is, ia, wf_Gkset%ngk( 1, ik), apwalm, evecfv( :, ist, 1), lmmaxapw, wfmt( :, :, ias))
            end do
          end do
          do igk = 1, wf_Gkset%ngk( 1, ik)
            wfir( igk) = evecfv( igk, ist, 1)/dsqrt( omega)
          end do
          call calc_zdata_rgrid( grid, ik, wfmt, wfir, wfongrid( ist, :), nosym=.true.)
          write(*,*) ist, fe, evalfv( ist, ik)
          wfongrid( ist, :) = conjg( wfongrid( ist, :))*rho*fe
        end if
      end do
      !call writematlab( transpose( wfongrid), 'psi')

      allocate( lp_jpvt( grid%npt))
      allocate( lp_tau( wf_nst))
      allocate( lp_work( 1))
      allocate( lp_rwork( 2*grid%npt))
      lp_jpvt = 0
      call zgeqp3( wf_nst, grid%npt, wfongrid, wf_nst, lp_jpvt, lp_tau, lp_work, -1, lp_rwork, lp_info)
      lp_lwork = lp_work( 1)
      deallocate( lp_work)
      allocate( lp_work( lp_lwork))
      call zgeqp3( wf_groups( wf_group)%nst, grid%npt, wfongrid, wf_groups( wf_group)%nst, lp_jpvt, lp_tau, lp_work, lp_lwork, lp_rwork, lp_info)
      write(*,*) lp_info
      do ir = 1, wf_groups( wf_group)%nwf
        write(*,'(i,3f13.6)') lp_jpvt( ir), grid%vpl( :, lp_jpvt( ir))
      end do
      deallocate( lp_tau, lp_work, lp_rwork)

      np2 = input%groundstate%nprad/2
      allocate( zeta( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%nwf))
      allocate( sval( wf_groups( wf_group)%nwf), &
                lsvec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst), &
                rsvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf), &
                eval( wf_groups( wf_group)%nst))
      allocate( zylm( lmmaxapw))
      allocate( xa( input%groundstate%nprad))
      allocate( ya( input%groundstate%nprad))
      allocate( c( input%groundstate%nprad))
      allocate( x( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))
      allocate( y( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      allocate( ur( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      
      do ik = 1, wf_kset%nkpt
        call wannier_getevec( ik, evecfv)
        call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)
        do ist = wf_groups( wf_group)%fst, wf_groups( wf_group)%lst
          zeta( ist, :) = zzero
          if( norm2(wf_groups( wf_group)%win_o) .gt. 1.d-20) then
            if( (evalfv( ist, ik) .lt. wf_groups( wf_group)%win_o(1)) .or. (evalfv( ist, ik) .gt. wf_groups( wf_group)%win_o(2))) then
              fe = 0.d0
            else
              fe = ( exp( sigma(1)*( wf_groups( wf_group)%win_i(1) - evalfv( ist, ik))) + 1.d0)*( exp( sigma(2)*( evalfv( ist, ik) - wf_groups( wf_group)%win_i(2))) + 1.d0)
              fe = 1.d0/fe
            end if
          else
            fe = 1.d0
          end if
          if( fe .gt. 1.d-6) then
            do ir = 1, wf_groups( wf_group)%nwf
              vl(:) = grid%vpl( :, lp_jpvt( ir))
              vc(:) = grid%vpc( :, lp_jpvt( ir))
              if( grid%mtpoint( lp_jpvt( ir))) then
                wfmt = zzero
                is = grid%atom( 1, lp_jpvt( ir))
                ia = grid%atom( 2, lp_jpvt( ir))
                ias = idxas( ia, is)

                call wavefmt( 1, input%groundstate%lmaxapw, is, ia, wf_Gkset%ngk( 1, ik), apwalm, evecfv( :, ist, 1), lmmaxapw, wfmt( :, :, ias))
                call ylm( vc, input%groundstate%lmaxapw, zylm)
                r = norm2( vc)
                do irmt = 1, nrmt( is)
                  if( spr( irmt, is) .ge. r) then
                    if( irmt .le. np2) then
                      ir0 = 1
                    else if( irmt .gt. nrmt( is)-np2) then
                      ir0 = nrmt( is) - input%groundstate%nprad + 1
                    else
                      ir0 = irmt - np2
                    end if
                    r = max( r, spr( 1, is))
                    zsum = zzero
                    do lm = 1, lmmaxapw
                      do j = 1, input%groundstate%nprad
                        i = ir0 + j - 1
                        xa( j) = dble( wfmt( lm, i, ias))
                        ya( j) = aimag( wfmt( lm, i, ias))
                      end do
                      ! interpolate radial part of wfmt
                      t1 = polynom( 0, input%groundstate%nprad, spr( ir0, is), xa, c, r)
                      t2 = polynom( 0, input%groundstate%nprad, spr( ir0, is), ya, c, r)
                      zsum = zsum + cmplx( t1, t2, 8)*zylm( lm)
                    end do ! lm
                    exit ! the loop over ir
                  end if
                end do ! ir

                t1 = twopi*dot_product( wf_kset%vkl( :, ik), dble( grid%iv( :, lp_jpvt( ir))))
                zeta( ist, ir) = conjg( zsum*cmplx( cos( t1), sin( t1), 8))
              else
                zsum = 0.d0
                do igk = 1, wf_Gkset%ngk( 1, ik)
                  t1 = twopi*dot_product( wf_Gkset%vgkl( :, igk, 1, ik), vl)
                  zsum = zsum + evecfv( igk, ist, 1)/dsqrt( omega)*cmplx( cos( t1), sin( t1), 8)
                end do
                zeta( ist, ir) = conjg( zsum)
              end if
            end do
            zeta( ist, :) = fe*zeta( ist, :)
          end if
        end do
        call zsvd( zeta, sval, lsvec, rsvec)
        ! for numerical stability
        !do i = 1, wf_nwf
        !  if( sval( i) .lt. 1.d-12) lsvec( :, i) = zzero
        !end do
        !write(*,'(100f13.6)') sval
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               lsvec, wf_groups( wf_group)%nst, &
               rsvec, wf_groups( wf_group)%nwf, zzero, &
               wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik), wf_groups( wf_group)%nst)
        if( norm2( wf_groups( wf_group)%win_o) .gt. 1.d-20) then
          ur = zzero
          !write(*,*) ik, wf_win_no( ik), wf_win_ni( ik)
          !write(*,*) wf_win_io( 1:wf_win_no( ik), ik)
          !write(*,*) wf_win_ii( 1:wf_win_ni( ik), ik)
          do ist = 1, wf_groups( wf_group)%win_no( ik)
            ur( ist, :) = wf_transform( wf_groups( wf_group)%win_io( ist, ik), wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik)
          end do
          call zgemm( 'n', 'c', wf_groups( wf_group)%win_no( ik), wf_groups( wf_group)%win_no( ik), wf_groups( wf_group)%nwf, zone, &
                 ur( 1:wf_groups( wf_group)%win_no( ik), :), wf_groups( wf_group)%win_no( ik), &
                 ur( 1:wf_groups( wf_group)%win_no( ik), :), wf_groups( wf_group)%win_no( ik), zzero, &
                 y( 1:wf_groups( wf_group)%win_no( ik), 1:wf_groups( wf_group)%win_no( ik)), wf_groups( wf_group)%win_no( ik))
          call zhediag( y( 1:wf_groups( wf_group)%win_no( ik), 1:wf_groups( wf_group)%win_no( ik)), &
                  eval( 1:wf_groups( wf_group)%win_no( ik)), lsvec( 1:wf_groups( wf_group)%win_no( ik), 1:wf_groups( wf_group)%win_no( ik)))
          y = zzero
          do ist = 1, wf_groups( wf_group)%nwf - wf_groups( wf_group)%win_ni( ik)
            y( 1:wf_groups( wf_group)%win_no( ik), ist) = lsvec( 1:wf_groups( wf_group)%win_no( ik), wf_groups( wf_group)%win_no( ik) - ist + 1)
          end do
          !write(*,*) "Y"
          !call plotmat( y(1:wf_groups( wf_group)%win_no( ik), 1:(wf_groups( wf_group)%nwf-wf_groups( wf_group)%win_ni( ik))))
          ur = zzero
          do ist = 1, wf_groups( wf_group)%win_ni( ik)
            ur( wf_groups( wf_group)%win_ii( ist, ik) - wf_groups( wf_group)%fst + 1, ist) = zone
          end do
          do ist = 1, wf_groups( wf_group)%win_no( ik)
            ur( wf_groups( wf_group)%win_io( ist, ik) - wf_groups( wf_group)%fst + 1, (wf_groups( wf_group)%win_ni( ik)+1):wf_groups( wf_group)%nwf) = y( ist, :)
          end do
          !write(*,*) "BLA"
          !call plotmat( ur)
          call zgemm( 'c', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
                 ur, wf_groups( wf_group)%nst, &
                 wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik), wf_groups( wf_group)%nst, zzero, &
                 x, wf_groups( wf_group)%nwf)
          call zsvd( x, sval, lsvec( 1:wf_groups( wf_group)%nwf, 1:wf_groups( wf_group)%nwf), rsvec)
          call zgemm( 'n', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
                 lsvec( 1:wf_groups( wf_group)%nwf, 1:wf_groups( wf_group)%nwf), wf_groups( wf_group)%nwf, &
                 rsvec, wf_groups( wf_group)%nwf, zzero, &
                 x, wf_groups( wf_group)%nwf)
          call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
                 ur, wf_groups( wf_group)%nst, &
                 x, wf_groups( wf_group)%nwf, zzero, &
                 wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik), wf_groups( wf_group)%nst)
          !write(*,*) "----------------------"
        end if
      end do

      deallocate( evecfv, apwalm, wfmt, wfir, wfongrid, zylm, lp_jpvt, zeta, sval, lsvec, rsvec, xa, ya, c)

      call wannier_loc
      return
    end subroutine wannier_scdm

end module mod_wannier_scdm
