module m_wannier_interpolate
    implicit none
    contains

subroutine wannier_interpolate_eval( eval1, nk2, kvl2, evalint)!, bandchar, lmax)
  use mod_wannier
  use m_wsweight
  use m_plotmat
  use mod_eigenvalue_occupancy

  implicit none
  integer, intent( in) :: nk2
  real(8), intent( in) :: eval1( wf_fst:wf_lst, wf_kset%nkpt)
  real(8), intent( in) :: kvl2( 3, nk2)
  real(8), intent( out) :: evalint( wf_fst:wf_lst, nk2)
  !real(4), intent( out), optional :: bandchar( natmtot, 0:lmax, wf_fst:wf_lst, nk2)
  !integer, intent( in), optional :: lmax
  
  integer :: nrpt, ia, ix, iy, iz, iknr, iknr2, ik, ir
  complex(8) :: ftweight
  real(8) :: t1, v1(3), v2(3)

  real(8), allocatable :: rptl(:,:), evaltmp(:)
  complex(8), allocatable :: auxmat(:,:), ueu(:,:,:), phase(:,:), hamilton(:,:,:), evecint(:,:,:)
  
  integer :: io1, io2, is, ias, l, m, lm, ngknr 
  real(8), allocatable :: apwint(:,:,:,:), loint(:,:), f1(:), f2(:), gf(:), cf(:,:)
  complex(8) :: zt
  integer, allocatable :: igkignr(:), cnt(:)
  real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:), sval(:)
  complex(8), allocatable :: apwalm(:,:,:,:,:), sfacgknr(:,:), evecfv(:,:,:), wfmt(:,:,:,:), rint(:,:), purint(:,:,:)
  complex(8), allocatable :: p1(:,:), p2(:,:)

  ! generate set of lattice vectors 
  nrpt = wf_kset%nkpt
  allocate( rptl( 3, nrpt))
  ia = 0
  do iz = -wf_kset%ngridk(3)/2, -wf_kset%ngridk(3)/2+wf_kset%ngridk(3)-1
    do iy = -wf_kset%ngridk(2)/2, -wf_kset%ngridk(2)/2+wf_kset%ngridk(2)-1
      do ix = -wf_kset%ngridk(1)/2, -wf_kset%ngridk(1)/2+wf_kset%ngridk(1)-1
        ia = ia + 1
        rptl( :, ia) = (/ dble( ix), dble( iy), dble( iz)/)
      end do
    end do
  end do
  
  ! calculate Hamlitonian matrix elements in Wannier representation 
  allocate( ueu( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))
  allocate( phase( wf_kset%nkpt, nk2))
  allocate( hamilton( wf_fst:wf_lst, wf_fst:wf_lst, nk2))
  phase = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ik, iy, auxmat)
#endif
  allocate( auxmat( wf_fst:wf_lst, wf_fst:wf_lst))
#ifdef USEOMP
!$OMP DO
#endif
  do iknr = 1, wf_kset%nkpt
    do iy = wf_fst, wf_lst
      ueu( iy, :, iknr) = wf_transform( iy, :, iknr)*eval1( iy, iknr)
    end do
    call zgemm( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
         wf_transform( :, :, iknr), wf_nst, &
         ueu( :, :, iknr), wf_nst, zzero, &
         auxmat, wf_nst)
    ueu( :, :, iknr) = auxmat
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( auxmat)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

  allocate( p1( nrpt, wf_kset%nkpt), p2( nrpt, nk2))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ik, ir, t1, v1, v2)
!$OMP DO
#endif
  do ir = 1, nrpt
    !call r3mv( input%structure%crystal%basevect, rptl( :, ir), v1)
    do ik = 1, nk2
      call ws_weight( rptl( :, ir), rptl( :, ir), kvl2( :, ik), p2( ir, ik), kgrid=.true.)
      !call r3mv( bvec, kvl2( :, ik), v2)
      !t1 = dot_product( v1, v2)
      !p2( ik, ir) = cmplx( cos( t1), sin( t1), 8)
    end do
    do iknr = 1, wf_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, iknr), p1( ir, iknr), kgrid=.true.)
      !t1 = dot_product( v1, wf_kset%vkc( :, iknr))
      !p1( iknr, ir) = cmplx( cos( t1), sin( t1), 8)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  call zgemm( 'C', 'N', wf_kset%nkpt, nk2, nrpt, zone, &
       p1, nrpt, &
       p2, nrpt, zzero, &
       phase, wf_kset%nkpt)

  phase = phase/wf_kset%nkpt
  do iy = wf_fst, wf_lst
    call zgemm( 'N', 'N', wf_nst, nk2, wf_kset%nkpt, zone, &
         ueu( iy, :, :), wf_nst, &
         phase, wf_kset%nkpt, zzero, &
         hamilton( iy, :, :), wf_nst)
  end do
  deallocate( ueu, p1, p2)

  ! interpolation
  allocate( evecint( wf_fst:wf_lst, wf_fst:wf_lst, nk2))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik)
!$OMP DO
#endif
  do ik = 1, nk2 
    call diaghermat( wf_nst, hamilton( :, :, ik), evalint( :, ik), evecint( :, :, ik))
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  deallocate( rptl, hamilton)

!  if( (.not. present( bandchar)) .or. (.not. present( lmax)) .or. (lmax .lt. 0)) then
!    write(*,*) "not present"
!    return
!  end if
!  
!  write(*,*) nstfv, nstsv
!
!  call readstate
!  call linengy
!  call genapwfr
!  call genlofr
!
!  allocate( wfmt( (lmax+1)**2, nrcmtmax, wf_fst:wf_lst, max( wf_kset%nkpt, nk2)))
!  allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
!  allocate( igkignr( ngkmax))
!  allocate( vgklnr( 3, ngkmax, nspinor), vgkcnr( 3, ngkmax, nspinor), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
!  allocate( sfacgknr( ngkmax, natmtot))
!  allocate( evecfv( nmatmax, nstfv, nspinor))
!
!  do is = 1, nspecies
!    do ia = 1, natoms( is)
!      ias = idxas( ia, is)
!      write(*,*) ias
!      allocate( auxmat( nrcmtmax, wf_fst:wf_lst))
!      do iknr = 1, wf_kset%nkpt
!        call gengpvec( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
!        call gensfacgp( ngknr, vgkcnr, ngkmax, sfacgknr)
!        call match( ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm(:, :, :, :, 1))
!        if( input%properties%wannier%input .eq. "groundstate") then
!          call getevecfv( wf_kset%vkl( :, iknr), vgklnr, evecfv)
!        else if( input%properties%wannier%input .eq. "hybrid") then
!          call getevecfv( wf_kset%vkl( :, iknr), vgklnr, evecfv)
!        else if( input%properties%wannier%input .eq. "gw") then
!          call getevecsvgw_new( "GW_EVECSV.OUT", iknr, wf_kset%vkl( :, iknr), nmatmax, nstfv, nspinor, evecfv)
!        else
!          call terminate
!        end if
!        do ix = wf_fst, wf_lst
!          call wavefmt( input%groundstate%lradstep, lmax, is, ia, ngknr, apwalm(:, :, :, :, 1), evecfv( :, ix, 1), (lmax+1)**2, wfmt( :, :, ix, iknr))
!        end do
!        do l = 0, lmax
!          do m = -l, l
!            lm = idxlm( l, m)
!            call zgemm( 'N', 'N', nrcmtmax, wf_nst, wf_nst, zone, &
!                 wfmt( lm, :, :, iknr), nrcmtmax, &
!                 wf_transform( :, :, iknr), wf_nst, zzero, &
!                 auxmat, nrcmtmax)
!            wfmt( lm, :, :, iknr) = auxmat(:,:)
!          end do
!        end do
!      end do
!      deallocate( auxmat)
!      allocate( auxmat( nrcmtmax, nk2))
!      do ix = wf_fst, wf_lst
!        do l = 0, lmax
!          do m = -l, l
!            lm = idxlm( l, m)
!            call zgemm( 'N', 'N', nrcmtmax, nk2, wf_kset%nkpt, zone, &
!                 wfmt( lm, :, ix, :), nrcmtmax, &
!                 phase, wf_kset%nkpt, zzero, &
!                 auxmat, nrcmtmax)
!            wfmt( lm, :, ix, 1:nk2) = auxmat(:,:)
!          end do
!        end do
!      end do
!      deallocate( auxmat)
!      allocate( auxmat( nrcmtmax, wf_fst:wf_lst))
!      do ik = 1, nk2
!        do l = 0, lmax
!          do m = -l, l
!            lm = idxlm( l, m)
!            call zgemm( 'N', 'N', nrcmtmax, wf_nst, wf_nst, zone, &
!                 wfmt( lm, :, :, ik), nrcmtmax, &
!                 evecint( :, :, ik), wf_nst, zzero, &
!                 auxmat, nrcmtmax)
!            wfmt( lm, :, :, ik) = auxmat(:,:)
!          end do
!        end do
!      end do
!      deallocate( auxmat)
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iz, l, m, lm, ir, f1, gf, cf)
!#endif
!      allocate( f1( nrcmtmax), gf( nrcmtmax), cf( 3, nrcmtmax))
!#ifdef USEOMP
!!$OMP DO
!#endif
!      do ik = 1, nk2
!        do iz = wf_fst, wf_lst
!          do l = 0, lmax
!            bandchar( ias, l, iz, ik) = 0.d0
!            do m = -l, l
!              lm = idxlm( l, m)
!              do ir = 1, nrcmt( is)
!                f1( ir) = dble( wfmt( lm, ir, iz, ik)*conjg( wfmt( lm, ir, iz, ik)))*rcmt( ir, is)**2
!              end do
!              call fderiv( -1, nrcmt( is), rcmt( :, is), f1, gf, cf)
!              bandchar( ias, l, iz, ik) = bandchar( ias, l, iz, ik) + gf( nrcmt( is))
!            end do
!          end do
!        end do
!      end do
!#ifdef USEOMP
!!$OMP END DO
!#endif
!      deallocate( f1, gf, cf)
!#ifdef USEOMP
!!$OMP END PARALLEL
!#endif
!    end do
!  end do
!
!  deallocate( wfmt, apwalm, igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr, sfacgknr, evecfv)
    
  return
end subroutine wannier_interpolate_eval
end module m_wannier_interpolate
