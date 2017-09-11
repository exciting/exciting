!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine wannier_interpolate_bandchar( lmax, is, ia, int_kset, evecint, phase, bc)
  use mod_wannier
  use modinput
  use modmain

  implicit none

  integer, intent( in) :: lmax, is, ia
  type( k_set), intent( in) :: int_kset
  complex(8), intent( in) :: evecint( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt)
  complex(8), intent( in) :: phase( wf_kset%nkpt, int_kset%nkpt)
  real(4), intent( out) :: bc( 0:lmax, wf_fst:wf_lst, int_kset%nkpt)

  integer :: lmmax
  integer :: l, m, lm, i, j, irc, ik, iq
  real(8) :: fr( nrcmtmax), gr( nrcmtmax), cf( 3, nrcmtmax)

  complex(8), allocatable :: wfmt_tmp(:,:), wfmt(:,:,:,:), evecfv(:,:), apwalm(:,:,:,:), wfmt_int(:,:)
  lmmax = (lmax+1)**2

  call readstate
  call readfermi
  call linengy
  call genapwfr
  call genlofr

  allocate( wfmt( lmmax, wf_fst:wf_lst, nrcmtmax, wf_kset%nkpt))

  wfmt = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, evecfv, apwalm, i, j, wfmt_tmp)
#endif
  allocate( wfmt_tmp( lmmax, nrcmtmax))
  allocate( evecfv( nmatmax, nstfv))
  allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
#ifdef USEOMP
!$OMP DO
#endif
  do ik = 1, wf_kset%nkpt
    if( input%properties%wannier%input .eq. "groundstate") then
      call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
    else if( input%properties%wannier%input .eq. "hybrid") then
      call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
    else if( input%properties%wannier%input .eq. "gw") then
      call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspinor, evecfv)
    else
      call terminate
    end if
    call match( wf_Gkset%ngk( 1, ik), &
                wf_Gkset%gkc( :, 1, ik), &
                wf_Gkset%tpgkc( :, :, 1, ik), &
                wf_Gkset%sfacgk( :, :, 1, ik), &
                apwalm( :, :, :, :))
    do i = wf_fst, wf_lst
      call wavefmt( input%groundstate%lradstep, lmax, is, ia, wf_Gkset%ngk( 1, ik), apwalm, evecfv(:, i), lmmax, wfmt_tmp)
      do j = wf_fst, wf_lst
        wfmt( :, j, :, ik) = wfmt( :, j, :, ik) + wf_transform( i, j, ik)*wfmt_tmp(:,:)
      end do
    end do
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( apwalm, evecfv, wfmt_tmp)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
  
  bc = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, l, m, lm, ik, wfmt_int, i, irc, fr, gr, cf)
#endif
  allocate( wfmt_int( wf_fst:wf_lst, nrcmtmax))
#ifdef USEOMP
!$OMP DO
#endif
  do iq = 1, int_kset%nkpt
    do l = 0, lmax
      do m = -l, l
        lm = idxlm( l, m)
        wfmt_int = zzero
        do ik = 1, wf_kset%nkpt
          call zgemm( 'T', 'N', wf_nst, nrcmtmax, wf_nst, phase( ik, iq), &
               evecint( :, :, iq), wf_nst, &
               wfmt( lm, :, :, ik), wf_nst, zone, &
               wfmt_int, wf_nst)
        end do
        do i = wf_fst, wf_lst
          do irc = 1, nrcmt( is)
            fr( irc) = dble( conjg( wfmt_int( i, irc))*wfmt_int( i, irc))*rcmt( irc, is)**2
          end do
          call fderiv( -1, nrcmt( is), rcmt( :, is), fr, gr, cf)
          !write(*,'(4I,F23.16)') iq, i, l, m, gr( nrcmt( is))
          bc( l, i, iq) = bc( l, i, iq) + gr( nrcmt( is))
        end do
      end do
    end do
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( wfmt_int)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
  deallocate( wfmt)
  return
end subroutine
