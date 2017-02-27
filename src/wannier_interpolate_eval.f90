! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bandstr
! !INTERFACE:
!
!
subroutine wannier_interpolate_eval( eval1, nk1, kvl1, eval2, nk2, kvl2, fst, lst)
  ! !USES:
  use mod_wannier
  use m_wsweight

  ! !DESCRIPTION:
  !
  ! !REVISION HISTORY:
  !EOP
  !BOC
  implicit none
  integer, intent( in) :: nk1, nk2, fst, lst
  real(8), intent( in) :: eval1( fst:lst, nk1)
  real(8), intent( in) :: kvl1( 3, nk1), kvl2( 3, nk2)
  real(8), intent( out) :: eval2( fst:lst, nk2)
  
  integer :: nrpt, ia, ix, iy, iz, iknr, ik, ir
  real(8) :: vc(3), energy
  complex(8) :: ftweight

  real(8), allocatable :: rptc(:,:), rptl(:,:), evaltmp(:), kdiff(:,:), eval1_tmp(:,:)
  complex(8), allocatable :: auxmat(:,:), auxmat2(:,:,:), wanme(:,:,:), evectmp(:,:)

  ! check band-range
  if( (fst .lt. wf_fst) .or. (lst .gt. wf_lst)) then
    write( *, '(" ERROR (wannier_interpolate_eval): The selected band-range for interpolation (",I3,"-",I3,") does not lie in the band-range for Wannier-functions (",I3,"-",I3,").")') fst, lst, wf_fst, wf_lst
    call terminate
  end if
  ! check k-points
  if( nk1 .ne. wf_kset%nkpt) then
    write( *, '(" ERROR (wannier_interpolate_eval): Different numbers of k-points.")')
    call terminate
  end if
  
  allocate( eval1_tmp( fst:lst, nk1))
  eval1_tmp( lst, :) = eval1( lst, :)
  eval1_tmp( lst-1, :) = eval1( lst, :)
  eval1_tmp( fst:(lst-2), :) = eval1( (fst+2):lst, :)

  ! generate set of lattice vectors 
  nrpt = nk1
  allocate( rptc( 3, nrpt), rptl( 3, nrpt))
  ia = 0
  do iz = -wf_kset%ngridk(3)/2, -wf_kset%ngridk(3)/2+wf_kset%ngridk(3)-1
    do iy = -wf_kset%ngridk(2)/2, -wf_kset%ngridk(2)/2+wf_kset%ngridk(2)-1
      do ix = -wf_kset%ngridk(1)/2, -wf_kset%ngridk(1)/2+wf_kset%ngridk(1)-1
        ia = ia + 1
        rptl( :, ia) = (/ dble( ix), dble( iy), dble( iz)/)
        call r3mv( input%structure%crystal%basevect, rptl( :, ia), rptc( :, ia))
      end do
    end do
  end do
  
  ! calculate Hamlitonian matrix elements in Wannier representation 
  allocate( auxmat( nrpt, nk1), auxmat2( nk1, fst:lst, fst:lst))
  allocate( wanme( nrpt, fst:lst, fst:lst))
  allocate( kdiff( 3, wf_kset%nkpt))
  do iknr = 1, nk1
    kdiff( 1, :) = wf_kset%vkl( 1, :) - kvl1( 1, iknr)
    kdiff( 2, :) = wf_kset%vkl( 2, :) - kvl1( 2, iknr)
    kdiff( 3, :) = wf_kset%vkl( 3, :) - kvl1( 3, iknr)
    ik = minloc( norm2( kdiff, 1), 1)
    if( norm2( kdiff( :, ik)) .gt. input%structure%epslat) then
      write( *, '(" ERROR (wannier_interpolate_eval): k-point does not belong to Wannier k-grid.")')
      write( *, '(" k-point: ",3F13.6)') kvl1( :, iknr)
      call terminate
    end if
    !write( *, '(3F13.6)') kvl1( :, iknr)
    do iy = fst, lst
      !write( *, '(I3,F23.16)') iy, eval1_tmp( iy, iknr)
      do ix = fst, lst
        auxmat2( iknr, ix, iy) = zzero
        do iz = fst, lst
          auxmat2( iknr, ix, iy) = auxmat2( iknr, ix, iy) + eval1( iz, iknr)*wf_transform( iz-fst+wf_fst, ix-fst+wf_fst, ik)*conjg( wf_transform( iz-fst+wf_fst, iy-fst+wf_fst, ik))
        end do
      end do
    end do
    call r3mv( bvec, wf_kset%vkl( :, ik), vc)
    do ir = 1, nrpt 
      auxmat( ir, iknr) = exp( -zi*dot_product( rptc( :, ir), vc))
    end do
  end do
  do iy = fst, lst
    call ZGEMM( 'N', 'N', nrpt, lst-fst+1, nk1, zone/nk1, &
         auxmat, nrpt, &
         auxmat2( :, :, iy), nk1, zzero, &
         wanme( :, :, iy), nrpt)
  end do
  deallocate( kdiff, auxmat2, rptc, auxmat)

  ! interpolation

#ifndef USEOMP
  allocate( auxmat( fst:lst, fst:lst), evectmp( fst:lst, fst:lst))
  allocate( evaltmp( fst:lst))
#endif
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, ir, auxmat, evaltmp, evectmp, ftweight, ix)
  allocate( auxmat( fst:lst, fst:lst), evectmp( fst:lst, fst:lst))
  allocate( evaltmp( fst:lst))
#endif
#ifdef USEOMP
!$OMP DO
#endif
  do ik = 1, nk2 
    auxmat(:,:) = zzero
    do ir = 1, nrpt
      call ws_weight( rptl( :, ir), rptl( :, ir), kvl2( :, ik), ftweight, .true.)
      auxmat(:,:) = auxmat(:,:) + wanme( ir, :, :)*conjg( ftweight)
    end do
    call diaghermat( lst-fst+1, auxmat, evaltmp, evectmp)
    do ix = fst, lst
      eval2( ix, ik) = evaltmp( ix)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( auxmat, evaltmp, evectmp)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
  deallocate( wanme, rptl)
  return
end subroutine wannier_interpolate_eval
!EOC
