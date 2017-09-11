module m_wannier_interpolate_fft
    implicit none
    contains

subroutine wannier_interpolate_fft( evalin, int_kset, evalout, evecout)
  use modmain
  use mod_wannier
  use m_wsweight
  use m_plotmat
  use mod_kqpts
  use mod_lattice
  use mod_eigenvalue_occupancy
  use mod_eigensystem
  use mod_atoms
  use mod_muffin_tin
  use mod_lattice
  use mod_constants

  implicit none
  real(8), intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
  type( k_set), intent( in) :: int_kset
  real(8), intent( out) :: evalout( wf_fst:wf_lst, int_kset%nkpt)
  complex(8), intent( out) :: evecout( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt)
  
  integer :: nrpt, ix, iy, iz, ik, iq, ir

  real(8), allocatable :: rptl(:,:)
  complex(8), allocatable :: auxmat(:,:), ueu(:,:,:), hamilton(:,:,:)
  complex(8), allocatable :: hreal(:,:,:)
  complex(8), allocatable :: prk(:,:), prq(:,:)

  !**********************************************
  ! interpolated eigenenergies and corresponding 
  ! eigenvectors in Wannier basis
  !**********************************************

  !call generate_k_vectors( int_kset, bvec, fact*wf_kset%ngridk, wf_kset%vkloff, wf_kset%isreduced)
  ! generate set of lattice vectors 
  nrpt = wf_kset%nkpt
  allocate( rptl( 3, nrpt))
  ir = 0
  do ix = -wf_kset%ngridk(1)/2, -wf_kset%ngridk(1)/2+wf_kset%ngridk(1)-1
    do iy = -wf_kset%ngridk(2)/2, -wf_kset%ngridk(2)/2+wf_kset%ngridk(2)-1
      do iz = -wf_kset%ngridk(3)/2, -wf_kset%ngridk(3)/2+wf_kset%ngridk(3)-1
        ir = ir + 1
        rptl( :, ir) = (/dble( ix), dble( iy), dble( iz)/)
      end do
    end do
  end do

  ! calculate phases for Fourier transform on Wigner-Seitz supercell
  allocate( prk( nrpt, wf_kset%nkpt), prq( nrpt, int_kset%nkpt))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iq, ir)
!$OMP DO
#endif
  do ir = 1, nrpt
    do iq = 1, int_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), int_kset%vkl( :, iq), prq( ir, iq), kgrid=.true.)
    end do
    do ik = 1, wf_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), prk( ir, ik), kgrid=.true.)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  
  ! calculate Hamlitonian in Wannier gauge 
  allocate( ueu( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iy, ir, auxmat)
#endif
  allocate( auxmat( wf_fst:wf_lst, wf_fst:wf_lst))
#ifdef USEOMP
!$OMP DO
#endif
  do ik = 1, wf_kset%nkpt
    auxmat = zzero
    do iy = wf_fst, wf_lst
      auxmat( iy, :) = wf_transform( iy, :, ik)*evalin( iy, ik)
    end do
    call zgemm( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
         wf_transform( :, :, ik), wf_nst, &
         auxmat, wf_nst, zzero, &
         ueu( :, :, ik), wf_nst)
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( auxmat)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

  ! calculate Hamiltonian in real space
  allocate( hreal( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))
  hreal = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, ir)
!$OMP DO
#endif
  do ir = 1, nrpt
    do ik = 1, wf_kset%nkpt
      hreal( :, :, ir) = hreal( :, :, ir) + conjg( prk( ir, ik))*ueu( :, :, ik)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  deallocate( ueu)
  hreal = hreal/nrpt

  ! calculate interpolated Hamiltonian in reciprocal space
  allocate( hamilton( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))
  hamilton = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ir)
!$OMP DO
#endif
  do iq = 1, int_kset%nkpt
    do ir = 1, nrpt
      hamilton( :, :, iq) = hamilton( :, :, iq) + prq( ir, iq)*hreal( :, :, ir)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

  ! diagonalize interpolated Hamiltonian
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq)
!$OMP DO
#endif
  do iq = 1, int_kset%nkpt 
    call diaghermat( wf_nst, hamilton( :, :, iq), evalout( :, iq), evecout( :, :, iq))
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  deallocate( rptl, hreal, hamilton, prk, prq)

  return
end subroutine wannier_interpolate_fft

end module m_wannier_interpolate_fft
