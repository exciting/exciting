module mod_wfint
  use modmain
  use mod_wannier
  use m_plotmat
  implicit none

! module variables
  type( k_set) :: wfint_kset                        ! k-point set on which the interpolation is performed
  type( Gk_set) :: wfint_Gkset                      ! G+k-point set on which the interpolation is performed
  logical :: wfint_initialized = .false.
  real(8) :: wfint_efermi                           ! interpolated fermi energy

  real(8), allocatable :: wfint_eval(:,:)           ! interpolated eigenenergies
  complex(8), allocatable :: wfint_evec(:,:,:)      ! interpolated eigenvectors
  complex(8), allocatable :: wfint_transform(:,:,:) ! corresponding expansion coefficients
  real(8), allocatable :: wfint_occ(:,:)            ! interpolated occupation numbers
  real(8), allocatable :: wfint_phase(:,:)          ! summed phase factors in interpolation
  complex(8), allocatable :: wfint_pkr(:,:)
  complex(8), allocatable :: wfint_pqr(:,:)

! methods
  contains
    !BOP
    ! !ROUTINE: wfint_init
    ! !INTERFACE:
    !
    subroutine wfint_init( int_kset, evalin_)
      ! !USES:
      ! !INPUT PARAMETERS:
      !   int_kset : k-point set on which the interpolation is performed on (in, type k_set)
      ! !DESCRIPTION:
      !   Sets up the interpolation grid and calculates the interpolated eigenenergies as well as 
      !   the corresponding expansion coefficients and phasefactors for the interpolated wavefunctions.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC
        type( k_set), intent( in) :: int_kset
        real(8), optional, intent( in) :: evalin_( wf_fst:wf_lst, wf_kset%nkpt)
    
        integer :: ik
        real(8), allocatable :: evalfv(:,:), evalin(:,:)

        if( wfint_initialized) call wfint_destroy

        wfint_kset = int_kset
    
        allocate( evalfv( nstfv, nspnfv))
        allocate( evalin( wf_fst:wf_lst, wf_kset%nkpt))
    
        if( present( evalin_)) then
          evalin = evalin_
          !do ik = 1, wf_kset%nkpt
          !  write(*,'(I,1000F13.6)') ik, evalin( :, ik)
          !end do
        else
          do ik = 1, wf_kset%nkpt
            call getevalfv( wf_kset%vkl( :, ik), evalfv)
            evalin( :, ik) = evalfv( wf_fst:wf_lst, 1)
          end do
        end if
    
        allocate( wfint_phase( wf_kset%nkpt, wfint_kset%nkpt))
        allocate( wfint_pkr( wf_kset%nkpt, wf_kset%nkpt))
        allocate( wfint_pqr( wfint_kset%nkpt, wf_kset%nkpt))
        allocate( wfint_eval( wf_fst:wf_lst, wfint_kset%nkpt))
        allocate( wfint_transform( wf_fst:wf_lst, wf_fst:wf_lst, wfint_kset%nkpt))
        call wfint_interpolate_eigsys( evalin)
    
        wfint_initialized = .true.
        deallocate( evalfv, evalin)

        return
    end subroutine wfint_init
    !EOC

!--------------------------------------------------------------------------------------
    
    !BOP
    ! !ROUTINE: wfint_init
    ! !INTERFACE:
    !
    subroutine wfint_destroy
      ! !USES:
      ! !INPUT PARAMETERS:
      !   int_kset : k-point set on which the interpolation is performed on (in, type k_set)
      ! !DESCRIPTION:
      !   Sets up the interpolation grid and calculates the interpolated eigenenergies as well as 
      !   the corresponding expansion coefficients and phasefactors for the interpolated wavefunctions.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC
    
        if( allocated( wfint_phase)) deallocate( wfint_phase)
        if( allocated( wfint_pkr)) deallocate( wfint_pkr)
        if( allocated( wfint_pqr)) deallocate( wfint_pqr)
        if( allocated( wfint_eval)) deallocate( wfint_eval)
        if( allocated( wfint_transform)) deallocate( wfint_transform)
        if( allocated( wfint_occ)) deallocate( wfint_occ)
        wfint_initialized = .false.
        return
    end subroutine wfint_destroy
    !EOC

!--------------------------------------------------------------------------------------
    
    !BOP
    ! !ROUTINE: wfint_interpolate_eigsys
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_eigsys( evalin)
      ! !USES:
      use m_wsweight
      ! !INPUT PARAMETERS:
      !   evalin : eigenenergies on original grid (in, real( wf_fst:wf_lst, wf_kset%nkpt))
      ! !DESCRIPTION:
      !   Calculates the interpolated eigenenergies as well as the corresponding expansion 
      !   coefficients and phasefactors for the interpolated wavefunctions.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC

      !use mod_kqpts
      !use mod_lattice
      !use mod_eigenvalue_occupancy
      !use mod_eigensystem
      !use mod_atoms
      !use mod_muffin_tin
      !use mod_lattice
      !use mod_constants
    
      implicit none
      real(8), intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
      
      integer :: nrpt, ix, iy, iz, ik, iq, ir
      complex(8) :: ftweight, z1
    
      real(8), allocatable :: rptl(:,:)
      complex(8), allocatable :: auxmat(:,:), ueu(:,:,:), hamilton(:,:,:)

      !**********************************************
      ! interpolated eigenenergies and corresponding 
      ! eigenvectors in Wannier basis
      !**********************************************
    
      ! generate set of lattice vectors 
      nrpt = wf_kset%nkpt
      allocate( rptl( 3, nrpt))
      ir = 0
      do iz = -wf_kset%ngridk(3)/2, -wf_kset%ngridk(3)/2+wf_kset%ngridk(3)-1
        do iy = -wf_kset%ngridk(2)/2, -wf_kset%ngridk(2)/2+wf_kset%ngridk(2)-1
          do ix = -wf_kset%ngridk(1)/2, -wf_kset%ngridk(1)/2+wf_kset%ngridk(1)-1
            ir = ir + 1
            rptl( :, ir) = (/ dble( ix), dble( iy), dble( iz)/)
          end do
        end do
      end do
      
      ! calculate Hamlitonian matrix elements in Wannier representation 
      allocate( ueu( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iy, auxmat)
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
    
      ! calculate phases for Fourier transform on Wigner-Seitz supercell
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iq, ir)
!$OMP DO
#endif
      do ir = 1, nrpt
        do iq = 1, wfint_kset%nkpt
          call ws_weight( rptl( :, ir), rptl( :, ir), wfint_kset%vkl( :, iq), wfint_pqr( iq, ir), kgrid=.true.)
        end do
        do ik = 1, wf_kset%nkpt
          call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), wfint_pkr( ik, ir), kgrid=.true.)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      allocate( auxmat( wfint_kset%nkpt, wf_kset%nkpt))
      call zgemm( 'N', 'C', wfint_kset%nkpt, wf_kset%nkpt, nrpt, zone, &
           wfint_pqr, wfint_kset%nkpt, &
           wfint_pkr, wf_kset%nkpt, zzero, &
           auxmat, wfint_kset%nkpt)
      wfint_phase = dble( transpose( auxmat))/wf_kset%nkpt
      deallocate( auxmat)
    
      ! calculate interpolated Hamiltonian
      allocate( hamilton( wf_fst:wf_lst, wf_fst:wf_lst, wfint_kset%nkpt))
#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iy)
!!$OMP DO
#endif
      do iy = wf_fst, wf_lst
        call zgemm( 'N', 'N', wf_nst, wfint_kset%nkpt, wf_kset%nkpt, zone, &
             ueu( iy, :, :), wf_nst, &
             cmplx( wfint_phase, 0, 8), wf_kset%nkpt, zzero, &
             hamilton( iy, :, :), wf_nst)
      end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
      deallocate( ueu)
    
      ! diagonalize interpolated Hamiltonian
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq)
!$OMP DO
#endif
      do iq = 1, wfint_kset%nkpt 
        call diaghermat( wf_nst, hamilton( :, :, iq), wfint_eval( :, iq), wfint_transform( :, :, iq))
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      deallocate( rptl, hamilton)
    
      return
    end subroutine wfint_interpolate_eigsys
    !EOC

!--------------------------------------------------------------------------------------
    
    !BOP
    ! !ROUTINE: wfint_interpolate_occupancy
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_occupancy
      ! !USES:
      use mod_eigenvalue_occupancy, only: occmax
      use mod_charge_and_moment, only: chgval
      ! !DESCRIPTION:
      !   Calclulates the interpolated occupation numbers for the wannierized bands and
      !   interpolated Fermi energy.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC

      integer, parameter :: maxit = 1000
      
      integer :: iq, ist, nvm, it
      real(8) :: e0, e1, fermidos, chg, x, t1
      
      real(8) :: sdelta, stheta
    
      allocate( wfint_occ( wf_fst:wf_lst, wfint_kset%nkpt))

      if( input%groundstate%stypenumber .ge. 0 ) then
        t1 = 1.d0/input%groundstate%swidth
        nvm = nint( chgval/occmax)
        if( (wf_fst .ne. 1) .or. (nvm .ge. wf_lst)) then
          write( *, '("Error (wfint_interpolate_occupancy): Please wannierize all occupied and at least one unoccupied band for the calculation of occupation numbers.")')
          call terminate
        end if
        ! check for insulator or semiconductor
        e0 = maxval( wfint_eval( nvm, :))
        e1 = minval( wfint_eval( nvm+1, :))
        wfint_efermi = 0.5*(e0 + e1)
    
        fermidos = 0.d0
        chg = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist, x) reduction(+: chg, fermidos)
!$OMP DO  
#endif
        do iq = 1, wfint_kset%nkpt
          do ist = wf_fst, wf_lst
            x = (wfint_eval( ist, iq) - wfint_efermi)*t1
            fermidos = fermidos + wfint_kset%wkpt( iq)*sdelta( input%groundstate%stypenumber, x)*t1
            wfint_occ( ist, iq) = occmax*stheta( input%groundstate%stypenumber, -x)
            chg = chg + wfint_kset%wkpt( iq)*wfint_occ( ist, iq)
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        fermidos = fermidos+occmax
        if( (e1 .ge. e0) .and. (abs( chg - chgval) .lt. input%groundstate%epsocc)) then
        !  write( *, '("Info (wannier_interpolate_occupancy): System has gap. Simplistic method used in determining efermi and occupation")')
        else
        ! metal found
          e0 = wfint_eval( 1, 1)
          e1 = e0
          do ist = wf_fst, wf_lst
            e0 = min( e0, minval( wfint_eval( ist, :)))
            e1 = max( e1, maxval( wfint_eval( ist, :)))
          end do
    
          it = 0
          do while( it .lt. maxit)
            wfint_efermi = 0.5*(e0 + e1)
            chg = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist, x) reduction(+: chg)
!$OMP DO  
#endif
            do iq = 1, wfint_kset%nkpt
              do ist = wf_fst, wf_lst
                x = (wfint_efermi - wfint_eval( ist, iq))*t1
                wfint_occ( ist, iq) = occmax*stheta( input%groundstate%stypenumber, x)
                chg = chg + wfint_kset%wkpt( iq)*wfint_occ( ist, iq)
              end do
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            if( chg .lt. chgval) then
              e0 = wfint_efermi
            else
              e1 = wfint_efermi
            end if
            if( (e1-e0) .lt. input%groundstate%epsocc) then
              it = maxit+1
            else
              it = it + 1
            end if
          end do
        end if
        if( it .eq. maxit) then
          write( *, '("Error (wfint_interpolate_occupancy): Fermi energy could not be found.")')
          call terminate
        end if
      else
        write( *, '("Error (wfint_interpolate_occupancy): Not implemented for this stype.")')
        call terminate
      end if
    
      return
    end subroutine wfint_interpolate_occupancy
    !EOC

!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfint_interpolate_bandgap
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_bandgap( factor_)
      ! !USES:
      ! !DESCRIPTION:
      !   Finds the position of VBM and CBM, the fundamental gap and the electron effective-mass tensor
      !
      ! !REVISION HISTORY:
      !   Created October 2017 (SeTi)
      !EOP
      !BOC

      integer, optional, intent( in) :: factor_

      integer, parameter :: maxit = 1000
      integer, parameter :: nline = 1000
      real(8), parameter :: infh = 1.d-3

      integer :: factor, ik, ist, i, j, n, imax, iv(3), ik0
      integer :: istvbm, istcbm
      integer :: lapack_info, lapack_ipiv(10)
      real(8) :: evbm, ecbm, kvbm(3), kcbm(3), dir(3,3), v(3)
      real(8) :: inv(10,10), coeff(10,1), energy(10,1), meff(3,3)
      real(8) :: a, amin, amax, eold, demax, d, length, l1, l2
      type( k_set) :: tmp_kset
      
      factor = 5
      if( present( factor_)) then
        if( factor_ .gt. 0) factor = factor_
      end if

      ! interpolate energies on fine grid
      call generate_k_vectors( tmp_kset, bvec, wf_kset%ngridk*factor, (/0.d0, 0.d0, 0.d0/), .true.)
      call wfint_init( tmp_kset)
      call wfint_interpolate_occupancy

      write(*,*) wfint_efermi
      ! find guess for VBM and CBM
      evbm = -1.d23
      ecbm = 1.d23
      do ik = 1, wfint_kset%nkpt
        do ist = wf_fst, wf_lst
          if( wfint_eval( ist, ik) .lt. wfint_efermi) then 
            if( wfint_eval( ist, ik) .gt. evbm) then
              evbm = wfint_eval( ist, ik)
              kvbm = wfint_kset%vkl( :, ik)
            end if
          else
            if( wfint_eval( ist, ik) .lt. ecbm) then
              ecbm = wfint_eval( ist, ik)
              kcbm = wfint_kset%vkl( :, ik)
            end if
          end if
        end do
      end do
      write(*,'(F23.16)') evbm
      write(*,'(3F23.16)') kvbm
      write(*,'(F23.16)') ecbm
      write(*,'(3F23.16)') kcbm
      write(*,*)

      call generate_k_vectors( tmp_kset, bvec, (/1, 1, nline+1/), (/0.d0, 0.d0, 0.d0/), .false.)

      ! find correct VBM
      n = 0
      d = 1.d0
      length = 1.0d0
      do while( d .gt. input%structure%epslat)
        n = n + 1
        if( mod( n, 3) .eq. 1) then
          dir = 0.d0
          dir(1,1) = 1.d0
          dir(2,2) = 1.d0
          dir(3,3) = 1.d0
        end if

        v = kvbm
        demax = 0.d0
        write(*,'(3F13.6)') dir(1,:)
        write(*,'(3F13.6)') dir(2,:)
        write(*,'(3F13.6)') dir(3,:)
        write(*,*)
        do i = 1, 3
          eold = evbm
          do ik = 0, nline
            tmp_kset%vkl( :, ik+1) = v + (-0.5d0*length + ik*length/nline)*dir(:,i)
          end do
          call wfint_init( tmp_kset)
          ik0 = -1
          do ik = 1, wfint_kset%nkpt
            do ist = wf_fst, wf_lst
              if( (wfint_eval( ist, ik) .lt. wfint_efermi) .and. (wfint_eval( ist, ik) .gt. evbm)) then
                evbm = wfint_eval( ist, ik)
                v = wfint_kset%vkl( :, ik)
                ik0 = ik-1
              end if
            end do
          end do
          if( evbm - eold .ge. demax) then
            demax = evbm - eold
            imax = i
          end if
          write(*,'(I,F23.16)') i, evbm - eold
        end do
        write(*,*) imax
        dir( :, imax) = v - kvbm
        dir( :, imax) = dir( :, imax)/norm2( dir( :, imax))

        write(*,'("-----------------")')
        do ik = 0, nline
          tmp_kset%vkl( :, ik+1) = kvbm + (-0.5d0*length + ik*length/nline)*dir(:,imax)
        end do
        call wfint_init( tmp_kset)
        ik0 = -1
        do ik = 1, wfint_kset%nkpt
          do ist = wf_fst, wf_lst
            if( (wfint_eval( ist, ik) .lt. wfint_efermi) .and. (wfint_eval( ist, ik) .gt. evbm)) then
              evbm = wfint_eval( ist, ik)
              v = wfint_kset%vkl( :, ik)
              istvbm = ist
              ik0 = ik-1
            end if
          end do
        end do
        !call r3frac( input%structure%epslat, v, iv)
        d = norm2( kvbm - v)
        kvbm = v
        write(*,'(F23.16)') evbm
        write(*,'(3F23.16)') kvbm
        write(*,'(F23.16)') d
        write(*,'("-----------------")')
        write(*,*)
        length = 0.5d0*length
      end do

      write(*,'("-----------------")')
      write(*,'("VBM at:   ",3F13.6)') kvbm
      write(*,'("E at VBM: ",F13.6," (band ",I3.3,")")') evbm, istvbm
      write(*,'("CBM at:   ",3F13.6)') kcbm
      write(*,'("E at CBM: ",F13.6," (band ",I3.3,")")') ecbm, istcbm
      write(*,'("fundamental gap: ",F13.6)') ecbm-evbm

      ! calculate electron effective mass tensor
      call generate_k_vectors( tmp_kset, bvec, (/1, 1, 10/), (/0.d0, 0.d0, 0.d0/), .false.)
      do i = 0, 20
        write(*,'(F23.16)') infh*0.5d0**i
        write(*,'(3F23.16)') kvbm
        write(*,*)
        l2 = 0.5d0**0*infh
        tmp_kset%vkl( :, 1) = kvbm + l2             *(/ 1,  0,  0/)
        tmp_kset%vkl( :, 2) = kvbm + l2             *(/ 0,  1,  0/)
        tmp_kset%vkl( :, 3) = kvbm + l2             *(/ 0,  0,  1/)
        tmp_kset%vkl( :, 4) = kvbm + l2/dsqrt( 2.d0)*(/ 1,  1,  0/)
        tmp_kset%vkl( :, 5) = kvbm + l2/dsqrt( 2.d0)*(/ 1,  0,  1/)
        tmp_kset%vkl( :, 6) = kvbm + l2/dsqrt( 2.d0)*(/ 0,  1,  1/)
        tmp_kset%vkl( :, 7) = kvbm + l2/dsqrt( 3.d0)*(/-1,  1,  1/)
        tmp_kset%vkl( :, 8) = kvbm + l2/dsqrt( 3.d0)*(/ 1, -1,  1/)
        tmp_kset%vkl( :, 9) = kvbm + l2/dsqrt( 3.d0)*(/ 1,  1, -1/)
        tmp_kset%vkl( :, 10) = kvbm
        call wfint_init( tmp_kset)
        evbm = wfint_eval( istvbm, 10)
        do ik = 1, 10
          energy( ik, 1) = wfint_eval( istvbm, ik)
          inv( ik, 1) = wfint_kset%vkl( 1, ik)**2
          inv( ik, 2) = wfint_kset%vkl( 2, ik)**2
          inv( ik, 3) = wfint_kset%vkl( 3, ik)**2
          inv( ik, 4) = 2.d0*wfint_kset%vkl( 1, ik)*wfint_kset%vkl( 2, ik)
          inv( ik, 5) = 2.d0*wfint_kset%vkl( 1, ik)*wfint_kset%vkl( 3, ik)
          inv( ik, 6) = 2.d0*wfint_kset%vkl( 2, ik)*wfint_kset%vkl( 3, ik)
          inv( ik, 7) = wfint_kset%vkl( 1, ik)
          inv( ik, 8) = wfint_kset%vkl( 2, ik)
          inv( ik, 9) = wfint_kset%vkl( 3, ik)
          inv( ik, 10) = 1.d0
          write(*,'(F23.16)') energy( ik, 1)
        end do
        call dgesv( 10, 1, inv, 10, lapack_ipiv, energy, 10, lapack_info)
        if( lapack_info .ne. 0) then
          write(*,'("Aborted! INFO: ",I)') lapack_info
          stop
        end if
        !energy = energy/infh/infh
        dir(1,1) = energy(1,1)
        dir(2,2) = energy(2,1)
        dir(3,3) = energy(3,1)
        dir(1,2) = energy(4,1)
        dir(2,1) = energy(4,1)
        dir(1,3) = energy(5,1)
        dir(3,1) = energy(5,1)
        dir(2,3) = energy(6,1)
        dir(3,2) = energy(6,1)
        v(1) = energy(7,1)
        v(2) = energy(8,1)
        v(3) = energy(9,1)
        call r3minv( 2.d0*dir, meff)
        call r3mv( meff, -v, dir(:,1))
        write(*,*)
        write(*,'(3F13.6)') meff(1,:)
        write(*,'(3F13.6)') meff(2,:)
        write(*,'(3F13.6)') meff(3,:)
        write(*,*)
        write(*,'(3F23.16)') dir(:,1)-kvbm
        write(*,*)
        write(*,'(F23.16)') energy(10,1)
        write(*,*)
        write(*,*)
        kvbm = dir(:,1)
      end do
      return
    end subroutine wfint_interpolate_bandgap
    !EOC

!--------------------------------------------------------------------------------------
    
    !BOP
    ! !ROUTINE: wfint_interpolate_occupancy
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_density
      ! !USES:
      ! !DESCRIPTION:
      !   Calclulates the interpolated occupation numbers for the wannierized bands and
      !   interpolated Fermi energy.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC

      integer :: ik, iq, igk, ir, irc, tp, ist, jst, is, ia, ias, ngknr, ig, ifg
      integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3, o1, o2, lmo1, lmo2, ilo1, ilo2
      integer :: lmaxapw, nlmomax
      real(8) :: gaunt, gnt, r1

      integer, allocatable :: idxgnt(:,:,:), nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:), lm2l(:)
      real(8), allocatable :: listgnt(:,:,:)

      complex(8), allocatable :: rhomt_int(:,:,:), rhoir_int(:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), match_combined(:,:,:)
      complex(8), allocatable :: dmat(:,:,:), dmatr(:,:,:,:), dmatq(:,:), wfrir(:,:,:), wfqir(:,:)

      complex(8), allocatable :: auxmat(:,:), auxmat4(:,:,:,:), zfft(:,:), auxmat3(:,:,:)
      
      external :: gaunt

      call timesec( r1)
      !write(*,*) r1
      lmaxapw = input%groundstate%lmaxapw
      ! count combined (l,m,o) indices and build index maps
      allocate( nlmo( nspecies))
      allocate( lmo2l( (lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2m( (lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2o( (lmaxapw + 1)**2*apwordmax, nspecies))
      nlmomax = 0
      do is = 1, nspecies
        nlmo( is) = 0
        do l1 = 0, lmaxapw
          do o1 = 1, apword( l1, is)
            do m1 = -l1, l1
              nlmo( is) = nlmo( is) + 1
              lmo2l( nlmo( is), is) = l1
              lmo2m( nlmo( is), is) = m1
              lmo2o( nlmo( is), is) = o1
            end do
          end do
        end do
        nlmomax = max( nlmomax, nlmo( is))
      end do
      if( nlmomax .gt. lmmaxapw*apwordmax) then
          write(*,*) "ERROR (wfint_interpolate_density): wrong nlmomax"
          write(*,*) nlmomax
          write(*,*) lmmaxapw*apwordmax
          write(*,*) lmaxapw, lmmaxapw, apwordmax
      end if

      allocate( lm2l( (lmaxapw + 1)**2))
      do l1 = 0, lmaxapw
        do m1 = -l1, l1
          lm1 = idxlm( l1, m1)
          lm2l( lm1) = l1
        end do
      end do

      ! build non-zero Gaunt list
      allocate( idxgnt( lmaxapw + 1, (lmaxapw + 1)**2, (lmaxapw + 1)**2), &
                listgnt( lmaxapw + 1, (lmaxapw + 1)**2, (lmaxapw + 1)**2))
      idxgnt(:,:,:) = 0
      listgnt(:,:,:) = 0.d0
      do l1 = 0, lmaxapw
        do m1 = -l1, l1
          lm1 = idxlm( l1, m1)

          do l2 = 0, lmaxapw
            do m2 = -l2, l2
              lm2 = idxlm( l2, m2)

              ik = 0
              do l3 = 0, lmaxapw
                do m3 = -l3, l3
                  lm3 = idxlm( l3, m3)
                  gnt = gaunt( l1, l2, l3, m1, m2, m3)
                  if( gnt .ne. 0.d0) then
                    ik = ik + 1
                    listgnt( ik, lm1, lm2) = gnt
                    idxgnt( ik, lm1, lm2) = lm3
                  end if
                end do
              end do

            end do
          end do
        
        end do
      end do

      call readstate
      call readfermi
      call linengy
      call genapwfr
      call genlofr
      call olprad

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      allocate( match_combined( nlmomax, ngkmax, natmtot))
      allocate( auxmat( nmatmax, wf_fst:wf_lst))
      allocate( auxmat4( wf_kset%nkpt, nlmomax+nlotot, wf_fst:wf_lst, natmtot))
      !allocate( auxmat3( wf_Gset%ngrtot, wf_kset%nkpt, wf_fst:wf_nst))
      allocate( zfft( wf_Gset%ngrtot, wf_fst:wf_lst))

      auxmat4 = zzero

      call timesec( r1)
      !write(*,*) r1

      ! build real space density matrix
      do ik = 1, wf_kset%nkpt
        ngknr = wf_Gkset%ngk( 1, ik)

        ! get matching coefficients
        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
          
        ! read eigenvector      
        if( input%properties%wannier%input .eq. "groundstate") then
          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
        else if( input%properties%wannier%input .eq. "hybrid") then
          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
        else if( input%properties%wannier%input .eq. "gw") then
          call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspnfv, evecfv)
        else
          call terminate
        end if

        match_combined = zzero
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            do lmo1 = 1, nlmo( is)
              l1 = lmo2l( lmo1, is)
              m1 = lmo2m( lmo1, is)
              o1 = lmo2o( lmo1, is)
              lm1 = idxlm( l1, m1)
              match_combined( lmo1, 1:ngknr, ias) = apwalm( 1:ngknr, o1, lm1, ias, 1)
            end do
          end do
        end do

        call zgemm( 'N', 'N', nmatmax, wf_nst, wf_nst, zone, &
             evecfv( :, wf_fst:wf_lst, 1), nmatmax, &
             wf_transform( :, :, ik), wf_nst, zzero, &
             auxmat, nmatmax)

        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            call zgemm( 'N', 'N', nlmo( is), wf_nst, ngknr, zone, &
                 match_combined( 1:nlmo( is), 1:ngknr, ias), nlmo( is), &
                 auxmat( 1:ngknr, :), ngknr, zzero, &
                 auxmat4( ik, 1:nlmo( is), :, ias), nlmo( is))
            do ilo1 = 1, nlorb( is)
              l1 = lorbl( ilo1, is)
              do m1 = -l1, l1
                lm1 = idxlm( l1, m1)
                auxmat4( ik, nlmo( is)+idxlo( lm1, ilo1, ias), :, ias) = auxmat( ngknr+idxlo( lm1, ilo1, ias), :)
              end do
            end do

          end do
        end do

!        zfft = zzero
!        do igk = 1, ngknr
!          ifg = igfft( wf_Gkset%igkig( igk, 1, ik))
!          zfft( ifg, :) = auxmat( igk, :)
!        end do
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist)
!!$OMP DO  
!#endif
!        do ist = wf_fst, wf_lst
!          call zfftifc( 3, ngrid, 1, zfft( :, ist))
!        end do
!#ifdef USEOMP                
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
!
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ig, r1, ifg)
!!$OMP DO  
!#endif
!        do ig = 1, wf_Gset%ngrtot
!          r1 = wf_kset%vkl( 1, ik)*wf_Gset%ivg( 1, ig)/ngrid(1) + &
!               wf_kset%vkl( 2, ik)*wf_Gset%ivg( 2, ig)/ngrid(2) + &
!               wf_kset%vkl( 3, ik)*wf_Gset%ivg( 3, ig)/ngrid(3)
!          ifg = igfft( ig)
!          auxmat3( ifg, ik, :) = zfft( ifg, :)*cmplx( cos( twopi*r1), sin( twopi*r1), 8)
!        end do      
!#ifdef USEOMP                
!!$OMP END DO
!!$OMP END PARALLEL
!#endif

      end do   

      call timesec( r1)
      !write(*,*) r1

      deallocate( evecfv, apwalm, auxmat, match_combined, zfft)
      allocate( dmatr( wf_kset%nkpt, nlmomax+nlotot, wf_fst:wf_lst, natmtot))

      dmatr = zzero

      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do ist = wf_fst, wf_lst
            call zgemm( 'C', 'N', wf_kset%nkpt, nlmo( is)+nlotot, wf_kset%nkpt, zone/wf_kset%nkpt, &
                 wfint_pkr, wf_kset%nkpt, &
                 auxmat4( :, 1:(nlmo( is)+nlotot), ist, ias), wf_kset%nkpt, zzero, &
                 dmatr( :, 1:(nlmo( is)+nlotot), ist, ias), wf_kset%nkpt)
          end do
        end do
      end do

      call timesec( r1)
      !write(*,*) r1

      deallocate( auxmat4)
      !allocate( wfrir( wf_Gset%ngrtot, wf_kset%nkpt, wf_fst:wf_nst))

!#ifdef USEOMP                
!!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist)
!!!$OMP DO  
!#endif
!      do ist = wf_fst, wf_lst
!        call zgemm( 'N', 'N', wf_Gset%ngrtot, wf_kset%nkpt, wf_kset%nkpt, zone/wf_kset%nkpt, &
!             auxmat3( :, :, ist), wf_Gset%ngrtot, &
!             conjg( wfint_pkr(:,:)), wf_kset%nkpt, zzero, &
!             wfrir( :, :, ist), wf_Gset%ngrtot)
!      end do
!#ifdef USEOMP                
!!!$OMP END DO
!!!$OMP END PARALLEL
!#endif

      call timesec( r1)
      !write(*,*) r1

      !deallocate( auxmat3)
      allocate( dmat( nlmomax+nlotot, nlmomax+nlotot, natmtot))
      allocate( rhoir_int( wf_Gset%ngrtot))

      dmat = zzero
      rhoir_int = zzero

      call timesec( r1)
      !write(*,*) r1

      call wfint_interpolate_occupancy

#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, is, ia, ias, ist, ir, auxmat, dmatq, zfft, wfqir), reduction(+:dmat, rhoir_int)
#endif
      allocate( auxmat( nlmomax+nlotot, wf_fst:wf_lst))
      allocate( dmatq( nlmomax+nlotot, wf_fst:wf_lst))
      allocate( zfft( wf_Gset%ngrtot, wf_fst:wf_lst))
      allocate( wfqir( wf_Gset%ngrtot, wf_fst:wf_lst))
#ifdef USEOMP                
!!$OMP DO  
#endif
      do iq = 1, wfint_kset%nkpt
        !write(*,*) iq
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            do ist = wf_fst, wf_lst
              call zgemv( 'T', wf_kset%nkpt, nlmomax+nlotot, zone, &
                   dmatr( :, :, ist, ias), wf_kset%nkpt, &
                   wfint_pqr( iq, :), 1, zzero, &
                   auxmat( :, ist), 1)
            end do
            call zgemm( 'N', 'N', nlmomax+nlotot, wf_nst, wf_nst, zone, &
                 auxmat, nlmomax+nlotot, &
                 wfint_transform( :, :, iq), wf_nst, zzero, &
                 dmatq, nlmomax+nlotot)
            do ist = wf_fst, wf_lst
              auxmat( :, ist) = dmatq( :, ist)*wfint_occ( ist, iq)
            end do
            call zgemm( 'N', 'C', nlmomax+nlotot, nlmomax+nlotot, wf_nst, cmplx( wfint_kset%wkpt( iq), 0, 8), &
                 dmatq, nlmomax+nlotot, &
                 auxmat, nlmomax+nlotot, zone, &
                 dmat( :, :, ias), nlmomax+nlotot)
          end do
        end do

!        do ist = wf_fst, wf_lst
!          call zgemv( 'N', wf_Gset%ngrtot, wf_kset%nkpt, zone, &
!               wfrir( :, :, ist), wf_Gset%ngrtot, &
!               wfint_pqr( iq, :), 1, zzero, &
!               zfft( :, ist), 1)
!        end do
!        call zgemm( 'N', 'N', wf_Gset%ngrtot, wf_nst, wf_nst, zone, &
!             zfft, wf_Gset%ngrtot, &
!             wfint_transform( :, :, iq), wf_nst, zzero, &
!             wfqir, wf_Gset%ngrtot)
!        do ist = wf_fst, wf_lst
!          do ir = 1, wf_Gset%ngrtot
!            rhoir_int( ir) = rhoir_int( ir) + wfint_kset%wkpt( iq)*wfint_occ( ist, iq)*(dble( wfqir( ir, ist))**2 + aimag( wfqir( ir, ist))**2)/omega
!          end do
!        end do
      end do
#ifdef USEOMP                
!!$OMP END DO
#endif
      deallocate( auxmat, dmatq, zfft, wfqir)
#ifdef USEOMP                
!!$OMP END PARALLEL
#endif

      deallocate( dmatr)

      call timesec( r1)
      !write(*,*) r1
      
      allocate( rhomt_int( lmmaxvr, nrmtmax, natmtot))

      rhomt_int = zzero
      
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do lmo1 = 1, nlmo( is)
            l1 = lmo2l( lmo1, is)
            m1 = lmo2m( lmo1, is)
            o1 = lmo2o( lmo1, is)
            lm1 = idxlm( l1, m1)
            ! APW-APW
            do lmo2 = 1, nlmo( is)
              l2 = lmo2l( lmo2, is)
              m2 = lmo2m( lmo2, is)
              o2 = lmo2o( lmo2, is)
              lm2 = idxlm( l2, m2)
              ik = 1
              do while( idxgnt( ik, lm1, lm2) .ne. 0)
                lm3 = idxgnt( ik, lm1, lm2)
                rhomt_int( lm3, :, ias) = rhomt_int( lm3, :, ias) + listgnt( ik, lm1, lm2)*dmat( lmo1, lmo2, ias)*&
                    apwfr( :, 1, o1, l1, ias)*apwfr( :, 1, o2, l2, ias)
                ik = ik + 1
              end do
            end do
            ! APW-LO
            do ilo2 = 1, nlorb( is)
              l2 = lorbl( ilo2, is)
              do m2 = -l2, l2
                lm2 = idxlm( l2, m2)
                ik = 1
                do while( idxgnt( ik, lm1, lm2) .ne. 0)
                  lm3 = idxgnt( ik, lm1, lm2)
                  rhomt_int( lm3, :, ias) = rhomt_int( lm3, :, ias) + listgnt( ik, lm1, lm2)*dmat( lmo1, nlmo( is)+idxlo( lm2, ilo2, ias), ias)*&
                      apwfr( :, 1, o1, l1, ias)*lofr( :, 1, ilo2, ias)
                  ik = ik + 1
                end do
              end do
            end do
          end do
          do ilo1 = 1, nlorb( is)
            l1 = lorbl( ilo1, is)
            do m1 = -l1, l1
              lm1 = idxlm( l1, m1)
              ! LO-LO
              do ilo2 = 1, nlorb( is)
                l2 = lorbl( ilo2, is)
                do m2 = -l2, l2
                  lm2 = idxlm( l2, m2)
                  ik = 1
                  do while( idxgnt( ik, lm1, lm2) .ne. 0)
                    lm3 = idxgnt( ik, lm1, lm2)
                    rhomt_int( lm3, :, ias) = rhomt_int( lm3, :, ias) + listgnt( ik, lm1, lm2)*dmat( nlmo( is)+idxlo( lm1, ilo1, ias), nlmo( is)+idxlo( lm2, ilo2, ias), ias)*&
                        lofr( :, 1, ilo1, ias)*lofr( :, 1, ilo2, ias)
                    ik = ik + 1
                  end do
                end do
              end do
              ! LO-APW
              do lmo2 = 1, nlmo( is)
                l2 = lmo2l( lmo2, is)
                m2 = lmo2m( lmo2, is)
                o2 = lmo2o( lmo2, is)
                lm2 = idxlm( l2, m2)
                ik = 1
                do while( idxgnt( ik, lm1, lm2) .ne. 0)
                  lm3 = idxgnt( ik, lm1, lm2)
                  rhomt_int( lm3, :, ias) = rhomt_int( lm3, :, ias) + listgnt( ik, lm1, lm2)*dmat( nlmo( is)+idxlo( lm1, ilo1, ias), lmo2, ias)*&
                      lofr( :, 1, ilo1, ias)*apwfr( :, 1, o2, l2, ias)
                  ik = ik + 1
                end do
              end do
            end do
          end do
        end do
      end do

      !call plotmat( rhomt_int(:,:,1))

      rhomt = dble( rhomt_int)
      rhoir = dble( rhoir_int)

      !call charge
      !write(*,'(100F13.6)') chgmt
      !write(*,'(F13.6)') chgir
      !write(*,*)

      call symrf( input%groundstate%lradstep, rhomt, rhoir)
      call rfmtctof( rhomt)
      !call gencore
      !call addrhocr

      call charge

      write(*,'(100F13.6)') chgmt
      !write(*,'(F13.6)') chgir
      !write(*,'(F13.6)') chgcalc
      !write(*,'(F13.6)') chgtot
      !write(*,*)

      !call rhonorm
      !call charge

      !write(*,'(100F13.6)') chgmt
      !write(*,'(F13.6)') chgir
      !write(*,'(F13.6)') chgcalc
      !write(*,'(F13.6)') chgtot
      !write(*,*)
      !write(*,*)
      
    end subroutine wfint_interpolate_density
    !EOC

!--------------------------------------------------------------------------------------
      
    subroutine wfint_interpolate_bandchar( lmax, bc)
      integer, intent( in) :: lmax
      real(8), intent( out) :: bc( 0:lmax, natmtot, wf_fst:wf_lst, wfint_kset%nkpt)

      integer :: ik, iq, ir, is, ia, ias, l, m, lm, o, ilo1, ilo2, lmmax, ngknr, maxdim
      integer :: lmcnt( 0:lmax, nspecies), o2idx( apwordmax, 0:lmax, nspecies), lo2idx( nlomax, 0:lmax, nspecies)
      real(8) :: f
      complex(8) :: z

      character(256) :: fname

      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), dmatk(:,:,:,:,:), dmatr(:,:,:,:,:), dmatq(:,:)
      real(8), allocatable :: rolpint(:,:,:,:), rholm(:,:,:,:)

      lmmax = (lmax + 1)**2

      fname = filext
      if( input%properties%wannier%input .eq. "hybrid") filext = '_PBE.OUT'
      call readstate
      filext = fname
      call readfermi
      call linengy
      call genapwfr
      call genlofr
      call olprad

      maxdim = 0
      do is = 1, nspecies
        o2idx( :, :, is) = 0
        lo2idx( :, :, is) = 0
        lmcnt( :, is) = 0
        do l = 0, lmax
          do o = 1, apword( l, is)
            lmcnt( l, is) = lmcnt( l, is) + 1
            o2idx( o, l, is) = lmcnt( l, is)
          end do
        end do
        do ilo1 = 1, nlorb( is)
          l = lorbl( ilo1, is)
          lmcnt( l, is) = lmcnt( l, is) + 1
          lo2idx( ilo1, l, is) = lmcnt( l, is)
        end do
        maxdim = max( maxdim, maxval( lmcnt( :, is)))
      end do

      ! radial overlap-intergrals
      allocate( rolpint( maxdim, maxdim, 0:lmax, natmtot))
      rolpint = 0.d0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          ! APW-APW integral
          do l = 0, lmax
            do o = 1, apword( l, is)
              rolpint( o2idx( o, l, is), o2idx( o, l, is), l, ias) = 1.d0
            end do
          end do
          do ilo1 = 1, nlorb( is)
            l = lorbl( ilo1, is)
            do o = 1, apword( l, is)
              if( (o2idx( o, l, is) .gt. 0) .and. (lo2idx( ilo1, l, is) .gt. 0)) then
                rolpint( o2idx( o, l, is), lo2idx( ilo1, l, is), l, ias) = oalo( o, ilo1, ias)
                rolpint( lo2idx( ilo1, l, is), o2idx( o, l, is), l, ias) = oalo( o, ilo1, ias)
              end if
            end do
            do ilo2 = 1, nlorb( is)
              if( lorbl( ilo2, is) .eq. l) then
                if( (lo2idx( ilo1, l, is) .gt. 0) .and. (lo2idx( ilo2, l, is) .gt. 0)) then
                  rolpint( lo2idx( ilo1, l, is), lo2idx( ilo2, l, is), l, ias) = ololo( ilo1, ilo2, ias)
                  rolpint( lo2idx( ilo2, l, is), lo2idx( ilo1, l, is), l, ias) = ololo( ilo2, ilo1, ias)
                end if
              end if
            end do
          end do
        end do
      end do
   
      allocate( dmatk( wf_fst:wf_lst, maxdim, lmmax, natmtot, wf_kset%nkpt))

      dmatk = zzero

      ! build k-point density coefficients
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, ngknr, apwalm, evecfv, is, ia, ias, l, m, lm, o, ilo1)
#endif
      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
#ifdef USEOMP                
!$OMP DO
#endif
      do ik = 1, wf_kset%nkpt
        ngknr = wf_Gkset%ngk( 1, ik)

        ! get matching coefficients
        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
          
#ifdef USEOMP                
!$OMP CRITICAL (readevec)
#endif
        ! read eigenvector      
        if( input%properties%wannier%input .eq. "groundstate") then
          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
        else if( input%properties%wannier%input .eq. "hybrid") then
          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
        else if( input%properties%wannier%input .eq. "gw") then
          call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspnfv, evecfv)
        else
          call terminate
        end if
#ifdef USEOMP                
!$OMP END CRITICAL (readevec)
#endif

        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            ! APW contribution
            do l = 0, lmax
              do m = -l, l
                lm = idxlm( l, m)
                do o = 1, apword( l, is)
                  call zgemv( 'T', ngknr, wf_nst, zone, &
                       evecfv( 1:ngknr, wf_fst:wf_lst, 1), ngknr, &
                       apwalm( 1:ngknr, o, lm, ias, 1), 1, zzero, &
                       dmatk( :, o2idx( o, l, is), lm, ias, ik), 1)
                end do
              end do
            end do
            ! LO contribution
            do ilo1 = 1, nlorb( is)
              l = lorbl( ilo1, is)
              do m = -l, l
                lm = idxlm( l, m)
                dmatk( :, lo2idx( ilo1, l, is), lm, ias, ik) = evecfv( ngknr+idxlo( lm, ilo1, ias), wf_fst:wf_lst, 1)
              end do
            end do
          end do
        end do
      end do
#ifdef USEOMP                
!$OMP END DO
#endif
      deallocate( evecfv, apwalm)
#ifdef USEOMP                
!$OMP END PARALLEL
#endif

      ! build R-point density coefficients
      allocate( dmatr( wf_fst:wf_lst, maxdim, lmmax, natmtot, wf_kset%nkpt))
      dmatr = zzero
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ir, ik, is, ia, ias, lm)
!$OMP DO
#endif
      do ir = 1, wf_kset%nkpt
        do ik = 1, wf_kset%nkpt
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              do lm = 1, lmmax
                call zgemm( 'T', 'N', wf_nst, maxdim, wf_nst, conjg( wfint_pkr( ik, ir))/wf_kset%nkpt, &
                     wf_transform( :, :, ik), wf_nst, &
                     dmatk( :, :, lm, ias, ik), wf_nst, zone, &
                     dmatr( :, :, lm, ias, ir), wf_nst)
              end do
            end do
          end do
        end do
      end do
#ifdef USEOMP                
!$OMP END DO
!$OMP END PARALLEL
#endif
      
      deallocate( dmatk)

      ! build q-point density coefficients
      allocate( dmatq( wf_fst:wf_lst, maxdim))
      bc = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ir, is, ia, ias, l, m, lm, ilo1, ilo2, dmatq) REDUCTION(+:bc)
!$OMP DO
#endif
      do iq = 1, wfint_kset%nkpt
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            do l = 0, lmax
              do m = -l, l
                lm = idxlm( l, m)
                ! calculate q-point density coefficient
                dmatq = zzero
                do ir = 1, wf_kset%nkpt
                  call zgemm( 'T', 'N', wf_nst, maxdim, wf_nst, wfint_pqr( iq, ir), &
                       wfint_transform( :, :, iq), wf_nst, &
                       dmatr( :, :, lm, ias, ir), wf_nst, zone, &
                       dmatq, wf_nst)
                end do
                ! add to bandcharacter
                do ilo1 = 1, lmcnt( l, is)
                  bc( l, ias, :, iq) = bc( l, ias, :, iq) + dble( conjg( dmatq( :, ilo1))*dmatq( :, ilo1))*rolpint( ilo1, ilo1, l, ias)
                  do ilo2 = ilo1+1, lmcnt( l, is)
                    bc( l, ias, :, iq) = bc( l, ias, :, iq) + 2.0d0*dble( conjg( dmatq( :, ilo1))*dmatq( :, ilo2))*rolpint( ilo1, ilo2, l, ias)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
#ifdef USEOMP                
!$OMP END DO
!$OMP END PARALLEL
#endif

      deallocate( dmatr, dmatq, rolpint)
      return
      
    end subroutine wfint_interpolate_bandchar

!--------------------------------------------------------------------------------------

    subroutine wfint_interpolate_me( mein, meout)
      complex(8), intent( in) :: mein( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt)
      complex(8), intent( out) :: meout( wf_fst:wf_lst, wf_fst:wf_lst, wfint_kset%nkpt)

      integer :: ir, ik, iq, ist
      complex(8), allocatable :: okwan(:,:,:), or(:,:,:), oqwan(:,:), ou(:,:)

      allocate( okwan( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))

#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, ir, ou)
#endif
      allocate( ou( wf_fst:wf_lst, wf_fst:wf_lst))
#ifdef USEOMP                
!$OMP DO  
#endif
      do ik = 1, wf_kset%nkpt
        call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
             mein( :, :, ik), wf_nst, &
             wf_transform( :, :, ik), wf_nst, zzero, &
             ou, wf_nst)
        call zgemm( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
             wf_transform( :, :, ik), wf_nst, &
             ou, wf_nst, zzero, &
             okwan( :, :, ik), wf_nst)
      end do
#ifdef USEOMP                
!$OMP END DO  
#endif
      deallocate( ou)
#ifdef USEOMP                
!$OMP END PARALLEL
#endif

      allocate( or( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))

      do ist = wf_fst, wf_lst
        call zgemm( 'N', 'N', wf_nst, wf_kset%nkpt, wf_kset%nkpt, zone/wf_kset%nkpt, &
             okwan( ist, :, :), wf_nst, &
             conjg( wfint_pkr), wf_kset%nkpt, zzero, &
             or( ist, :, :), wf_nst)
      end do

      deallocate( okwan)

#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ir, ou, oqwan)
#endif
      allocate( ou( wf_fst:wf_lst, wf_fst:wf_lst))
      allocate( oqwan( wf_fst:wf_lst, wf_fst:wf_lst))
#ifdef USEOMP                
!$OMP DO  
#endif
      do iq = 1, wfint_kset%nkpt
        oqwan = zzero 
        do ir = 1, wf_kset%nkpt
          oqwan(:,:) = oqwan(:,:) + or( :, :, ir)*wfint_pqr( iq, ir)
        end do
        call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
             oqwan, wf_nst, &
             wfint_transform( :, :, iq), wf_nst, zzero, &
             ou, wf_nst)
        call zgemm( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
             wfint_transform( :, :, iq), wf_nst, &
             ou, wf_nst, zzero, &
             meout( :, :, iq), wf_nst)
      end do
#ifdef USEOMP                
!$OMP END DO  
#endif
      deallocate( ou, oqwan)
#ifdef USEOMP                
!$OMP END PARALLEL
#endif

      deallocate( or)

    end subroutine wfint_interpolate_me

!--------------------------------------------------------------------------------------

    subroutine wfint_interpolate_eigvec
      integer :: ik, iq, is, ia, ias, l, m, lm, o, lmo, nlmomax, ilo1, ilo2
      integer :: i1, i2, i3, i4, iv(3)
      integer :: ngk, ngq, igk, igq
      integer :: ngkmax_old

      integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:)
      complex(8), allocatable :: evecfv(:,:,:), matchk(:,:,:), matchq(:,:,:)
      complex(8), allocatable :: cuv(:,:), uv(:,:)
      complex(8), allocatable :: olp_qk(:,:), olp_al_k(:,:), olp_al_q(:,:), z_q(:,:), z_q_H(:,:)
      complex(8), allocatable :: lsvec(:,:), rsvec(:,:)
      real(8), allocatable :: olp_ll(:,:), sval(:)

      complex(8), allocatable :: apwalm(:,:,:,:), sfacgp(:,:)
      real(8), allocatable :: gpc(:), tpgpc(:,:)
      
      character(256) :: fname

      ! build G+k-point set for interpolation
      call generate_Gk_vectors( wfint_Gkset, wfint_kset, wf_Gset, wf_Gkset%gkmax)
      if( allocated( wfint_evec)) deallocate( wfint_evec)
      ngkmax_old = ngkmax
      ngkmax = max( wf_Gkset%ngkmax, wfint_Gkset%ngkmax)

      ! count combined (l,m,o) indices and build index maps
      allocate( nlmo( nspecies))
      allocate( lmo2l( (input%groundstate%lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2m( (input%groundstate%lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2o( (input%groundstate%lmaxapw + 1)**2*apwordmax, nspecies))
      nlmomax = 0
      do is = 1, nspecies
        nlmo( is) = 0
        do l = 0, input%groundstate%lmaxapw
          do o = 1, apword( l, is)
            do m = -l, l
              nlmo( is) = nlmo( is) + 1
              lmo2l( nlmo( is), is) = l
              lmo2m( nlmo( is), is) = m
              lmo2o( nlmo( is), is) = o
            end do
          end do
        end do
        nlmomax = max( nlmomax, nlmo( is))
      end do

      ! generate radial overlaps
      ! WARNING: make sure, the correct state is read in advance according to wannier%input type
      !call genapwfr
      !call genlofr
      call olprad

      allocate( wfint_evec( wfint_Gkset%ngkmax+nlotot, wf_fst:wf_lst, wfint_kset%nkpt))
      wfint_evec = zzero

      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
      allocate( gpc (ngkmax))
      allocate( tpgpc (2, ngkmax))
      allocate( sfacgp (ngkmax, natmtot))

      allocate( evecfv( nmatmax, nstsv, nspinor))
      allocate( matchk( nlmomax, wf_Gkset%ngkmax, natmtot))
      allocate( matchq( nlmomax, wfint_Gkset%ngkmax, natmtot))
      allocate( olp_qk( wfint_Gkset%ngkmax+nlotot, wf_Gkset%ngkmax+nlotot))
      allocate( olp_al_k( nlotot, wf_Gkset%ngkmax))
      allocate( olp_al_q( wfint_Gkset%ngkmax, nlotot))
      allocate( olp_ll( nlotot, nlotot))
      allocate( uv( wf_fst:wf_lst, wf_fst:wf_lst))
      allocate( cuv( nmatmax, wf_fst:wf_lst))
      allocate( z_q( wfint_Gkset%ngkmax+nlotot, wf_fst:wf_lst))
      allocate( z_q_h( wf_fst:wf_lst, wfint_Gkset%ngkmax+nlotot))
      allocate( sval( wf_fst:wf_lst))
      allocate( lsvec( wf_fst:wf_lst, wf_fst:wf_lst))
      allocate( rsvec( 1:(wfint_Gkset%ngkmax+nlotot), 1:(wfint_Gkset%ngkmax+nlotot)))

      ! build k- and q-independent LO-LO overlap
      olp_ll = 0.d0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do ilo1 = 1, nlorb( is)
            l = lorbl( ilo1, is)
            do ilo2 = 1, nlorb( is)
              if( l .eq. lorbl( ilo2, is)) then
                i1 = idxlm( l, -l)
                i3 = idxlo( i1, ilo1, ias)
                i4 = idxlo( i1, ilo2, ias)
                do i2 = idxlm( l, -l), idxlm( l, l)
                  olp_ll( i3+i2-i1, i4+i2-i1) = ololo( ilo1, ilo2, ias)
                end do
                !do m = -l, l
                !  lm = idxlm( l, m)
                !  i1 = idxlo( lm, ilo1, ias)
                !  i2 = idxlo( lm, ilo2, ias)
                !  olp_ll( i1, i2) = ololo( ilo1, ilo2, ias)
                !end do
              end if
            end do
          end do
        end do
      end do

      do iq = 1, wfint_kset%nkpt
        z_q = zzero
        ngq = wfint_Gkset%ngk( 1, iq)
        ! find matching coefficients for G+q
        gpc( 1:ngq) = wfint_Gkset%gkc( 1:ngq, 1, iq)
        tpgpc( :, 1:ngq) = wfint_Gkset%tpgkc( :, 1:ngq, 1, iq)
        sfacgp( 1:ngq, :) = wfint_Gkset%sfacgk( 1:ngq, :, 1, iq)
        call match( ngq, gpc, tpgpc, sfacgp, apwalm)
        ! build q-dependent APW-LO overlap and combined matching-coefficients
        matchq = zzero
        olp_al_q = zzero
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            do ilo1 = 1, nlorb( is)
              l = lorbl( ilo1, is)
              i1 = idxlm( l, -l)
              i2 = idxlm( l, l)
              i3 = idxlo( i1, ilo1, ias)
              i4 = idxlo( i2, ilo1, ias)
              do o = 1, apword( l, is)
                olp_al_q( 1:ngq, i3:i4) = olp_al_q( 1:ngq, i3:i4) + conjg( apwalm( 1:ngq, o, i1:i2, ias))*oalo( o, ilo1, ias)
              end do
            end do
            do lmo = 1, nlmo( is)
              l = lmo2l( lmo, is)
              m = lmo2m( lmo, is)
              o = lmo2o( lmo, is)
              lm = idxlm( l, m)
              matchq( lmo, 1:ngq, ias) = apwalm( 1:ngq, o, lm, ias)
            end do
          end do
        end do
          
        do ik = 1, wf_kset%nkpt
          if( abs( wfint_phase( ik, iq)) .gt. input%structure%epslat) then
            write(*,*) ik
            olp_qk = zzero
            ngk = wf_Gkset%ngk( 1, ik)
            ! find matching coefficients for G+k
            gpc( 1:ngk) = wf_Gkset%gkc( 1:ngk, 1, ik)
            tpgpc( :, 1:ngk) = wf_Gkset%tpgkc( :, 1:ngk, 1, ik)
            sfacgp( 1:ngk, :) = wf_Gkset%sfacgk( 1:ngk, :, 1, ik)
            call match( ngk, gpc, tpgpc, sfacgp, apwalm)
            ! build k-dependent APW-LO overlap and combined matching-coefficients
            matchk = zzero
            olp_al_k = zzero
            do is = 1, nspecies
              do ia = 1, natoms( is)
                ias = idxas( ia, is)
                do ilo1 = 1, nlorb( is)
                  l = lorbl( ilo1, is)
                  do m = -l, l
                  lm = idxlm( l, m)
                  i1 = idxlo( lm, ilo1, ias)
                    do o = 1, apword( l, is)
                      olp_al_k( i1, 1:ngk) = olp_al_k( i1, 1:ngk) + apwalm( 1:ngk, o, lm, ias)*oalo( o, ilo1, ias)
                    end do
                  end do
                end do
                do lmo = 1, nlmo( is)
                  l = lmo2l( lmo, is)
                  m = lmo2m( lmo, is)
                  o = lmo2o( lmo, is)
                  lm = idxlm( l, m)
                  matchk( lmo, 1:ngk, ias) = apwalm( 1:ngk, o, lm, ias)
                end do
              end do
            end do
            ! build q-k-overlap
            do is = 1, nspecies
              do ia = 1, natoms( is)
                ias = idxas( ia, is)
                call zgemm( 'C', 'N', ngq, ngk, nlmo( is), zone, &
                     matchq( 1:nlmo( is), 1:ngq, ias), nlmo( is), &
                     matchk( 1:nlmo( is), 1:ngk, ias), nlmo( is), zone, &
                     olp_qk( 1:ngq, 1:ngk), ngq)
              end do
            end do
            if( norm2( wf_kset%vkl( :, ik) - wfint_kset%vkl( :, iq)) .lt. input%structure%epslat) then
              write(*,'(2I,3F13.6)') iq, ik, wf_kset%vkl( :, ik)
              do igk = 1, ngk
                do igq = 1, ngq
                  iv(:) = -wf_Gset%ivg( :, wf_Gkset%igkig( igk, 1, ik)) + wf_Gset%ivg( :, wfint_Gkset%igkig( igq, 1, iq))
                  i1 = ivgig( iv(1), iv(2), iv(3))
                  if( (i1 .gt. 0) .and. (i1 .le. ngvec)) then
                    olp_qk( igq, igk) = olp_qk( igq, igk) + cfunig( i1)
                  end if
                end do
              end do
            end if
            i1 = ngq+1
            i2 = ngq+nlotot
            i3 = ngk+1
            i4 = ngk+nlotot
            olp_qk( i1:i2, i3:i4) = zone*olp_ll(:,:)
            olp_qk( 1:ngq, i3:i4) = olp_al_q( 1:ngq, :)
            olp_qk( i1:i2, 1:ngk) = olp_al_k( :, 1:ngk)
            
            ! read eigenvector      
            if( input%properties%wannier%input .eq. "groundstate") then
              call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
            else if( input%properties%wannier%input .eq. "hybrid") then
              call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
            else if( input%properties%wannier%input .eq. "gw") then
              call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspinor, evecfv)
            else
              call terminate
            end if
            write( fname, '("evec/evec_cal_",I3.3)') ik
            call writematlab( evecfv( 1:i4, wf_fst:wf_lst, 1), fname)

            ! sum up Zq
            call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
                 wf_transform( :, :, ik), wf_nst, &
                 wfint_transform( :, :, iq), wf_nst, zzero, &
                 uv, wf_nst)
            call plotmat( uv)
            call zgemm( 'N', 'N', i4, wf_nst, wf_nst, zone, &
                 evecfv( 1:i4, wf_fst:wf_lst, 1), i4, &
                 uv, wf_nst, zzero, &
                 cuv( 1:i4, :), i4)
            call zgemm( 'N', 'N', i2, wf_nst, i4, cmplx( wfint_phase( ik, iq), 0, 8), &
                 olp_qk( 1:i2, 1:i4), i2, &
                 cuv( 1:i4, :), i4, zone, &
                 z_q( 1:i2, :), i2)
            if( norm2( wf_kset%vkl( :, ik) - wfint_kset%vkl( :, iq)) .lt. input%structure%epslat) then
              write( fname, '("olp/olp",3I3.3)') nint( wf_kset%vkl( :, ik)*1000)
              call writematlab( olp_qk( 1:i2, 1:i4), fname)
              call plotmat( uv, matlab=.true.)
              write(*,*) wfint_phase( ik, iq)
            end if
          end if
        end do !ik
        
        ! calculate interpolated eigenvectors
        i2 = ngq+nlotot
        z_q_h = conjg( transpose( z_q))
        call zgesdd_wrapper( z_q_h, wf_nst, i2, sval, lsvec, rsvec( 1:i2, 1:i2))
        lsvec = conjg( transpose( lsvec))
        rsvec = conjg( transpose( rsvec))
        do i1 = 1, wf_nst
          if( sval( i1) .lt. input%structure%epslat) then
            write(*,'(2I,F23.16)') iq, i1+wf_fst-1, sval( i1)
            rsvec( :, i1) = zzero
          else
            rsvec( :, i1) = rsvec( :, i1)/sval( i1)
          end if
        end do
        call zgemm( 'N', 'N', i2, wf_nst, wf_nst, zone, &
             rsvec( 1:i2, 1:wf_nst), i2, &
             lsvec, wf_nst, zzero, &
             wfint_evec( 1:i2, :, iq), i2)
        write( fname, '("evec/evec_int_",I3.3)') iq
        call writematlab( wfint_evec( 1:i2, :, iq), fname)
        !call zgemm( 'N', 'N', wf_nst, wf_nst, i2, zone, &
        !     z_q_h( :, 1:i2), wf_nst, &
        !     wfint_evec( 1:i2, :, iq), i2, zzero, &
        !     lsvec, wf_nst)
        !write(*,*) iq
        !call plotmat( lsvec)
      end do !iq

      deallocate( evecfv, apwalm, matchk, matchq, olp_qk, olp_al_k, olp_al_q, olp_ll, uv, cuv, z_q, z_q_h, sval, lsvec, rsvec, gpc, tpgpc, sfacgp) 
      ngkmax = ngkmax_old

    end subroutine wfint_interpolate_eigvec


!--------------------------------------------------------------------------------------
! helper functions

    subroutine wfint_putwfmt( ik, ias, wfmttp)
      use m_getunit
    
      integer, intent( in) :: ik, ias
      complex(8), intent( in) :: wfmttp( lmmaxvr, nrcmtmax, wf_fst:wf_lst)
    
      integer :: un, recl, offset
      character(256) :: filename
    
      inquire( iolength=recl) wf_kset%vkl( :, ik), ias, wf_fst, wf_lst, lmmaxvr, nrcmtmax, natmtot, wfmttp
    
      filename = 'WANNIER_WFMT.TMP'
      call getunit( un)
      open( un, file=filename, action='write', form='unformatted', access='direct', recl=recl)
      offset = (ik-1)*natmtot + ias
      write( un, rec=offset) wf_kset%vkl( :, ik), ias, wf_fst, wf_lst, lmmaxvr, nrcmtmax, natmtot, wfmttp
      close( un)
    
      return
    end subroutine wfint_putwfmt
    
    subroutine wfint_putwfir( ik, wfir)
      use m_getunit
    
      integer, intent( in) :: ik
      complex(8), intent( in) :: wfir( wf_Gset%ngrtot, wf_fst:wf_lst)
    
      integer :: un, recl
      character(256) :: filename
    
      inquire( iolength=recl) wf_kset%vkl( :, ik), wf_fst, wf_lst, wf_Gset%ngrtot, wfir
    
      filename = 'WANNIER_WFIR.TMP'
      call getunit( un)
      open( un, file=filename, action='write', form='unformatted', access='direct', recl=recl)
      write( un, rec=ik) wf_kset%vkl( :, ik), wf_fst, wf_lst, wf_Gset%ngrtot, wfir
      close( un)
    
      return
    end subroutine wfint_putwfir
    
    subroutine wfint_getwfmt( ik, ias, wfmttp)
      use m_getunit
    
      integer, intent( in) :: ik, ias
      complex(8), intent( out) :: wfmttp( lmmaxvr*nrcmtmax, wf_fst:wf_lst)
    
      integer :: i, un, recl, offset, fst, lst, ias_, lmmaxvr_, nrcmtmax_, natmtot_
      real(8) :: vl(3)
      character(256) :: filename
      logical :: exist
    
      inquire( iolength=recl) wf_kset%vkl( :, ik), ias, wf_fst, wf_lst, lmmaxvr, nrcmtmax, natmtot, wfmttp
    
      filename = 'WANNIER_WFMT.TMP'
      call getunit( un)
    
      do i = 1, 100
        inquire( file=filename, exist=exist)
        if( exist) then
          open( un, file=filename, action='read', form='unformatted', access='direct', recl=recl)
          exit
        else
          call system( 'sync')
          write(*,*) "Waiting for other process to write"
          call sleep( 1)
        end if
      end do
    
      offset = (ik-1)*natmtot + ias
      read( un, rec=offset) vl, ias_, fst, lst, lmmaxvr_, nrcmtmax_, natmtot_, wfmttp
      if( norm2( vl - wf_kset%vkl( :, ik)) .gt. input%structure%epslat) then
        write(*, '("Error (wannier_getwfir): differing vectors for k-point ",I8)') ik
        Write (*, '(" current	   : ", 3G18.10)') wf_kset%vkl (:, ik)
        Write (*, '(" WANNIER_WFMT.TMP : ", 3G18.10)') vl
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( (fst .ne. wf_fst) .or. (lst .ne. wf_lst)) then
        write(*, '("Error (wannier_getwfir): invalid band ranges")')
        Write (*, '(" current	   : ", 2I8)') wf_fst, wf_lst
        Write (*, '(" WANNIER_WFMT.TMP : ", 2I8)') fst, lst
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( ias_ .ne. ias) then
        write(*, '("Error (wannier_getwfir): differing atom-index")')
        Write (*, '(" current	   : ", I8)') ias
        Write (*, '(" WANNIER_WFMT.TMP : ", I8)') ias_
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( nrcmtmax_ .ne. nrcmtmax) then
        write(*, '("Error (wannier_getwfir): invalid number of radial points")')
        Write (*, '(" current	   : ", I8)') nrcmtmax
        Write (*, '(" WANNIER_WFMT.TMP : ", I8)') nrcmtmax_
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( lmmaxvr_ .ne. lmmaxvr) then
        write(*, '("Error (wannier_getwfir): invalid number of angular points")')
        Write (*, '(" current	   : ", I8)') lmmaxvr
        Write (*, '(" WANNIER_WFMT.TMP : ", I8)') lmmaxvr_
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( natmtot_ .ne. natmtot) then
        write(*, '("Error (wannier_getwfir): invalid number of atoms")')
        Write (*, '(" current	   : ", I8)') natmtot
        Write (*, '(" WANNIER_WFMT.TMP : ", I8)') natmtot_
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      close( un)
    
      return
    end subroutine wfint_getwfmt
    
    subroutine wfint_getwfir( ik, wfir)
      use m_getunit
    
      integer, intent( in) :: ik
      complex(8), intent( out) :: wfir( wf_Gset%ngrtot, wf_fst:wf_lst)
    
      integer :: i, un, recl, fst, lst, ng
      real(8) :: vl(3)
      character(256) :: filename
      logical :: exist
    
      inquire( iolength=recl) wf_kset%vkl( :, ik), wf_fst, wf_lst, wf_Gset%ngrtot, wfir
    
      filename = 'WANNIER_WFIR.TMP'
      call getunit( un)
    
      do i = 1, 100
        inquire( file=filename, exist=exist)
        if( exist) then
          open( un, file=filename, action='read', form='unformatted', access='direct', recl=recl)
          exit
        else
          call system( 'sync')
          write(*,*) "Waiting for other process to write"
          call sleep( 1)
        end if
      end do
    
      open( un, file=filename, action='read', form='unformatted', access='direct', recl=recl)
      read( un, rec=ik) vl, fst, lst, ng, wfir
      if( norm2( vl - wf_kset%vkl( :, ik)) .gt. input%structure%epslat) then
        write(*, '("Error (wannier_getwfir): differing vectors for k-point ",I8)') ik
        Write (*, '(" current	   : ", 3G18.10)') wf_kset%vkl (:, ik)
        Write (*, '(" WANNIER_WFIR.TMP : ", 3G18.10)') vl
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( (fst .ne. wf_fst) .or. (lst .ne. wf_lst)) then
        write(*, '("Error (wannier_getwfir): invalid band ranges")')
        Write (*, '(" current	   : ", 2I8)') wf_fst, wf_lst
        Write (*, '(" WANNIER_WFIR.TMP : ", 2I8)') fst, lst
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      if( ng .ne. wf_Gset%ngrtot) then
        write(*, '("Error (wannier_getwfir): invalid number of spatial points")')
        Write (*, '(" current	   : ", I8)') wf_Gset%ngrtot
        Write (*, '(" WANNIER_WFIR.TMP : ", I8)') ng
        Write (*, '(" file	  : ", a      )') filename
        stop
      end if
      close( un)
    
      return
    end subroutine wfint_getwfir
    
    subroutine wfint_destroywf
      use m_getunit
    
      integer :: un
      logical :: exist
    
      inquire( file='WANNIER_WFMT.TMP', exist=exist)
      if( exist) then
        call getunit( un)
        open( un, file='WANNIER_WFMT.TMP')
        close( un, status='delete')
      end if
    
      inquire( file='WANNIER_WFIR.TMP', exist=exist)
      if( exist) then
        call getunit( un)
        open( un, file='WANNIER_WFIR.TMP')
        close( un, status='delete')
      end if
    
      return
    end subroutine wfint_destroywf

end module mod_wfint
