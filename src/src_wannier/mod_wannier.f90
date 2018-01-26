module mod_wannier
  use mod_APW_LO
  use mod_atoms 
  use mod_kpoint
  use mod_constants
  use mod_muffin_tin
  use mod_Gvector
  use mod_Gkvector
  use mod_kpointset
  use mod_eigensystem, only: nmatmax, idxlo, oalo, ololo
  use mod_spin, only: nspinor, nspnfv
  use mod_lattice, only: bvec
  use mod_eigenvalue_occupancy, only: nstfv, nstsv, efermi
  use modgw, only: kset
  use m_getunit
  use m_ematqk
  use mod_pwmat
  use m_plotmat
  use modmpi
  use mod_misc
  use mod_symmetry, only: nsymcrys

  implicit none

! module variables
  integer :: wf_nprojtot                                ! total number of projection functions
  integer :: wf_nprojused                               ! number of local-orbitals used for projection
  integer :: wf_fst                                     ! lowest band used for generation
  integer :: wf_lst                                     ! highest band used for generation
  integer :: wf_nst                                     ! (maximum) number of bands used for generation
  integer :: wf_nwf                                     ! number of Wannier functions to generate
  integer :: wf_info                                    ! unit to write WANNIERINFO.OUT
  logical :: wf_initialized = .false.
  real(8) :: wf_t0
  type( k_set) :: wf_kset
  type( k_set) :: wf_kset_red
  type( G_set) :: wf_Gset
  type( Gk_set) :: wf_Gkset
  character(265) :: wf_filename
  
  integer, allocatable :: wf_projst(:,:)                ! l, m, species, atom of local-orbitals for projection
  integer, allocatable :: wf_projused(:)                ! flag wether a local-obital was used for projection
  complex(8), allocatable :: wf_transform(:,:,:)        ! unitary transformation matrices
  complex(8), allocatable :: wf_projection(:,:,:)       ! full overlap of wafefunctions and local-orbitals
  complex(8), allocatable :: wf_opf(:,:)                ! expansion coefficients for optimal projection functions
  complex(8), allocatable :: wf_pwmat(:,:,:,:)          ! original plane-wave matrix elements for neighboring k-points
  complex(8), allocatable :: wf_m0(:,:,:,:)             ! plane-wave matrix elements for neighboring k-points after subspace selection
  complex(8), allocatable :: wf_m(:,:,:,:)              ! transformed plane-wave matrix elements for neighboring k-points
  real(8), allocatable :: wf_centers(:,:)               ! centers of Wannier functions
  real(8), allocatable :: wf_omega(:)                   ! total localization functional for each Wannier function
  real(8), allocatable :: wf_omega_i(:)                 ! gauge independent localization functional for each Wannier function
  real(8), allocatable :: wf_omega_d(:)                 ! diagonal localization functional for each Wannier function
  real(8), allocatable :: wf_omega_od(:)                ! off-diagonal localization functional for each Wannier function
  real(8), allocatable :: wf_sheet(:,:,:)
  real(8), allocatable :: wf_rguide(:,:)

  ! disentanglement
  logical :: wf_disentangle = .false.
  real(8) :: wf_win_i(2), wf_win_o(2)                   ! boundaries of inner/outer window
  integer, allocatable :: wf_win_ii(:,:), wf_win_io(:,:)! bands inside inner/outer window per k-point
  integer, allocatable :: wf_win_ni(:), wf_win_no(:)    ! number of bands in inner/outer window per k-point
  complex(8), allocatable :: wf_subspace(:,:,:)         ! selected subspace for band-disentanglement
  complex(8), allocatable :: wf_sub_transform(:,:,:)
  real(8), allocatable :: wf_sub_eval(:,:)              ! Hamiltonian eigenvalues in subspace basis

  ! geometry
  integer :: wf_n_nshells                               ! number of shells used for gradient calculation
  integer :: wf_n_ntot                                  ! total numbers of neighbors in used shells
  real(8), allocatable :: wf_n_dist(:)                  ! distance of neighbors in given shell
  integer, allocatable :: wf_n_n(:)                     ! number of neighbors in given shell
  integer, allocatable :: wf_n_ik(:,:,:)                ! k-point index of given neighbor in given shell for given k-point
  integer, allocatable :: wf_n_ik2(:,:,:)               ! k-point index of given neighbor in given shell for given k-point
  integer, allocatable :: wf_n_usedshells(:)            ! index of used shells
  integer, allocatable :: wf_n_ns2n(:,:)                ! map from shell and neighbor in shell to total neighbor index
  integer, allocatable :: wf_n_n2neg(:,:)               ! map from shell and neighbor in shell to its negative
  real(8), allocatable :: wf_n_vl(:,:,:), wf_n_vc(:,:,:)! vector of given neighbor in given shell in lattice and cartesian coordinates
  real(8), allocatable :: wf_n_wgt(:)                   ! geometric weight of given shell

! methods
  contains
    !BOP
    ! !ROUTINE: wannier_init
    ! !INTERFACE:
    !
    subroutine wannier_init
      ! !USES:
      use mod_symmetry
      ! !DESCRIPTION:
      !   Reads local-orbitals from species files which are indicated to be used
      !   as projection functions for the generation of Wannier functions to the
      !   module {\tt mod\_wannier} and constructs necessary geometry (nearest 
      !   neighbors and geometric weights) as well as the overlap of the
      !   wavefunctions and the local-orbitals.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      ! local variables
      real(8) :: t0, t1

      character(256) :: fname
      
      write(*,*) "wannier_init"
      call init0
      call init1

      !********************************************************************
      ! open INFO file
      !********************************************************************
      write( wf_filename, '("WANNIER")')
      call getunit( wf_info)
      open( wf_info, file=trim( wf_filename)//"_INFO"//trim(filext), action='WRITE', form='FORMATTED')
      call timesec( wf_t0)

      !********************************************************************
      ! load used k-grid
      !********************************************************************
      call wannier_setkpts

      !********************************************************************
      ! check input parameters
      !********************************************************************
      call wannier_readinput

      !********************************************************************
      ! find geometry
      !********************************************************************
      call wannier_geometry
      
      !********************************************************************
      ! analyze energy windows for disentanglement
      !********************************************************************
      !if( maxval( abs( input%properties%wannier%outerwindow)) .gt. 0) then
      !  if( input%properties%wannier%input .eq. "gw") then
      !    if( (input%properties%wannier%nst .le. 0) .or. (input%properties%wannier%nst .gt. input%gw%nbgw - input%gw%ibgw + 1)) then
      !      write( *, '(" ERROR (wannier_init): number of Wannier functions (nst) out of range (",I4,":",I4,").")'), &
      !          1, input%gw%nbgw - input%gw%nbgw + 1
      !      call terminate
      !    end if
      !  else
      !    if( (input%properties%wannier%nst .le. 0) .or. (input%properties%wannier%nst .gt. nstfv)) then
      !      write( *, '(" ERROR (wannier_init): number of Wannier functions (nst) out of range (",I4,":",I4,").")'), &
      !          1, nstfv
      !      call terminate
      !    end if
      !  end if
      !  wf_nwf = input%properties%wannier%nst
      !  wf_fst = minval( input%properties%wannier%outerwindow)
      !  wf_lst = maxval( input%properties%wannier%outerwindow)
      !  wf_nst = wf_lst - wf_fst + 1
      !  wf_win_o(1) = wf_fst
      !  wf_win_o(2) = wf_lst
      !  wf_win_no = wf_win_o(2) - wf_win_o(1) + 1
      !  if( wf_win_o(2) .eq. wf_win_o(1)) then
      !    write( *, '(" ERROR (wannier_init): The given outer energy-window is too narrow.")')
      !    call terminate
      !  end if
      !  if( maxval( abs( input%properties%wannier%innerwindow)) .gt. 0) then
      !    wf_win_i(1) = minval( input%properties%wannier%innerwindow)
      !    wf_win_i(2) = maxval( input%properties%wannier%innerwindow)
      !    wf_win_ni = wf_win_i(2) - wf_win_i(1) + 1
      !    if( wf_win_i(2) .eq. wf_win_i(1)) then
      !      write( *, '(" ERROR (wannier_init): The given inner energy-window is too narrow.")')
      !      call terminate
      !    end if
      !    if( wf_win_ni .gt. wf_nwf) then
      !      write( *, '(" ERROR (wannier_init): The inner energy-window must not contain more than ,"I4", bands.")'), wf_nwf
      !      call terminate
      !    end if
      !    if( (wf_win_i(1) .lt. wf_win_o(1)) .or. (wf_win_i(2) .gt. wf_win_o(2))) then
      !      write( *, '(" ERROR (wannier_init): The inner energy-window is not fully contained inside the outer energy-window.")')
      !      call terminate
      !    end if
      !  else
      !    wf_win_i = 0
      !    wf_win_ni = 0
      !  end if
      !  if( wf_win_o(2) - wf_win_o(1) + 1 .lt. wf_nwf) then
      !    write( *, '(" ERROR (wannier_init): The outer energy-window must contain at least ",I4," bands.")') wf_nwf
      !    call terminate
      !  end if
      !end if
      !if( input%properties%wannier%method .eq. "disentangle") wf_disentangle = .true.
      
      !write( wf_info, '(" calculate projection overlap matrices...")')
      call timesec( t0)

      call wannier_genradfun
      !********************************************************************
      ! auto remove linear dependent projectors
      !********************************************************************
      if( (input%properties%wannier%method .ne. "fromfile") .and. (input%properties%wannier%method .ne. "maxfromfile")) then
        call wannier_removeldprojf( eps=1.d-4)
      end if
        
      call wannier_writeinfo_lo
      call wannier_writeinfo_geometry
      call wannier_writeinfo_task( input%properties%wannier%method)

      !********************************************************************
      ! build full projection overlap matrices
      !********************************************************************
      select case (input%properties%wannier%method)
        case( "pro")
          call wannier_projection
        case( "opf")
          call wannier_projection
        case( "opfmax")
          call wannier_projection
        case( "promax")
          call wannier_projection
        case( "disentangle")
          call wannier_projection
        case default
      end select

      !********************************************************************
      ! calculate plane-wave matrix elements (if necesseary)
      !********************************************************************
      select case (input%properties%wannier%method)
        case( "opf")
          call wannier_emat
        case( "opfmax")
          call wannier_emat
        case( "promax")
          call wannier_emat
        case( "maxfromfile")
          call wannier_emat
        case( "disentangle")
          call wannier_emat
        case default
      end select
      
      wf_initialized = .true.

      if( wf_disentangle) call wannier_subspace
      
      return
    end subroutine wannier_init
    !EOC

    !BOP
    ! !ROUTINE: wannier_emat
    ! !INTERFACE:
    !
    subroutine wannier_emat
      ! !USES:
      ! !DESCRIPTION:
      !   Generates plane-wave matrix elements for neighboring k-points
      !   needed in order to calculate well/maximally localized Wannier
      !   functions.
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      ! local variables
      integer :: i, is, n, isym, k1, k2, iknr, ngknr, maxn, idxn, cntk
      real(8) :: t0, t1
      logical :: success

      ! allocatable arrays
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)
      
      write(*,*) "wannier_emat..."
      write( wf_info, '(" calculate plane-wave matrix-elements...")')
      call timesec( t0)

      ! check for existing file
      success = .true.
      inquire( file=trim( wf_filename)//"_EMAT"//trim( filext), exist=success)
      if( success) then
        call wannier_reademat( success)
        if( success) then
          call timesec( t1)
          write( wf_info, '(" ...plane-wave matrix-elements read from file.")')
          write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
          write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
          write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
        end if
      end if

      if( .not. success) then          
        if( allocated( wf_pwmat)) deallocate( wf_pwmat)
        allocate( wf_pwmat( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt, wf_n_ntot))

        wf_m0(:,:,:,:) = zzero

        !write(*,'("shape of M0:     ",4I6)') shape( wf_pwmat)
        !write(*,'("shape of evec:   ",3I6)') nmatmax, nstfv, nspinor
        !write(*,'("size of evec:    ",F13.6," GB")') nmatmax*nstfv*nspinor*16*1.d-9
        !write(*,'("shape of apwalm: ",4I6)') ngkmax, apwordmax, lmmaxapw, natmtot
        !write(*,'("size of apwalm:  ",F13.6," GB")') ngkmax*apwordmax*lmmaxapw*natmtot*16*1.d-9
        
        call pwmat_init( input%groundstate%lmaxapw, 8)

        do i = 1, wf_n_nshells
          is = wf_n_usedshells( i) 
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            call pwmat_init_qg( wf_n_vl( :, n, is), (/0, 0, 0/))
            cntk = 0
            write(*,'(1a1,"  neighbor ",i2," of ",i2,": ",i3,"%",$)') char(13), idxn, wf_n_ntot, nint( 100.d0*cntk/wf_kset%nkpt)
#ifdef USEOMP
!!$omp parallel default( shared) private( iknr, evecfv1, evecfv2)
#endif
        allocate( evecfv1( nmatmax, nstfv, nspinor))
        allocate( evecfv2( nmatmax, nstfv, nspinor))
#ifdef USEOMP
!!$omp do
#endif
            do iknr = 1, wf_kset%nkpt   
              !write(*,'(I)') iknr
              !call myfindkpt( wf_kset%vkl( :, iknr), wf_kset_red, isym, k1)
              !call myfindkpt( wf_kset_red%vkl( :, k1), wf_kset, isym, k1)

              !call myfindkpt( wf_kset%vkl( :, wf_n_ik( n, is, iknr)), wf_kset_red, isym, k2)
              !call myfindkpt( wf_kset_red%vkl( :, k2), wf_kset, isym, k2)

              !call pwmat_init_qg( wf_kset%vkl( :, k2) - wf_kset%vkl( :, k1), (/0, 0, 0/))
              k1 = iknr
              k2 = wf_n_ik( n, is, k1)
#ifdef USEOMP
!!$omp critical( readevec)
#endif
              ! read eigenvectors
              call wannier_getevec( k1, evecfv1)
              call wannier_getevec( k2, evecfv2)
#ifdef USEOMP
!!$omp end critical( readevec)
#endif
              ! generate plane-wave matrix elements
              call pwmat_genpwmat( wf_kset%vkl( :, k1), wf_kset%vkc( :, k1), wf_fst, wf_lst, wf_fst, wf_lst, &
                   evecfv1( :, wf_fst:wf_lst, 1), &
                   evecfv2( :, wf_fst:wf_lst, 1), &
                   wf_pwmat( :, :, iknr, idxn))
#ifdef USEOMP
!!$omp atomic update
#endif
              cntk = cntk + 1
#ifdef USEOMP
!!$omp end atomic
#endif
              write(*,'(1a1,"o neighbor ",i2," of ",i2,": ",i3,"%",$)') char(13), idxn, wf_n_ntot, nint( 100.d0*cntk/wf_kset%nkpt)
            end do
#ifdef USEOMP
!!$omp end do
#endif
        deallocate( evecfv1, evecfv2)
#ifdef USEOMP
!!$omp end parallel
#endif
            write(*,*)
          end do
        end do

        !write(*,*) "destroy"
        call pwmat_destroy
        !write(*,*) "write"
        call wannier_writeemat
        call timesec( t1)
        write( wf_info, '(" ...plane-wave matrix-elements calculated.")')
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
        write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
      end if

      write(*,*) "wannier_emat", wf_fst, wf_lst, wf_nst, wf_nwf
      if( .not. wf_disentangle) then
        if( allocated( wf_m0)) deallocate( wf_m0)
        allocate( wf_m0( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))
        do i = 1, wf_nwf
          do n = 1, wf_nwf
            wf_m0( i, n, :, :) = wf_pwmat( wf_fst+i-1, wf_fst+n-1, :, :)
          end do
        end do
      end if
      !do i = 1, wf_n_nshells
      !  is = wf_n_usedshells( i) 
      !  do n = 1, wf_n_n( is)
      !    idxn = wf_n_ns2n( n, is)
      !    do iknr = 1, wf_kset%nkpt
      !      k1 = iknr
      !      !call myfindkpt( wf_kset%vkl( :, iknr), wf_kset_red, isym, k1)
      !      !call myfindkpt( wf_kset_red%vkl( :, k1), wf_kset, isym, k1)

      !      !call myfindkpt( wf_kset%vkl( :, wf_n_ik( n, is, iknr)), wf_kset_red, isym, k2)
      !      !call myfindkpt( wf_kset_red%vkl( :, k2), wf_kset, isym, k2)
      !      !k2 = wf_n_ik( n, is, k1)
      !      write(*,'(I)') idxn
      !      write(*,'(I,3F13.6)') k1, wf_kset%vkl( :, k1)
      !      !write(*,'(I,3F13.6)') k2, wf_kset%vkl( :, k2)
      !      call plotmat( wf_m0( :, :, iknr, idxn))
      !      write(*,*)
      !      !call plotmat( wf_m0( :, :, iknr, idxn), matlab=.true.)
      !      !write(*,*)
      !    end do
      !  end do
      !end do

      write( wf_info, *)
      call flushifc( wf_info)
      !stop
      return
      !EOC
    end subroutine wannier_emat
    !EOP

    !BOP
    ! !ROUTINE: wannier_gen_pro
    ! !INTERFACE:
    !
    subroutine wannier_gen_pro
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      !   fst       : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   lst       : N, number of bands used for generation (in,integer)
      !   nproj     : N, number of projection orbitals used for generation (in,integer)
      !   loproj    : indices of local-orbitals used as projection functions
      !               (in, integer(nproj))
      ! !DESCRIPTION:
      !   Generates unitary transformation matrices $U_{\vec{k}}$ via the
      !   projection method from which a set of smooth Bloch-like states 
      !   $$ \left|\Psi_{m,\vec{k}}^W\right\rangle = \sum\limits_{n=n_1}^{n_2}
      !   \left(U_{\vec{k}}\right)_{nm} \left|\Psi_{n,\vec{k}}^H\right\rangle $$
      !   out of the $N=n_2-n_1+1$ Hamiltonian eigenstates $\Psi_{n,\vec{k}}^H$
      !   in the band range $[n_1,n_2]$ can be constructed and used for the
      !   generation of localized Wannier functions
      !   $$ \left| \vec{R}n\right\rangle = \frac{1}{N_k} \sum\limits_{\vec{k}}
      !   e^{-i\vec{k}\cdot\vec{R}} \left| \Psi_{n,\vec{k}}^W\right\rangle $$
      !   The matrices are stored in the module variable {\tt wf\_transform
      !   (complex(nst,nproj,nkptnr))}.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      ! local variables
      integer :: iproj, iknr, i, j
      real(8) :: t0, t1
      real(8), allocatable :: sval(:), ravg(:,:), omega(:)
      complex(8), allocatable :: projm(:,:), lsvec(:,:), rsvec(:,:)

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(" perform simple projection step...")')
      call timesec( t0)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      !********************************************************************
      ! build transformation matrices
      !********************************************************************

      allocate( projm( wf_nwf, wf_nwf), sval( wf_nwf), lsvec( wf_nwf, wf_nwf), rsvec( wf_nwf, wf_nwf))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_nwf, wf_nwf, wf_kset%nkpt))
         
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, i, j, projm, sval, lsvec, rsvec)
!$OMP DO  
#endif
      do iknr = 1, wf_kset%nkpt
        j = 0
        do i = 1, wf_nprojtot
          if( wf_projused( i) .eq. 1) then
            j = j+1
            projm( :, j) = wf_projection( wf_fst:wf_lst, i, iknr)
          end if
        end do
        call zgesdd_wrapper( projm, wf_nwf, wf_nwf, sval, lsvec, rsvec)
        ! for numerical stability
        do i = 1, wf_nwf
          if( sval( i) .lt. 1.d-12) lsvec( :, i) = zzero
        end do
        call ZGEMM( 'N', 'N', wf_nwf, wf_nwf, wf_nwf, zone, &
             lsvec, wf_nwf, &
             rsvec, wf_nwf, zzero, &
             wf_transform( :, :, iknr), wf_nwf)
      end do 
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      call wannier_loc
      call wannier_writefile
      deallocate( projm, sval, lsvec, rsvec)
      call timesec( t1)
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(" ...simple projection step performed.")')
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega)
      write( wf_info, *)
      call flushifc( wf_info)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      return
    end subroutine wannier_gen_pro
    !EOC
    
    !BOP
    ! !ROUTINE: wannier_gen_opf
    ! !INTERFACE:
    !
    subroutine wannier_gen_opf( fst, lst, nopf, bands, maxit, nowrite)
      ! !USES:
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

      integer, optional, intent( in) :: fst, lst, nopf, bands(:), maxit
      logical, optional, intent( in) :: nowrite
      ! local variables
      integer :: opf_fst, opf_lst, opf_nopf
      integer :: iknr, is, n, nx, iproj, i, j, l, ix, minit, ntot, idxn, ist, jst, a
      real(8) :: lambda, sum_n_wgt, phi, phimean, uncertainty, theta, evalmin, omegamin, limit, x1
      real(8) :: opf_min_q(3,3), opf_min_p(3,1), opf_min_c, opf_min_sol(3), eval(3), revec(3,3), opf_min_big(6,6), opf_min_aux(3,3), r, tp(2)
      real(8) :: t0, t1, t2
      complex(8) :: opf_min_z(3,1), evec(3,3), evec_big(6,6), eval_big(6), r1, r2, tmp
      logical :: opf_nowrite

      ! allocatable arrays
      integer, allocatable :: opf_bands(:)
      complex(8), allocatable :: opf_transform(:,:,:), opf_projm(:,:,:), lsvec(:,:), rsvec(:,:), opf_mixing(:,:), opf_x(:,:,:), opf_x0(:,:,:), opf_m(:,:)
      complex(8), allocatable :: auxmat(:,:), projm(:,:), lsvec2(:,:), rsvec2(:,:)
      real(8), allocatable :: sval(:), opf_t(:), sval2(:), phihist(:)

      opf_nowrite = .false.
      if( present( nowrite)) opf_nowrite = nowrite
      opf_fst = wf_fst
      if( present( fst)) opf_fst = fst
      opf_lst = wf_lst
      if( present( lst)) opf_lst = lst
      opf_nopf = wf_nwf
      if( present( nopf)) opf_nopf = nopf
      allocate( opf_bands( opf_nopf))
      do i = 1, opf_nopf
        opf_bands( i) = opf_fst + i - 1
      end do
      if( present( bands)) opf_bands = bands
      if( wf_nwf .ne. opf_nopf) write(*,*) " Warning (wannier_gen_opf): differing number of Wannier functions"

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
      ! build rectangular transformation matrices
      !********************************************************************
      allocate( opf_transform( opf_nopf, wf_nprojused, wf_kset%nkpt))
      allocate( opf_projm( opf_nopf, wf_nprojused, wf_kset%nkpt))
      allocate( lsvec( opf_nopf, opf_nopf), rsvec( wf_nprojused, wf_nprojused), sval( opf_nopf))
      write(*,'("used orbitals: ",I)') wf_nprojused
      write(*,'("band range:    ",3I)') opf_fst, opf_lst, opf_nopf
      write(*,'("bands:         ",200I)') opf_bands
      if( wf_disentangle) then
        allocate( auxmat( wf_fst:wf_lst, wf_nprojused))
        do iknr = 1, wf_kset%nkpt
          j = 0
          auxmat = zzero
          do i = 1, wf_nprojtot
            if( wf_projused( i) .eq. 1) then
              j = j+1
              auxmat( :, j) = wf_projection( wf_fst:wf_lst, i, iknr)
            end if
          end do
          call zgemm( 'c', 'n', opf_nopf, wf_nprojused, wf_nst, zone, &
               wf_subspace( :, :, iknr), wf_nst, &
               auxmat, wf_nst, zzero, &
               opf_projm( :, :, iknr), opf_nopf)
        end do
        deallocate( auxmat)
      else
        do iknr = 1, wf_kset%nkpt
          j = 0
          do i = 1, wf_nprojtot
            if( wf_projused( i) .eq. 1) then
              j = j+1
              do ist = 1, opf_nopf
                opf_projm( ist, j, iknr) = wf_projection( opf_bands( ist), i, iknr)
              end do
              !if( iknr .eq. 1) write(*,*) j, i
            end if
          end do
        end do
      end if
      do iknr = 1, wf_kset%nkpt
        call zgesdd_wrapper( opf_projm( :, :, iknr), opf_nopf, wf_nprojused, sval, lsvec, rsvec)
        ! for numerical stability
        !do ist = 1, wf_nst
        !  if( sval( ist) .lt. 1.d-12) lsvec( :, ist) = zzero
        !end do
        
        call zgemm( 'n', 'n', opf_nopf, wf_nprojused, opf_nopf, zone, &
             lsvec, opf_nopf, &
             rsvec, wf_nprojused, zzero, &
             opf_transform( :, :, iknr), opf_nopf)
      end do

      !call writematlab( opf_projm(:,:,12), 'pmopf')
      !call writematlab( opf_transform(:,:,12), 'tmopf')

      sum_n_wgt = 0.d0
      do i = 1, wf_n_nshells
        sum_n_wgt = sum_n_wgt + 2.d0*wf_n_n( wf_n_usedshells( i))*wf_n_wgt( wf_n_usedshells( i))
      end do
      write(*,'("sum wgt = ",F13.6)') sum_n_wgt
      lambda = sum_n_wgt*input%properties%wannier%lambdaopf
      !minit = nint( opf_nopf*(wf_nprojused - 0.5d0*(opf_nopf + 1)))
      !minit = min( minit, 5)
      minit = 2
      !limit = max( 1.d-3, input%properties%wannier%uncertainty)
      limit = 5.d-4

      !********************************************************************
      ! build enlarged overlap and constraint matrices
      !********************************************************************
      allocate( opf_x( wf_nprojused, wf_nprojused, wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_x0( wf_nprojused, wf_nprojused, wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_t( wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_m( opf_nopf, opf_nopf))
      ! build enlarged overlap matrices
      allocate( auxmat( opf_nopf, wf_nprojused))
      i = 0
      do iknr = 1, wf_kset%nkpt
        do j = 1, wf_n_nshells
          is = wf_n_usedshells( j)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            i = i + 1
            !do ist = 1, opf_nst
            !  do jst = 1, opf_nst
            !    opf_m( ist, jst) = wf_m0( opf_bands( ist), opf_bands( jst), iknr, idxn)
            !  end do
            !end do
            call zgemm( 'n', 'n', opf_nopf, wf_nprojused, opf_nopf, zone, &
                 wf_m0( :, :, iknr, idxn), opf_nopf, &
                 opf_transform( :, :, wf_n_ik( n, is, iknr)), opf_nopf, zzero, &
                 auxmat, opf_nopf)
            call zgemm( 'c', 'n', wf_nprojused, wf_nprojused, opf_nopf, zone, &
                 opf_transform( :, :, iknr), opf_nopf, &
                 auxmat, opf_nopf, zzero, &
                 opf_x0( :, :, i), wf_nprojused)
            opf_t( i) = -2.d0*wf_n_wgt( is)/wf_kset%nkpt
          end do
        end do
      end do
      ! build constraint matrices
      do iknr = 1, wf_kset%nkpt
        i = i + 1
        call zgemm( 'c', 'n', wf_nprojused, wf_nprojused, opf_nopf, zone, &
             opf_projm( :, :, iknr), opf_nopf, &
             opf_projm( :, :, iknr), opf_nopf, zzero, &
             opf_x0( :, :, i), wf_nprojused)
        do is = 1, wf_nprojused
          opf_x0( is, is, i) = opf_x0( is, is, i) - zone
        end do
        opf_t( i) = lambda/wf_kset%nkpt
      end do
      nx = i

      !********************************************************************
      ! minimize Lagrangian
      !********************************************************************
      if( allocated( wf_opf)) deallocate( wf_opf)
      allocate( wf_opf( wf_nprojused, opf_nopf))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( opf_nopf, opf_nopf, wf_kset%nkpt))
      wf_transform = zzero
      do i = 1, opf_nopf
        wf_transform( i, i, :) = zone
      end do

      allocate( opf_mixing( wf_nprojused, wf_nprojused))
      allocate( projm( opf_nopf, opf_nopf))
      allocate( sval2( opf_nopf))
      allocate( lsvec2( opf_nopf, opf_nopf))
      allocate( rsvec2( opf_nopf, opf_nopf))
      allocate( phihist( minit))
      deallocate( lsvec, rsvec)
      allocate( lsvec(3,3), rsvec(3,3))
      omegamin = 1.d13
      opf_mixing = zzero
      opf_x = opf_x0
      do i = 1, wf_nprojused
        opf_mixing( i, i) = zone
      end do
      n = 0
      phihist = 0.d0
      uncertainty = 1.d0
      call wannier_loc( initial=.true.)
      ! start minimization
      call timesec( t2)
      write(*,*) "opf: start minimization"
      write(*,'("nx = ",I6)') nx
      write(*,'("omega = ",4F13.6)') sum( wf_omega), sum( wf_omega_i), sum( wf_omega_od), sum( wf_omega_d)
      do while( (n .lt. minit) .or. (uncertainty .gt. limit))
        n = n + 1
        if( present( maxit)) then
          if( n .gt. maxit) exit
        end if
        do i = 1, opf_nopf
          do j = i+1, wf_nprojused
            opf_min_q = 0.d0
            opf_min_p = 0.d0
            opf_min_c = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, opf_min_z) reduction(+: opf_min_q, opf_min_p)
!$OMP DO  
#endif
            do is = 1, nx
              opf_min_z(1,1) = 0.5d0*(opf_x( i, i, is) - opf_x( j, j, is))
              opf_min_z(2,1) = -0.5d0*(opf_x( i, j, is) + opf_x( j, i, is))
              opf_min_z(3,1) = 0.5d0*zi*(opf_x( i, j, is) - opf_x( j, i, is))
              opf_min_q = opf_min_q + opf_t( is)*dble( matmul( opf_min_z, transpose( conjg( opf_min_z))))
              opf_min_p = opf_min_p + opf_t( is)*dble( conjg( opf_x( i, i, is) + opf_x( j, j, is))*opf_min_z)
              !opf_min_c = 0.25d0*opf_t( is)*abs( opf_x( i, i, is) + opf_x( j, j, is))**2
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            call diagsymmat( 3, opf_min_q, eval, revec)
            if( j .le. opf_nopf) then
              opf_min_sol(:) = revec( :, 1)
              if( opf_min_sol(1) .lt. 0.d0) opf_min_sol = -opf_min_sol
            else
              opf_min_big = 0.d0
              opf_min_big( 1:3, 1:3) = 2*opf_min_q
              opf_min_big( 1:3, 4:6) = 0.25d0*matmul( opf_min_p, transpose( opf_min_p)) - matmul( opf_min_q, opf_min_q)
              opf_min_big( 4, 1) = 1.d0
              opf_min_big( 5, 2) = 1.d0
              opf_min_big( 6, 3) = 1.d0
              call diaggenmat( 6, cmplx( opf_min_big, 0, 8), eval_big, evec_big)
              evalmin = maxval( abs( eval_big))
              do is = 1, 6
                if( abs( aimag( eval_big( is))) .lt. input%structure%epslat) evalmin = min( evalmin, dble( eval_big( is)))
              end do
              do is = 1, 3
                opf_min_q( is, is) = opf_min_q( is, is) - evalmin
              end do
              if( minval( eval(:) - evalmin) .lt. 1.d-6) then
                call zgesdd_wrapper( cmplx( opf_min_q, 0.d0, 8), 3, 3, eval, lsvec, rsvec)
                lsvec = conjg( transpose( lsvec))
                rsvec = conjg( transpose( rsvec))
                do is = 1, 3
                  if( eval( is) .gt. 1.d-6) then
                    rsvec( :, is) = rsvec( :, is)/eval( is)
                  else
                    rsvec( :, is) = zzero
                  end if
                end do
                opf_min_aux = dble( matmul( rsvec, lsvec))
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
              do a = 1, wf_nprojused
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
                do a = 1, wf_nprojused
                  tmp = opf_x( a, i, is)
                  opf_x( a, i, is) = tmp*r1 - opf_x( a, j, is)*conjg( r2)
                  opf_x( a, j, is) = opf_x( a, j, is)*r1 + tmp*r2
                end do
                do a = 1, wf_nprojused
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
        wf_opf(:,:) = opf_mixing( :, 1:opf_nopf)
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, projm, sval2, lsvec2, rsvec2)
!$OMP DO  
#endif
        do iknr = 1, wf_kset%nkpt
          call zgemm( 'n', 'n', opf_nopf, opf_nopf, wf_nprojused, zone, &
                 opf_projm( :, :, iknr), opf_nopf, &
                 wf_opf, wf_nprojused, zzero, &
                 projm, opf_nopf)
          !write(*,*) iknr
          !call plotmat( projm)
          call zgesdd_wrapper( projm, opf_nopf, opf_nopf, sval2, lsvec2, rsvec2)
          ! for numerical stability
          !do ist = 1, wf_nst
          !  if( sval2( ist) .lt. 1.d-12) lsvec2( :, ist) = zzero
          !end do
          !write(*,*) iknr
          call zgemm( 'n', 'n', opf_nopf, opf_nopf, opf_nopf, zone, &
                 lsvec2, opf_nopf, &
                 rsvec2, opf_nopf, zzero, &
                 wf_transform( :, :, iknr), opf_nopf)
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        call wannier_loc( totonly=.false.)

        phi = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, iknr) reduction(+:phi)
!$OMP DO  
#endif
        do is = 1, nx
          do ist = 1, opf_nopf
            phi = phi + opf_t( is)*opf_x( ist, ist, is)*conjg( opf_x( ist, ist, is))
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        phi = phi + opf_nopf*sum_n_wgt

!        if( sum( wf_omega) .gt. phi) then
!          input%properties%wannier%lambdaopf = input%properties%wannier%lambdaopf/2
!          write(*,'("auto-adjust lambda to ",F13.6," and restart")') input%properties%wannier%lambdaopf
!          omegamin = 1.d13
!          opf_mixing = zzero
!          opf_x = opf_x0
!          do i = 1, wf_nprojused
!            opf_mixing( i, i) = zone
!          end do
!          n = 0
!          phihist = 0.d0
!          uncertainty = 1.d0
!          opf_x = opf_x0
!          is = wf_kset%nkpt*wf_n_ntot + 1
!          ist = wf_kset%nkpt*( wf_n_ntot + 1)
!          opf_t( is:ist) = opf_t( is:ist)/2
!        end if

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

        ! convergence analysis
!        if( n .eq. 1) phihist(:) = sum( wf_omega)
!        phihist = cshift( phihist, -1)
!        phihist(1) = sum( wf_omega)
!        phimean = sum( phihist(:))/minit
!        if( n .gt. 1) then
!          uncertainty = sqrt( sum( (phihist(:)-phimean)**2)/(min( n, minit)-1))/abs( sum( wf_omega))
!        else
!          uncertainty = 1.d0
!        end if

        call timesec( t1)
        write(*,'(I7,3x,40F23.16)') n, t1-t2, phi, sum( wf_omega), sum( wf_omega_i), sum( wf_omega_od), sum( wf_omega_d), uncertainty
        !write(*,'(I7,3x,40F23.16)') n, t1-t2, phi, sum( wf_omega), uncertainty

      end do
      call writematlab( wf_opf, 'opf')
      if( .not. opf_nowrite) call wannier_writefile
      call timesec( t1)
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(" ...optimized projection functions (OPF) calculated.")')
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"minimum iterations: ",T40,7x,I6)') minit
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') n
      write( wf_info, '(5x,"cutoff uncertainty: ",T40,E13.6)') limit
      write( wf_info, '(5x,"achieved uncertainty: ",T40,E13.6)') uncertainty
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega)
      write( wf_info, *)
      call flushifc( wf_info)
      call writematlab( wf_opf, "opf")

!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      deallocate( lsvec2, rsvec2, projm, opf_x, auxmat, phihist, opf_mixing, opf_x0, opf_t, opf_transform, lsvec, rsvec, sval, sval2, opf_bands)
      return
      !EOC
    end subroutine wannier_gen_opf
    !EOP

    !BOP
    ! !ROUTINE: wannier_gen_max
    ! !INTERFACE:
    !
    subroutine wannier_maxloc
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      !   loproj    : indices of local-orbitals used as projection functions
      !               (in, integer(nband))
      ! !DESCRIPTION:
      !   Does the same thing as {\tt genwf} does but the transformation
      !   matrices are used for the generation of maximally localized Wannier
      !   functions. The matrices are computed in a self consistent loop. 
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      ! parameters
      integer :: minit, maxit, maxcg, maxls
      real(8) :: mixing

      ! local variables
      integer :: ix, iy, iz, i, n, is, iknr, ngknr, k1, k2, z0, z1, idxn, ncg, newls
      real(8) :: step, alpha, grad, dg1, dg2, nwgt, last_step, noise_level, eps1, eps2
      real(8) :: omegastart, omega, omegamean, uncertainty, omega_new, omega_old, omega_pred, omega_new1, omega_new2
      real(8) :: t0, t1, t2, v1(3), v2(3), r1, r2, r3, a, b, c, d, grad_norm, lower, upper
      real(8) :: lspar(2,3), m3(3,3), m3i(3,3), v3(3)
      complex(8) :: c1, c2
      logical :: success

      ! external
      complex(8) :: zdotc

      ! allocatable arrays
      integer, allocatable :: nn(:), integerlist(:), counter(:)
      real(8), allocatable :: ndist(:), nvl(:,:,:), nvc(:,:,:), bwgt(:)
      real(8), allocatable :: omegahist(:), gradhist(:), eval(:,:), eval2(:)
      complex(8), allocatable :: auxmat(:,:), evec(:,:), mlwf_r(:,:), mlwf_t(:,:), mlwf_grad(:,:,:), mlwf_grad_last(:,:,:), mlwf_transform(:,:,:), mlwf_dir(:,:,:), mlwf_dir_last(:,:,:)
      complex(8), allocatable :: noise(:,:)

      minit = input%properties%wannier%minit
      if( minit .le. 0) minit = min( 100, wf_nwf*(wf_nwf+1)/2)
      maxit = input%properties%wannier%maxit
      maxcg = min( 100, wf_nwf*wf_nwf)
      maxcg = min( 20, wf_nwf*(wf_nwf+1)/2)
      maxls = 20
      newls = 1
      mixing = input%properties%wannier%sl
      noise_level = input%properties%wannier%noise

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      allocate( auxmat( wf_nwf, wf_nwf))
      allocate( noise( wf_nwf, wf_nwf))

      write( wf_info, '(" minimize localization functional Omega...")')
      call timesec( t0)

      !********************************************************************
      ! initialize M matrices and Wannier centers
      !********************************************************************
      call wannier_loc
      wf_rguide = wf_centers
      call wannier_phases( firstcall=.true.)
      !wf_sheet = 0.d0
        
      !********************************************************************
      ! minimize localization functional
      !********************************************************************
      ! mixing parameter for self consistent minimization
      step = mixing
      last_step = step
      eps1 = 0.2
      eps2 = 2.0

      ! start minimization loop
      allocate( mlwf_transform( wf_nwf, wf_nwf, wf_kset%nkpt))
      allocate( evec( wf_nwf, wf_nwf), eval( wf_nwf, wf_kset%nkpt), eval2( wf_nwf))
      allocate( mlwf_r( wf_nwf, wf_nwf), &
                mlwf_t( wf_nwf, wf_nwf), &
                mlwf_grad( wf_nwf, wf_nwf, wf_kset%nkpt), &
                mlwf_grad_last( wf_nwf, wf_nwf, wf_kset%nkpt), &
                mlwf_dir( wf_nwf, wf_nwf, wf_kset%nkpt), &
                mlwf_dir_last( wf_nwf, wf_nwf, wf_kset%nkpt))
      allocate( omegahist( minit), gradhist( minit))
      allocate( integerlist( minit))
      allocate( counter( wf_kset%nkpt))
      do iz = 1, minit
        integerlist( iz) = iz
      end do
      iz = 0
      ncg = 0
      omegastart = sum( wf_omega)
      write(*,*) omegastart
      omega = omegastart
      omega_old = omega
      omegahist = 0.d0
      gradhist = 1.d0
      grad = 1.d0
      success = .false.
      uncertainty = 1.d0

      ! combined neighbor weights
      nwgt = 0.d0
      do i = 1, wf_n_nshells
        is = wf_n_usedshells( i)
        nwgt = nwgt + wf_n_wgt( is)*wf_n_n( is)
      end do
      nwgt = nwgt*4.0d0

      !************************************************************
      !* start of minimization
      !************************************************************
      call timesec( t2)
      do while( .not. success)
        iz = iz + 1
    
        call wannier_phases
        !do i = 1, wf_nwf
        !  write(*,'(I,5x,3F13.6,5x,3F13.6)') i, wf_centers( :, i), wf_rguide( :, i)
        !end do
        mlwf_grad = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, i, is, n, idxn, ix, a, mlwf_r, mlwf_t, k2)
!$OMP DO
#endif
        do iknr = 1, wf_kset%nkpt
          !-------------------------------------------------------
          !- calculate gradient
          !-------------------------------------------------------
          do i = 1, wf_n_nshells
            is = wf_n_usedshells( i)
            do n = 1, wf_n_n( is)
              idxn = wf_n_ns2n( n, is)
              k2 = wf_n_ik2( n, is, iknr)
              ! calculating R and T
              mlwf_r = wf_m( :, :, iknr, idxn)
              mlwf_t = wf_m( :, :, iknr, idxn)
              do ix = 1, wf_nwf
                a = 0.d0
                if( abs( wf_m( ix, ix, iknr, idxn)) .gt. 1.d-10) then
                  !a = atan2( dble( aimag( wf_m( ix, ix, iknr, idxn))), dble( wf_m( ix, ix, iknr, idxn)))
                  !a = aimag( log( wf_m( ix, ix, iknr, idxn)))
                  a = aimag( log( exp( zi*wf_sheet( ix, iknr, idxn))*wf_m( ix, ix, iknr, idxn))) - wf_sheet( ix, iknr, idxn)
                  !if( (a + dot_product( wf_n_vc( :, n, is), wf_centers( :, ix))) .gt. pi) a = a - twopi
                  !if( (a + dot_product( wf_n_vc( :, n, is), wf_centers( :, ix))) .lt. -pi) a = a + twopi
                  !a = aimag( log( exp( zi*wf_sheet( ix, iknr, idxn))*wf_m( ix, ix, iknr, idxn))) - wf_sheet( ix, iknr, idxn)
                  mlwf_r( :, ix) = mlwf_r( :, ix)*conjg( wf_m( ix, ix, iknr, idxn))
                  mlwf_t( :, ix) = mlwf_t( :, ix)/wf_m( ix, ix, iknr, idxn)*(a + dot_product( wf_n_vc( :, n, is), wf_centers( :, ix)))
                else
                  mlwf_r( :, ix) = zzero
                  mlwf_t( :, ix) = zzero
                end if
              end do
              mlwf_r = mlwf_r - conjg( transpose( mlwf_r))
              mlwf_t = mlwf_t + conjg( transpose( mlwf_t))
              ! calculating gradient
              mlwf_grad( :, :, iknr) = mlwf_grad( :, :, iknr) - wf_n_wgt( is)*( 0.5*mlwf_r(:,:) + 0.5*zi*mlwf_t(:,:))

              ! calculating R and T
              mlwf_r = conjg( transpose( wf_m( :, :, k2, idxn)))
              mlwf_t = conjg( transpose( wf_m( :, :, k2, idxn)))
              do ix = 1, wf_nwf
                a = 0.d0
                if( abs( wf_m( ix, ix, k2, idxn)) .gt. 1.d-10) then
                  !a = atan2( dble( aimag( wf_m( ix, ix, k2, idxn))), dble( wf_m( ix, ix, k2, idxn)))
                  !a = aimag( log( wf_m( ix, ix, k2, idxn)))
                  a = aimag( log( exp( zi*wf_sheet( ix, k2, idxn))*wf_m( ix, ix, k2, idxn))) - wf_sheet( ix, k2, idxn)
                  !if( (a + dot_product( wf_n_vc( :, n, is), wf_centers( :, ix))) .gt. pi) a = a - twopi
                  !if( (a + dot_product( wf_n_vc( :, n, is), wf_centers( :, ix))) .lt. -pi) a = a + twopi
                  !a = aimag( log( exp( zi*wf_sheet( ix, iknr, idxn))*wf_m( ix, ix, iknr, idxn))) - wf_sheet( ix, iknr, idxn)
                  mlwf_r( :, ix) = mlwf_r( :, ix)*wf_m( ix, ix, k2, idxn)
                  mlwf_t( :, ix) = -mlwf_t( :, ix)/conjg( wf_m( ix, ix, k2, idxn))*(a + dot_product( wf_n_vc( :, n, is), wf_centers( :, ix)))
                else
                  mlwf_r( :, ix) = zzero
                  mlwf_t( :, ix) = zzero
                end if
              end do
              mlwf_r = mlwf_r - conjg( transpose( mlwf_r))
              mlwf_t = mlwf_t + conjg( transpose( mlwf_t))
              ! calculating gradient
              mlwf_grad( :, :, iknr) = mlwf_grad( :, :, iknr) - wf_n_wgt( is)*( 0.5*mlwf_r(:,:) + 0.5*zi*mlwf_t(:,:))
            end do
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        !stop
        mlwf_grad = 4.0d0*mlwf_grad/wf_kset%nkpt
        grad_norm = sqrt( dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_grad, 1)))

        !-------------------------------------------------------
        !- calculate CG parameter
        !-------------------------------------------------------
        !if( (iz .gt. 1) .and. ( abs( omega_old - omega) .gt. input%properties%wannier%uncertainty)) then              
        if( (iz .gt. 1)) then              
          ! Hestenes-Stiefel
          if( input%properties%wannier%solver .eq. "hs") then
            r1 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_grad - mlwf_grad_last, 1))
            r2 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_dir_last, 1, mlwf_grad - mlwf_grad_last, 1))

          ! Fletcher-Reeves
          else if( input%properties%wannier%solver .eq. "fr") then
            r1 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_grad, 1))
            r2 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad_last, 1, mlwf_grad_last, 1))

          ! Polak-Ribier
          else if( input%properties%wannier%solver .eq. "pr") then
            r1 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_grad - mlwf_grad_last, 1))
            r2 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad_last, 1, mlwf_grad_last, 1))
          
          ! steepest descent
          else
            r1 = 0.d0
            r2 = 1.d0
          end if

          alpha = max( 0.d0, r1/r2)
          ncg = ncg + 1
        else
          alpha = 0.d0
          ncg = 0
        end if

        !-------------------------------------------------------
        !- add direction
        !-------------------------------------------------------
        mlwf_dir = -mlwf_grad + alpha*mlwf_dir_last

        r1 = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt

        ! reset if necessary
        if((mod( ncg, maxcg) .eq. 0) .or. (r1 .ge. 0.d0)) then
          mlwf_dir = -mlwf_grad
          r1 = -grad_norm*grad_norm/nwgt
          ncg = 0
        end if
        
        mlwf_grad_last = mlwf_grad
        mlwf_dir_last = mlwf_dir

#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr)
!!$OMP DO
#endif
        ! diagonalize update directions
        do iknr = 1, wf_kset%nkpt
          call diaghermat( wf_nwf, zi*mlwf_dir( :, :, iknr)/nwgt, eval( :, iknr), mlwf_grad( :, :, iknr))
        end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
        omega_old = sum( wf_omega)
        step = max( mixing, last_step)
        !step = mixing
        iy = 0
        if( input%properties%wannier%ls .and. (mod( iz-1, newls) .eq. 0)) then
          !-------------------------------------------------------
          !- perform line-search
          !-------------------------------------------------------

          lspar( 1, 1) = 0.d0
          lspar( 2, 1) = omega_old

          call wannier_update( step*eps2, eval, mlwf_grad, new=.true.)
          call wannier_loc( totonly=.true.)
          lspar( 1, 3) = step*eps2
          lspar( 2, 3) = sum( wf_omega)
          iy = iy + 1

          call wannier_update( step, eval, mlwf_grad)
          call wannier_loc( totonly=.true.)
          lspar( 1, 2) = step
          lspar( 2, 2) = sum( wf_omega)
          omega = lspar( 2, 2)
          iy = iy + 1

          ! minimum of parabolic fit
          m3( 1, :) = (/lspar( 1, 1)**2, lspar( 1, 1), 1.d0/)
          m3( 2, :) = (/lspar( 1, 2)**2, lspar( 1, 2), 1.d0/)
          m3( 3, :) = (/lspar( 1, 3)**2, lspar( 1, 3), 1.d0/)
          call r3minv( m3, m3i)
          call r3mv( m3i, lspar( 2, :)-minval( lspar( 2, :)), v3)
          a = -0.5d0*v3(2)/v3(1)
          !write(*,'(I,2F13.6,4(2F13.6,"  --  "))') iy, step, lspar( 2, 2), lspar( :, 1), lspar( :, 2), lspar( :, 3), a, -0.25d0*v3(2)**2/v3(1)+v3(3)

          ! condition 1 and 2 fulfilled
          if( (lspar( 2, 2) .le. omega_old + lspar( 1, 2)*eps1*r1) .and. (lspar( 2, 3) .gt. omega_old + lspar( 1, 3)*eps1*r1)) then
            ! do nothing

          ! condition 1 violated
          ! step too large
          else if( lspar( 2, 2) .gt. omega_old + lspar( 1, 2)*eps1*r1) then
            b = maxval( lspar( 1, :))     !maximum step length
            omega = lspar( 2, 2)
            step = lspar( 1, 2)

            do while( (iy .lt. 20) .and. (omega .gt. omega_old + step*eps1*r1))
              if( (a .le. 0.d0) .or. (a .gt. b) .or. (v3(1) .lt. 0.d0)) then
                step = step/eps2
              else
                step = a
              end if
              call wannier_update( step, eval, mlwf_grad)
              call wannier_loc( totonly=.true.)
              omega = sum( wf_omega)
              iy = iy + 1
              n = maxloc( lspar( 2, :), 1)  !position of Omega_max
              if( omega .le. lspar( 2, n)) then
                lspar( 1, n) = step
                lspar( 2, n) = omega
              end if
              m3( 1, :) = (/lspar( 1, 1)**2, lspar( 1, 1), 1.d0/)
              m3( 2, :) = (/lspar( 1, 2)**2, lspar( 1, 2), 1.d0/)
              m3( 3, :) = (/lspar( 1, 3)**2, lspar( 1, 3), 1.d0/)
              call r3minv( m3, m3i)
              call r3mv( m3i, lspar( 2, :)-minval( lspar( 2, :)), v3)
              a = -0.5d0*v3(2)/v3(1)
              b = maxval( lspar( 1, :))     !maximum step length
              !write(*,'(I,2F13.6,4(2F13.6,"  --  "))') iy, step, omega, lspar( :, 1), lspar( :, 2), lspar( :, 3), a, -0.25d0*v3(2)**2/v3(1)+v3(3)
              if( abs( a/step-1.d0) .lt. 1.d-3) exit
            end do
            last_step = step

          ! condition 2 violated
          ! step too short
          else if( lspar( 2, 3) .le. omega_old + lspar( 1, 3)*eps1*r1) then
            b = minval( lspar( 1, :))     !minimum step length
            omega = lspar( 2, 3)
            step = lspar( 1, 3)

            do while( (iy .lt. 20) .and. (omega .le. omega_old + step*eps1*r1))
              if( (a .le. b) .or. (v3(1) .lt. 0.d0)) then
                step = step*eps2
              else
                step = a
              end if
              call wannier_update( step, eval, mlwf_grad)
              call wannier_loc( totonly=.true.)
              omega = sum( wf_omega)
              iy = iy + 1
              n = maxloc( lspar( 2, :), 1)  !position of Omega_max
              if( omega .le. lspar( 2, n)) then
                lspar( 1, n) = step
                lspar( 2, n) = omega
              end if
              m3( 1, :) = (/lspar( 1, 1)**2, lspar( 1, 1), 1.d0/)
              m3( 2, :) = (/lspar( 1, 2)**2, lspar( 1, 2), 1.d0/)
              m3( 3, :) = (/lspar( 1, 3)**2, lspar( 1, 3), 1.d0/)
              call r3minv( m3, m3i)
              call r3mv( m3i, lspar( 2, :)-minval( lspar( 2, :)), v3)
              a = -0.5d0*v3(2)/v3(1)
              b = minval( lspar( 1, :))     !minimum step length
              !write(*,'(I,2F13.6,4(2F13.6,"  --  "))') iy, step, omega, lspar( :, 1), lspar( :, 2), lspar( :, 3), a, -0.25d0*v3(2)**2/v3(1)+v3(3)
              if( abs( a/step-1.d0) .lt. 1.d-3) exit
            end do
            last_step = step
          else
            write(*,*) "WTF"
            write(*,*) omega_old, omega_new1, omega_new2
            call terminate
          end if
          ! eventually check whether we have realy decreased omega. 
          ! Otherwhise go back to the smalest value tested. But change something.
          if( omega .ge. omega_old) then
            n = minloc( lspar( 2, :), 1)
            if( lspar( 1, n) .eq. 0.d0) then
              lspar( 2, n) = 1.d100
              n = minloc( lspar( 2, :), 1)
            end if
            call wannier_update( lspar( 1, n), eval, mlwf_grad)
            call wannier_loc( totonly=.true.)
            omega = sum( wf_omega)
            iy = iy + 1
            !write(*,'("CAUTION: ",3F13.6)') lspar( 1, n), lspar( 2, n), omega
          end if
        else
          call wannier_update( step, eval, mlwf_grad, new=.true.)
          iy = iy + 1
          call wannier_loc( totonly=.true.)
          omega = sum( wf_omega)
        end if

        if( noise_level .gt. input%properties%wannier%uncertainty) then
          noise = zone
          call random_seed()
          do ix = 2, wf_nwf
            do iy = 1, ix-1
              call random_number( r1)
              c1 = cmplx( cos( twopi*r1), sin( twopi*r1), 8)
              noise( ix, iy) = c1
              noise( iy, ix) = -conjg( c1)
            end do
          end do
          call diaghermat( wf_nwf, zi*noise, eval2, auxmat)
          do ix = 1, wf_nwf
            mlwf_dir( :, ix, 1) = exp( -zi*noise_level*eval2( ix))*auxmat( :, ix)
          end do
          call zgemm( 'N', 'C', wf_nwf, wf_nwf, wf_nwf, zone, &
               mlwf_dir( :, :, 1), wf_nwf, &
               auxmat, wf_nwf, zzero, &
               noise, wf_nwf)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, auxmat)
!$OMP DO  
#endif
          do iknr = 1, wf_kset%nkpt
            call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
                   wf_transform( :, :, iknr), wf_nwf, &
                   noise, wf_nwf, zzero, &
                   auxmat, wf_nwf)
            wf_transform( :, :, iknr) = auxmat
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        end if
        
        !write(*,'(3F13.6)') wf_centers( :, 12)
        
        !-------------------------------------------------------
        !- convergence analysis
        !-------------------------------------------------------
        if( iz .eq. 1) omegahist(:) = 0.d0
        omegahist = cshift( omegahist, -1)
        gradhist = cshift( gradhist, -1)
        omegahist(1) = omega
        omegamean = sum( omegahist(:))/min( minit, iz)
        if( iz .eq. 1) then
          uncertainty = 1.d0
          grad = 1.d0
        else
          uncertainty = sqrt( sum( (omegahist(:)-omegamean)**2)/(min( minit, iz)-1))/abs( omega)
          grad = dot_product( dble( integerlist( 1:min( minit, iz))), omegahist( 1:min( minit, iz))-omegamean) - (min( minit, iz)+1)*0.5d0*sum( omegahist( 1:min( minit, iz))-omegamean)
          grad = grad/sum( (dble( integerlist( 1:min( minit, iz)))-(min( minit, iz)+1)*0.5d0)**2)/abs( omega)
          gradhist(1) = grad
        end if
        v1 = (/gradhist(1), gradhist(3), gradhist(5)/)
        v2 = (/gradhist(2), gradhist(4), gradhist(6)/)

        !if( (minval( v1)*maxval( v1) .ge. 0.d0) .and. &
        !    (minval( v2)*maxval( v2) .ge. 0.d0) .and. &
        !    (minval( v1)+maxval( v1) .ge. 0.d0) .and. &
        !    (minval( v2)+maxval( v2) .le. 0.d0)) success = .true.
        if( omega .gt. omegastart) success = .true.
        if( uncertainty .le. input%properties%wannier%uncertainty) success = .true.
        if( maxval( gradhist) .le. input%properties%wannier%uncertainty) success = .true.
        if( iz .lt. minit) success = .false.
        if( grad_norm .lt. input%properties%wannier%uncertainty) success = .true.
        if( iz .ge. maxit) success = .true.
        call timesec( t1)
        write(*,'(I4,7F23.16,2I4)') iz, t1-t2, omega, omega_old-omega, grad_norm, step, r1, uncertainty, ncg, iy!, maxval( gradhist)!, omega2 !omegai+omegad+omegaod
      end do
      !************************************************************
      !* end of minimization
      !************************************************************

      write(*,*)
      if( omega .gt. omegastart) then
        write(*, '("ERROR (genmlwf): Localization functional diverged. Procedure aborted after ",I4," loops.")') iz
      else if( iz .ge. maxit) then
        write(*, '("ERROR (genmlwf): Not converged after ",I6," cycles.")') maxit
      else
        write(*,'(" success: Convergence reached after ",I4," cycles.")') iz 
        write(*,'(" Localization gain: ",I3,"%")') nint( 100d0*(omegastart-omega)/omega)
      end if
        
      call wannier_loc
      call wannier_writefile
  
      deallocate( auxmat, eval, evec, mlwf_r, mlwf_t, mlwf_grad, mlwf_grad_last, mlwf_dir, mlwf_dir_last, mlwf_transform, omegahist, gradhist)
      call timesec( t1)
      write( wf_info, '(" ...localization functional Omega minimized.")')
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"minimum/maximum iterations: ",T40,I6,"/",I6)') minit, maxit
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') iz
      write( wf_info, '(5x,"aimed uncertainty: ",T40,E13.6)') input%properties%wannier%uncertainty
      write( wf_info, '(5x,"achieved uncertainty: ",T40,E13.6)') uncertainty
      write( wf_info, '(5x,"norm of gradient: ",T40,E13.6)') grad_norm
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega)
      write( wf_info, '(5x,"localization gain: ",T40,7x, I5,"%")') nint( 100d0*(omegastart-sum( wf_omega))/sum( wf_omega))
      write( wf_info, *)
      call flushifc( wf_info)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      return
      !EOC
    end subroutine wannier_maxloc
    !EOP

    subroutine wannier_update( step, direval, direvec, new)
      real(8), intent( in) :: step, direval( wf_nwf, wf_kset%nkpt)
      complex(8), intent( in) :: direvec( wf_nwf, wf_nwf, wf_kset%nkpt)
      logical, optional, intent( in) :: new

      integer :: ik, ist
      real(8) :: used_step
      real(8), save :: last_step = 0.d0
      complex(8) :: auxmat1( wf_nwf, wf_nwf), auxmat2( wf_nwf, wf_nwf)
      logical :: update_step

      update_step = .true.
      if( present( new)) update_step = .not. new
      if( update_step) then
        used_step = step - last_step
      else
        used_step = step
      end if
      !write(*,'("update: ",3F13.6)') last_step, step, used_step 
      last_step = step

#ifdef USEOMP
!$omp parallel default( shared) private( ik, ist, auxmat1, auxmat2)
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        do ist = 1, wf_nwf
          auxmat1( :, ist) = exp( -zi*used_step*direval( ist, ik))*direvec( :, ist, ik) 
        end do
        call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
               wf_transform( :, :, ik), wf_nwf, &
               auxmat1, wf_nwf, zzero, &
               auxmat2, wf_nwf)
        call zgemm( 'n', 'c', wf_nwf, wf_nwf, wf_nwf, zone, &
               auxmat2, wf_nwf, &
               direvec( :, :, ik), wf_nwf, zzero, &
               wf_transform( :, :, ik), wf_nwf)
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      return
    end subroutine wannier_update

    subroutine wannier_projonwan
      !!!!!
      ! correct wf_nst and wf_nwf here
      !!!!!
      use m_wsweight
      integer :: ik, is, ia, ias, lmaxapw, l, l_, m, m_, lm, lm_, o, o_, lmo, lmo_, ir, ilo, ilo_, ngknr, ist, jst, ig, igk, ifg, nlmomax, nshell, nrpt
      real(8) :: x, fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)
      
      integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:)
      real(8), allocatable :: rptl(:,:)
      complex(8), allocatable :: wanfmt(:,:,:,:), wanfir(:,:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:), match_combined(:,:), wfmt(:,:), wfir(:,:), cfun(:), zfft(:), wswgt(:)
      complex(8), allocatable :: addmt(:,:), addir(:,:)

      lmaxapw = input%groundstate%lmaxapw
      nlmomax = (lmaxapw + 1)**2*apwordmax
      nshell = 1
      nrpt = (1 + 2*nshell)**3

      allocate( rptl( 3, nrpt))
      allocate( wswgt( nrpt))
      ir = 0
      do l = -nshell, nshell
        do m = -nshell, nshell
          do o = -nshell, nshell
            ir = ir + 1
            rptl( :, ir) = dble( (/o, m, l/))
          end do
        end do
      end do

      allocate( nlmo( nspecies))
      allocate( lmo2l( nlmomax, nspecies), &
                lmo2m( nlmomax, nspecies), &
                lmo2o( nlmomax, nspecies))
      allocate( wanfmt( nlmomax+nlotot, wf_fst:wf_lst, natmtot, nrpt))
      allocate( wanfir( wf_Gset%ngrtot, wf_fst:wf_lst, nrpt))
      
      call readstate
      call readfermi
      call linengy
      call genapwfr
      call genlofr
      call olprad
      call genidxlo

      call wannier_readfun( nshell, wanfmt, wanfir, nlmo, lmo2l, lmo2m, lmo2o)

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
      allocate( match_combined( nlmomax, ngkmax))
      allocate( wfmt( nlmomax+nlotot, wf_fst:wf_lst))
      allocate( wfir( wf_Gset%ngrtot, wf_fst:wf_lst))
      allocate( cfun( wf_Gset%ngrtot))
      allocate( zfft( wf_Gset%ngrtot))
      allocate( addmt( wf_fst:wf_lst, wf_nst))
      allocate( addir( wf_fst:wf_lst, wf_nst))

      write(*,*) nlmomax, nlotot
      write(*,*) shape( wanfmt)
      write(*,*) shape( wfmt)

      wf_projection = zzero
      wf_projused = 0
      wf_projused( 1:wf_nst) = 1

      cfun = zzero
      do ig = 1, wf_Gset%ngrtot
        cfun( igfft( ig)) = cfunig( ig)
      end do
      call zfftifc( 3, ngrid, 1, cfun)

      do ik = 1, wf_kset%nkpt
        ngknr = wf_Gkset%ngk( 1, ik)
        write(*,*) ik

        do ir = 1, nrpt
          call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), wswgt( ir), kgrid=.true.)
        end do

        ! get matching coefficients
        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)
          
        ! read eigenvector      
        if( input%properties%wannier%input .eq. "gs") then
          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
        else if( input%properties%wannier%input .eq. "hybrid") then
          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
        else if( input%properties%wannier%input .eq. "gw") then
          call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspnfv, evecfv)
        else
          call terminate
        end if

        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)

            match_combined = zzero
            wfmt = zzero

            do lmo = 1, nlmo( is)
              l = lmo2l( lmo, is)
              m = lmo2m( lmo, is)
              o = lmo2o( lmo, is)
              lm = idxlm( l, m)
              match_combined( lmo, 1:ngknr) = apwalm( 1:ngknr, o, lm, ias)
            end do

            call zgemm( 'N', 'N', nlmo( is), wf_nst, ngknr, zone, &
                 match_combined( 1:nlmo( is), 1:ngknr), nlmo( is), &
                 evecfv( 1:ngknr, wf_fst:wf_lst, 1), ngknr, zzero, &
                 wfmt( 1:nlmo( is), :), nlmo( is))

            do ilo = 1, nlorb( is)
              l = lorbl( ilo, is)
              do m = -l, l
                lm = idxlm( l, m)
                wfmt( nlmo( is)+idxlo( lm, ilo, ias), :) = evecfv( ngknr+idxlo( lm, ilo, ias), wf_fst:wf_lst, 1)
              end do
            end do

            do lmo = 1, nlmo( is)
              l = lmo2l( lmo, is)
              m = lmo2m( lmo, is)
              o = lmo2o( lmo, is)
              lm = idxlm( l, m)
              do ir = 1, nrpt
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!$OMP DO
#endif    
                do ist = wf_fst, wf_lst
                  do jst = wf_fst, wf_lst
                    wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
                        conjg( wfmt( lmo, ist))*wanfmt( lmo, jst, ias, ir)*wswgt( ir)
                  end do
                end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
              end do
              do ilo_ = 1, nlorb( is)
                l_ = lorbl( ilo_, is)
                do m_ = -l_, l_
                  lm_ = idxlm( l_, m_)
                  if( lm .eq. lm_) then
                    do ir = 1, nrpt
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!$OMP DO
#endif    
                      do ist = wf_fst, wf_lst
                        do jst = wf_fst, wf_lst
                          wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
                              conjg( wfmt( lmo, ist))*wanfmt( nlmo( is)+idxlo( lm_, ilo_, ias), jst, ias, ir)*cmplx( oalo( o, ilo_, ias), 0, 8)*wswgt( ir)
                        end do
                      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
                    end do
                  end if
                end do
              end do
            end do
            do ilo = 1, nlorb( is)
              l = lorbl( ilo, is)
              do m = -l, l
                lm = idxlm( l, m)
                do lmo_ = 1, nlmo( is)
                  l_ = lmo2l( lmo_, is)
                  m_ = lmo2m( lmo_, is)
                  o_ = lmo2o( lmo_, is)
                  lm_ = idxlm( l_, m_)
                  if( lm .eq. lm_) then
                    do ir = 1, nrpt
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!$OMP DO
#endif    
                      do ist = wf_fst, wf_lst
                        do jst = wf_fst, wf_lst
                          wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
                              conjg( wfmt( nlmo( is)+idxlo( lm, ilo, ias), ist))*wanfmt( lmo_, jst, ias, ir)*cmplx( oalo( o_, ilo, ias), 0, 8)*wswgt( ir)
                        end do
                      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
                    end do
                  end if
                end do
                do ilo_ = 1, nlorb( is)
                  l_ = lorbl( ilo_, is)
                  do m_ = -l_, l_
                    lm_ = idxlm( l_, m_)
                    if( lm .eq. lm_) then
                      do ir = 1, nrpt
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!$OMP DO
#endif    
                        do ist = wf_fst, wf_lst
                          do jst = wf_fst, wf_lst
                            wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
                                conjg( wfmt( nlmo( is)+idxlo( lm, ilo, ias), ist))*wanfmt( nlmo( is)+idxlo( lm_, ilo_, ias), jst, ias, ir)*cmplx( ololo( ilo, ilo_, ias), 0, 8)*wswgt( ir)
                          end do
                        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
                      end do
                    end if
                  end do
                end do
              end do
            end do                      
              
          end do
        end do

        wfir = zzero
        do igk = 1, ngknr
          ifg = igfft( wf_Gkset%igkig( igk, 1, ik))
          wfir( ifg, :) = evecfv( igk, wf_fst:wf_lst, 1)
        end do
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist)
!$OMP DO  
#endif
        do ist = wf_fst, wf_lst
          call zfftifc( 3, ngrid, 1, wfir( :, ist))
        end do
#ifdef USEOMP                
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ig, x, ifg)
!$OMP DO  
#endif
        do ig = 1, wf_Gset%ngrtot
          x = wf_kset%vkl( 1, ik)*wf_Gset%ivg( 1, ig)/ngrid(1) + &
              wf_kset%vkl( 2, ik)*wf_Gset%ivg( 2, ig)/ngrid(2) + &
              wf_kset%vkl( 3, ik)*wf_Gset%ivg( 3, ig)/ngrid(3)
          ifg = igfft( ig)
          wfir( ifg, :) = wfir( ifg, :)*cmplx( cos( twopi*x), sin( twopi*x), 8)
        end do      
#ifdef USEOMP                
!$OMP END DO
!$OMP END PARALLEL
#endif

        do ir = 1, nrpt
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst, ig, zfft)
!$OMP DO  
#endif
          do ist = wf_fst, wf_lst
            do jst = wf_fst, wf_lst
              do ig = 1, wf_Gset%ngrtot
                zfft( ig) = conjg( wfir( ig, ist))*wanfir( ig, jst, ir)*cfun( ig)
              end do
              call zfftifc( 3, ngrid, -1, zfft)
              wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
                  zfft( igfft( wf_Gset%ivgig( 0, 0, 0)))*wswgt( ir)
            end do
          end do
#ifdef USEOMP                
!$OMP END DO
!$OMP END PARALLEL
#endif
        end do

      end do   
    
      deallocate( wanfmt, wanfir, wfmt, wfir, evecfv, apwalm, match_combined, cfun, zfft, nlmo, lmo2l, lmo2m, lmo2o)

      return
    end subroutine wannier_projonwan
    
    subroutine wannier_gen_fromfile
      logical :: success
      real(8), allocatable :: ravg(:,:), omega(:)

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call wannier_readfile( success)
      if( .not. success) then
        write( wf_info, '(" Failed to read Wannier functions from file. Recalculate them.")')
        write( wf_info, *)
        call wannier_emat
        call wannier_gen_opf
        call wannier_maxloc
      else
        write( wf_info, '(" file successsfully read")')
        write( wf_info, *)
      end if
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      return
    end subroutine wannier_gen_fromfile
    
    subroutine wannier_subspace
      integer :: ik, ikb, i, nk, mk, n, is, no, ni, idxn, ist, jst, fst, lst, it, maxit
      real(8) :: mixing, maxdiff

      integer, allocatable :: idxo(:), cntb(:)
      real(8), allocatable :: sval(:), eval(:), evalmem(:,:)
      complex(8), allocatable :: projm(:,:), z(:,:,:), evec(:,:), lsvec(:,:), rsvec(:,:), auxmat(:,:), auxmat2(:,:), m(:,:,:,:), subspace(:,:,:), subspace_mem(:,:,:)

      write( *, '(2F13.6)') wf_win_i
      write( *, '(2F13.6)') wf_win_o
      write(*,*) wf_fst, wf_lst, wf_nst, wf_nwf
      write(*,*)
      !do ik = 1, wf_kset%nkpt
      !  write(*,*) wf_win_ni( ik), wf_win_ii( 1:wf_win_ni( ik), ik)
      !  write(*,*) wf_win_no( ik), wf_win_io( 1:wf_win_ni( ik), ik)
      !  write(*,*)
      !end do
      !stop

      allocate( idxo( wf_nwf))
      allocate( cntb( wf_fst:wf_lst))
      cntb = 0
      idxo = 0
      do ik = 1, wf_kset%nkpt
        if( wf_win_ni( ik) .gt. 0) then
          do i = 1, wf_win_ni( ik)
            cntb( wf_win_ii( i, ik)) = cntb( wf_win_ii( i, ik)) + 1
          end do
        else
          do i = 1, wf_win_no( ik)
            cntb( wf_win_io( i, ik)) = cntb( wf_win_io( i, ik)) + 1
          end do
        end if
      end do
      i = 0
      do while( (i .lt. wf_nwf) .and. (maxval( cntb) .gt. 0))
        n = wf_fst + maxloc( cntb, 1) - 1
        i = i + 1
        idxo( i) = n
        cntb( n) = 0
      end do
      do while( i .lt. wf_nwf)
        n = maxval( idxo) + 1
        i = i + 1
        idxo( i) = n
      end do
      write(*,'("used for initial OPF: ", 100i4)') idxo
        
      maxit = 10000

      ! search for nwf OPFs
      if( allocated( wf_m0)) deallocate( wf_m0)
      allocate( wf_m0( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))
      do n = 1, wf_nwf
        do i = 1, wf_nwf
          wf_m0( n, i, :, :) = wf_pwmat( idxo( n), idxo( i), :, :)
        end do
      end do

      call wannier_loc( initial=.true.)
      write(*,'(4F13.6)') sum( wf_omega), sum( wf_omega_i), sum( wf_omega_d), sum( wf_omega_od)
      wf_disentangle = .false.
      call wannier_gen_opf( nopf=wf_nwf, bands=idxo, nowrite=.true.)

      allocate( projm( wf_nst, wf_nprojused), m( wf_nst, wf_nst, wf_kset%nkpt, wf_n_ntot))
      allocate( lsvec( wf_nst, wf_nst), rsvec( wf_nwf, wf_nwf), sval( wf_nst))
      allocate( subspace( wf_nst, wf_nwf, wf_kset%nkpt))
      allocate( subspace_mem( wf_nst, wf_nwf, wf_kset%nkpt))
      allocate( auxmat( wf_nst, wf_nwf))
      allocate( auxmat2( wf_nst, wf_nst))
      m = zzero
      
      ! get initial subspace
      do ik = 1, wf_kset%nkpt
        nk = wf_win_no( ik) + wf_win_ni( ik)
        !write(*,*) ik, nk
        !write(*,'(I5," - ",100I4)') wf_win_ni( ik), wf_win_ii( 1:wf_win_ni( ik), ik)
        !write(*,'(I5," - ",100I4)') wf_win_no( ik), wf_win_io( 1:wf_win_no( ik), ik)
        projm( :, :) = zzero
        subspace( :, :, ik) = zzero
        jst = 0
        do i = 1, wf_nprojtot
          if( wf_projused( i) .eq. 1) then
            jst = jst + 1
            do n = 1, wf_win_ni( ik)
              projm( n, jst) = wf_projection( wf_win_ii( n, ik), i, ik)
            end do
            do n = 1, wf_win_no( ik)
              projm( wf_win_ni( ik)+n, jst) = wf_projection( wf_win_io( n, ik), i, ik)
            end do
          end if
        end do
        
        call zgemm( 'n', 'n', nk, wf_nwf, wf_nprojused, zone, &
               projm( 1:nk, :), nk, &
               wf_opf, wf_nprojused, zzero, & 
               auxmat( 1:nk, :), nk)
        call zgesdd_wrapper( auxmat( 1:nk, :), &
               nk, wf_nwf, &
               sval( 1:wf_nwf), &
               lsvec( 1:nk, 1:nk), &
               rsvec)
        call zgemm( 'n', 'n', nk, wf_nwf, wf_nwf, zone, &
               lsvec( 1:nk, 1:wf_nwf), nk, &
               rsvec, wf_nwf, zzero, &
               auxmat( 1:nk, :), nk)
        if( wf_win_ni( ik) .gt. 0) then
          call zgemm( 'n', 'c', nk, nk, wf_nwf, zone, &
                 auxmat( 1:nk, :), nk, &
                 auxmat( 1:nk, :), nk, zzero, &
                 auxmat2( 1:nk, 1:nk), nk)
          mk = wf_nwf - wf_win_ni( ik)
          auxmat2( 1:wf_win_ni( ik), :) = zzero
          auxmat2( :, 1:wf_win_ni( ik)) = zzero
          call diaghermat( nk, auxmat2( 1:nk, 1:nk), sval( 1:nk), lsvec( 1:nk, 1:nk))
          do n = 1, mk
            !write(*,'("## ",I,F13.6)') wf_nwf-n, sval( nk-n)
            subspace( 1:wf_win_no( ik), n, ik) = lsvec( (wf_win_ni( ik)+1):nk, nk-n+1)
          end do
          !write(*,*) wf_win_ni( ik), nk, mk
        else
          subspace( 1:nk, :, ik) = auxmat( 1:nk, :)
        end if
        !call plotmat( subspace( :, :, ik))
      end do

      m = zzero
      do ik = 1, wf_kset%nkpt
        do i = 1, wf_win_no( ik)
          do n = 1, wf_win_no( ik)
            m( i, n, ik, :) = wf_pwmat( wf_win_io( i, ik), wf_win_io( n, ik), ik, :)
          end do
        end do
      end do
 
      allocate( evec( wf_nst, wf_nst), eval( wf_nst), evalmem( wf_nst, wf_kset%nkpt))
      allocate( z( wf_nst, wf_nst, wf_kset%nkpt))
      z = zzero
      maxdiff = 1.d0
      evalmem = 0.d0
      it = 0
      write(*,*) "start disentanglement"
      do while( (maxdiff .gt. 1.d-3) .and. (it .lt. maxit))
        it = it + 1
        maxdiff = 0.d0
        mixing = 0.5d0
        if( it .eq. 1) mixing = 1.d0
#ifdef USEOMP
!!$omp parallel default( shared) private( ik, i, is, n, idxn, ikb, eval, evec, auxmat)
!!$omp do
#endif
        do ik = 1, wf_kset%nkpt
          nk = wf_win_no( ik)
          z( :, :, ik) = (1.d0 - mixing)*z( :, :, ik)
          if( nk .ge. 1) then
            do i = 1, wf_n_nshells
              is = wf_n_usedshells( i)
              do n = 1, wf_n_n( is)
                idxn = wf_n_ns2n( n, is)
                ikb = wf_n_ik( n, is, ik)
                mk = wf_win_no( ikb)
                if( mk .ge. 1) then
                  ni = wf_nwf - wf_win_ni( ikb)
                  call zgemm( 'n', 'n', nk, ni, mk, zone, &
                         m( 1:nk, 1:mk, ik, idxn), nk, &
                         subspace( 1:mk, 1:ni, ikb), mk, zzero, &
                         auxmat( 1:nk, 1:ni), nk)
                  call zgemm( 'n', 'c', nk, nk, ni, cmplx( mixing*wf_n_wgt( is), 0, 8), &
                         auxmat( 1:nk, 1:ni), nk, &
                         auxmat( 1:nk, 1:ni), nk, zone, &
                         z( 1:nk, 1:nk, ik), nk)
                end if

                ikb = wf_n_ik2( n, is, ik)
                mk = wf_win_no( ikb)
                if( mk .ge. 1) then
                  ni = wf_nwf - wf_win_ni( ikb)
                  call zgemm( 'c', 'n', nk, ni, mk, zone, &
                         m( 1:mk, 1:nk, ikb, idxn), mk, &
                         subspace( 1:mk, 1:ni, ikb), mk, zzero, &
                         auxmat( 1:nk, 1:ni), nk)
                  call zgemm( 'n', 'c', nk, nk, ni, cmplx( mixing*wf_n_wgt( is), 0, 8), &
                         auxmat( 1:nk, 1:ni), nk, &
                         auxmat( 1:nk, 1:ni), nk, zone, &
                         z( 1:nk, 1:nk, ik), nk)
                end if
              end do
            end do
          end if
          ni = wf_nwf - wf_win_ni( ik)
          if( ni .ge. 1) then
            call diaghermat( nk, z( 1:nk, 1:nk, ik), eval( 1:nk), evec( 1:nk, 1:nk))
            subspace( 1:nk, 1:ni, ik) = evec( 1:nk, (nk-ni+1):nk)
          end if
          !subspace_mem( :, :, ik) = evec( :, (no-ni+1):no)
          !write(*,'(100F13.6)') eval
#ifdef USEOMP
!!$omp atomic update
#endif
          maxdiff = max( maxdiff, norm2( eval( (nk-ni+1):nk) - evalmem( (nk-ni+1):nk, ik)))
#ifdef USEOMP
!!$omp end atomic
#endif
          evalmem( 1:nk, ik) = eval( 1:nk)
        end do
#ifdef USEOMP
!!$omp end do
!!$omp end parallel
#endif
        !subspace = subspace_mem
        if( mod( it, 100) .eq. 0) write(*,*) it, maxdiff
      end do

      !do ik = 1, wf_kset%nkpt
      !  write(*,*) ik
      !  write(*,'(I5," - ",100I4)') wf_win_ni( ik), wf_win_ii( 1:wf_win_ni( ik), ik)
      !  write(*,'(I5," - ",100I4)') wf_win_no( ik), wf_win_io( 1:wf_win_no( ik), ik)
      !  call plotmat( subspace( :, :, ik))
      !end do
      !stop

      if( allocated( wf_subspace)) deallocate( wf_subspace)
      allocate( wf_subspace( wf_fst:wf_lst, wf_nwf, wf_kset%nkpt))
      wf_subspace = zzero
      do ik = 1, wf_kset%nkpt
        do n = 1, wf_win_ni( ik)
          wf_subspace( wf_win_ii( n, ik), n, ik) = zone
        end do
        do n = 1, wf_win_no( ik)
          do i = 1, wf_nwf - wf_win_ni( ik)
            wf_subspace( wf_win_io( n, ik), wf_win_ni( ik)+i, ik) = subspace( n, i, ik)
          end do
        end do
        !write(*,*) ik
        !write(*,'(I5," - ",100I4)') wf_win_ni( ik), wf_win_ii( 1:wf_win_ni( ik), ik)
        !write(*,'(I5," - ",100I4)') wf_win_no( ik), wf_win_io( 1:wf_win_no( ik), ik)
        !call plotmat( wf_subspace( :, :, ik))
      end do

      call wannier_diagonalize_subspace

      if( allocated( wf_m0)) deallocate( wf_m0)
      allocate( wf_m0( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))
      deallocate( auxmat)
      allocate( auxmat( wf_fst:wf_lst, wf_nwf))
      do i = 1, wf_n_nshells
        is = wf_n_usedshells( i)
        do n = 1, wf_n_n( is)
          idxn = wf_n_ns2n( n, is)
          do ik = 1, wf_kset%nkpt
            ikb = wf_n_ik( n, is, ik)
            call zgemm( 'n', 'n', wf_nst, wf_nwf, wf_nst, zone, &
                 wf_pwmat( :, :, ik, idxn), wf_nst, &
                 wf_subspace( :, :, ikb), wf_nst, zzero, &
                 auxmat, wf_nst)
            call zgemm( 'c', 'n', wf_nwf, wf_nwf, wf_nst, zone, &
                 wf_subspace( :, :, ik), wf_nst, &
                 auxmat, wf_nst, zzero, &
                 wf_m0( :, :, ik, idxn), wf_nwf)
          end do
        end do
      end do

      deallocate( auxmat, m, subspace, evec, evalmem, eval, lsvec, rsvec, sval, idxo, projm, z, subspace_mem)
      
      write(*,*) it
      ! deallocate them since they have the wrong shape
      if( allocated( wf_m)) deallocate( wf_m)
      if( allocated( wf_centers)) deallocate( wf_centers)
      if( allocated( wf_omega)) deallocate( wf_omega)
      if( allocated( wf_omega_i)) deallocate( wf_omega_i)
      if( allocated( wf_omega_d)) deallocate( wf_omega_d)
      if( allocated( wf_omega_od)) deallocate( wf_omega_od)
      call wannier_loc( initial=.true.)
      write(*,'(4F13.6)') sum( wf_omega), sum( wf_omega_i), sum( wf_omega_d), sum( wf_omega_od)

      wf_disentangle = .true.

      call wannier_gen_opf( maxit=1)
      call wannier_maxloc

      !fst = wf_fst
      !lst = wf_fst + wf_nwf - 1
      !wf_lst = lst
      !wf_nst = wf_nwf
      !wf_m0 = wf_pwmat( fst:lst, fst:lst, :, :)
      !
      !allocate( auxmat( wf_nwf, wf_nwf))
      !allocate( projm( wf_nwf, wf_nwf))
      !allocate( lsvec( wf_nwf, wf_nwf))
      !allocate( rsvec( wf_nwf, wf_nwf))
      !allocate( sval( wf_nwf))
      !do ik = 1, wf_kset%nkpt
      !  !call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
      !  !     wf_sub_transform( :, :, ik), wf_nwf, &
      !  !     wf_transform( :, :, ik), wf_nwf, zzero, &
      !  !     auxmat, wf_nwf)
      !  call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
      !       wf_subspace( fst:lst, :, ik), wf_nwf, &
      !       wf_transform( :, :, ik), wf_nwf, zzero, &
      !       projm, wf_nwf)
      !  call zgesdd_wrapper( projm, wf_nwf, wf_nwf, sval, lsvec, rsvec)
      !  call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
      !       lsvec, wf_nwf, &
      !       rsvec, wf_nwf, zzero, &
      !       wf_transform( :, :, ik), wf_nwf)
      !end do
      !
      !call wannier_loc
      !call wannier_maxloc
    end subroutine wannier_subspace

    subroutine wannier_diagonalize_subspace
      use m_getunit

      integer :: ik, iq, ist, nkpqp, fstqp, lstqp, un, recl, nk, nkequi, isymequi( wf_kset%nkpt), ikequi( wf_kset%nkpt), fst, lst
      real(8) :: eval( wf_fst:wf_lst, wf_kset%nkpt), vl(3), efermiqp, efermiks
      character(256) :: fname, fxt
      complex(8) :: hamilton( wf_nwf, wf_nwf), auxmat( wf_fst:wf_lst, wf_nwf), evec( wf_nwf, wf_nwf)
      logical :: exist

      real(8), allocatable :: evalfv(:,:)

      call wannier_geteval( evalfv, fst, lst)

      if( allocated( wf_sub_eval)) deallocate( wf_sub_eval)
      allocate( wf_sub_eval( wf_nwf, wf_kset%nkpt))
      if( allocated( wf_sub_transform)) deallocate( wf_sub_transform)
      allocate( wf_sub_transform( wf_nwf, wf_nwf, wf_kset%nkpt))

      do ik = 1, wf_kset%nkpt
        do ist = wf_fst, wf_lst
          auxmat( ist, :) = wf_subspace( ist, :, ik)*evalfv( ist, ik)
        end do
        call zgemm( 'c', 'n', wf_nwf, wf_nwf, wf_nst, zone, &
               wf_subspace( :, :, ik), wf_nst, &
               auxmat, wf_nst, zzero, &
               hamilton, wf_nwf)
        call diaghermat( wf_nwf, hamilton, wf_sub_eval( :, ik), wf_sub_transform( :, :, ik))
        call zgemm( 'n', 'n', wf_nst, wf_nwf, wf_nwf, zone, &
               wf_subspace( :, :, ik), wf_nst, &
               wf_sub_transform( :, :, ik), wf_nwf, zzero, &
               auxmat, wf_nst)
        wf_subspace( :, :, ik) = auxmat
        write(*,*) ik
        write(*,'(100F13.6)') wf_sub_eval( :, ik)
        write(*,'(100F13.6)') evalfv( wf_fst:wf_lst, ik)
        write(*,'(I5," - ",100I4)') wf_win_ni( ik), wf_win_ii( 1:wf_win_ni( ik), ik)
        write(*,'(I5," - ",100I4)') wf_win_no( ik), wf_win_io( 1:wf_win_no( ik), ik)
        call plotmat( wf_subspace( :, :, ik))
        write(*,*)
      end do

    end subroutine wannier_diagonalize_subspace

    subroutine wannier_loc( totonly, initial)
      use mod_lattice, only : ainv
      logical, optional, intent( in) :: totonly, initial

      integer :: iknr, is, n, i, j, k, idxn
      complex(8), allocatable :: auxmat(:,:)
      real(8), allocatable :: logsum(:,:), log2sum(:,:), abssum(:,:), abs2sum(:,:)
      real(8) :: tmp, o, oi, od, ood, v3(3)
      logical :: tot, init

      tot = .false.
      if( present( totonly)) tot = totonly
      init = .false.
      if( present( initial)) init = initial
      
      allocate( auxmat( wf_nwf, wf_nwf))
      allocate( logsum( wf_nwf, wf_n_ntot))
      allocate( log2sum( wf_nwf, wf_n_ntot))
      allocate( abssum( wf_nwf, wf_n_ntot))
      allocate( abs2sum( wf_nwf, wf_n_ntot))

      if( .not. wf_initialized) call wannier_init
      if( .not. allocated( wf_m0)) then
        write(*,*) "calc emat"
        call wannier_emat
      end if
      if( .not. allocated( wf_m)) allocate( wf_m( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))
      if( .not. allocated( wf_sheet)) allocate( wf_sheet( wf_nwf, wf_kset%nkpt, wf_n_ntot))
      if( .not. allocated( wf_centers)) then
        allocate( wf_centers( 3, wf_nwf))
        wf_centers = 0.d0
      end if
      if( .not. allocated( wf_omega)) allocate( wf_omega( wf_nwf))
      if( .not. tot) then
        if( .not. allocated( wf_omega_i)) allocate( wf_omega_i( wf_nwf))
        if( .not. allocated( wf_omega_d)) allocate( wf_omega_d( wf_nwf))
        if( .not. allocated( wf_omega_od)) allocate( wf_omega_od( wf_nwf))
      end if

      wf_sheet = 0.d0
      !do iknr = 1, wf_kset%nkpt
      !  do i = 1, wf_n_nshells 
      !    is = wf_n_usedshells( i)
      !    do n = 1, wf_n_n( is)
      !      idxn = wf_n_ns2n( n, is)
      !      k = wf_n_ik( n, is, iknr)
      !      do j = 1, wf_nwf
      !        wf_sheet( j, iknr, idxn) = dot_product( wf_kset%vkc( :, k), wf_centers( :, j))
      !        !write(*,*) iknr, idxn, j, wf_sheet( j, iknr, idxn)
      !      end do
      !    end do
      !  end do
      !end do
      if( .not. init) then
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, is, i, n, idxn, auxmat)
!$OMP DO
#endif
        do iknr = 1, wf_kset%nkpt
          do i = 1, wf_n_nshells 
            is = wf_n_usedshells( i)
            do n = 1, wf_n_n( is)
              idxn = wf_n_ns2n( n, is)
              call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
                   wf_m0( :, :, iknr, idxn), wf_nwf, &
                   wf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nwf, zzero, &
                   auxmat, wf_nwf)
              call ZGEMM( 'c', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
                   wf_transform( :, :, iknr), wf_nwf, &
                   auxmat, wf_nwf, zzero, &
                   wf_m( :, :, iknr, idxn), wf_nwf)
            end do
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      else
        wf_m = wf_m0
        write(*,*) "wannier_loc: init"
      end if

      logsum = 0.d0
      log2sum = 0.d0
      abssum = 0.d0
      abs2sum = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( j, i, n, idxn, iknr, tmp, is)
!$OMP DO COLLAPSE(2)
#endif
      !do i = 1, wf_n_nshells 
      !  is = wf_n_usedshells( i)
      !  do n = 1, wf_n_n( is)
      !    idxn = wf_n_ns2n( n, is)
        do idxn = 1, wf_n_ntot
          do j = 1, wf_nwf
            do iknr = 1, wf_kset%nkpt
              tmp = 0.d0
              if( abs( wf_m( j, j, iknr, idxn)) .gt. 1.d-10) then
                !tmp = aimag( log( wf_m( j, j, iknr, idxn)))
                tmp = aimag( log( exp( zi*wf_sheet( j, iknr, idxn))*wf_m( j, j, iknr, idxn))) - wf_sheet( j, iknr, idxn)
                !if( (tmp + dot_product( wf_n_vc( :, n, is), wf_centers( :, j))) .gt. pi) tmp = tmp - twopi
                !if( (tmp + dot_product( wf_n_vc( :, n, is), wf_centers( :, j))) .lt. -pi) tmp = tmp + twopi
                !tmp = aimag( log( exp( zi*wf_sheet( j, iknr, idxn))*wf_m( j, j, iknr, idxn))) - wf_sheet( j, iknr, idxn)
                logsum( j, idxn) = logsum( j, idxn) + tmp
                log2sum( j, idxn) = log2sum( j, idxn) + tmp*tmp
                abssum( j, idxn) = abssum( j, idxn) + dble( wf_m( j, j, iknr, idxn)*conjg( wf_m( j, j, iknr, idxn)))
              end if
              if( .not. tot) then
                do k = 1, wf_nwf
                  abs2sum( j, idxn) = abs2sum( j, idxn) + 0.5d0*( dble( wf_m( j, k, iknr, idxn)*conjg( wf_m( j, k, iknr, idxn))) + &
                                                                     dble( wf_m( k, j, iknr, idxn)*conjg( wf_m( k, j, iknr, idxn))))
                end do
              end if
            end do
          end do
        end do
      !end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      logsum = logsum/wf_kset%nkpt
      log2sum = log2sum/wf_kset%nkpt
      abssum = abssum/wf_kset%nkpt
      abs2sum = abs2sum/wf_kset%nkpt

      wf_centers = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( i, n, j, is, idxn)
!$OMP DO
#endif
      do j = 1, wf_nwf
        do i = 1, wf_n_nshells 
          is = wf_n_usedshells( i)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            wf_centers( :, j) = wf_centers( :, j) - 2.d0*wf_n_wgt( is)*logsum( j, idxn)*wf_n_vc( :, n, is)
          end do
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      !if( .not. tot) then
      !  do j = 1, wf_nwf
      !    call r3mv( ainv, wf_centers( :, j), v3)
      !    write(*,'(3F13.5)') v3
      !  end do
      !end if

      wf_omega = 0.d0
      wf_omega_i = 0.d0
      wf_omega_d = 0.d0
      wf_omega_od = 0.d0

      if( tot) then
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( i, is, n, idxn, j)
!$OMP DO
#endif
        do j = 1, wf_nwf
          do i = 1, wf_n_nshells 
            do n = 1, wf_n_n( wf_n_usedshells( i))
              is = wf_n_usedshells( i)
              idxn = wf_n_ns2n( n, is)
              ! total spread
              wf_omega( j) = wf_omega( j) + 2.d0*wf_n_wgt( is)*( 1.d0 - abssum( j, idxn) + log2sum( j, idxn))
            end do
          end do
          wf_omega( j) = wf_omega( j) - dot_product( wf_centers( :, j), wf_centers( :, j))
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      else
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( i, is, n, idxn, j)
!$OMP DO
#endif
        do j = 1, wf_nwf
          do i = 1, wf_n_nshells 
            do n = 1, wf_n_n( wf_n_usedshells( i))
              is = wf_n_usedshells( i)
              idxn = wf_n_ns2n( n, is)
              ! total spread
              wf_omega( j) = wf_omega( j) + 2.d0*wf_n_wgt( is)*( 1.d0 - abssum( j, idxn) + log2sum( j, idxn))
              ! gauge independent spread
              wf_omega_i( j) = wf_omega_i( j) + 2.d0*wf_n_wgt( is)*( 1.d0 - abs2sum( j, idxn))
              ! diagonal spread
              wf_omega_d( j) = wf_omega_d( j) + 2.d0*wf_n_wgt( is)*log2sum( j, idxn)
              ! off-diagonal spread
              wf_omega_od( j) = wf_omega_od( j) + 2.d0*wf_n_wgt( is)*( abs2sum( j, idxn) - abssum( j, idxn))
            end do
          end do
          wf_omega( j) = wf_omega( j) - dot_product( wf_centers( :, j), wf_centers( :, j))
          wf_omega_d( j) = wf_omega_d( j) - dot_product( wf_centers( :, j), wf_centers( :, j))
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      end if

      deallocate( auxmat, logsum, log2sum, abssum, abs2sum)
      return
    end subroutine wannier_loc

    subroutine wannier_phases( firstcall)
      logical, optional, intent( in) :: firstcall

      integer :: ik, ist, jst, idxn, i, n, s, i1, i2
      real(8) :: smat(3,3), svec(3), sinv(3,3), xx( wf_n_ntot), xx0
      complex(8) :: csum( wf_n_ntot)
      logical :: initial

      real(8) :: r3mdet

      initial = .false.
      if( present( firstcall)) initial = firstcall

      if( .not. allocated( wf_rguide)) allocate( wf_rguide( 3, wf_nwf))
      if( .not. allocated( wf_sheet)) allocate( wf_sheet( wf_nwf, wf_kset%nkpt, wf_n_ntot))

      do ist = 1, wf_nwf
        do idxn = 1, wf_n_ntot
          csum( idxn) = zzero
          do ik = 1, wf_kset%nkpt
            csum( idxn) = csum( idxn) + wf_m( ist, ist, ik, idxn)
          end do
        end do
        smat = 0.d0
        svec = 0.d0
        do i = 1, wf_n_nshells 
          s = wf_n_usedshells( i)
          do n = 1, wf_n_n( s)
            idxn = wf_n_ns2n( n, s)
            if( idxn .le. 3) then
              xx( idxn) = -aimag( log( csum( idxn)))
            else
              xx0 = dot_product( wf_n_vc( :, n, s), wf_rguide( :, ist))
              xx( idxn) = xx0 - aimag( log( csum( idxn)*exp( zi*xx0)))
            end if
            do i1 = 1, 3
              do i2 = 1, 3
                smat( i1, i2) = smat( i1, i2) + wf_n_vc( i1, n, s)*wf_n_vc( i2, n, s)
              end do
              svec( i1) = svec( i1) + wf_n_vc( i1, n, s)*xx( idxn)
            end do
            if( idxn .ge. 3) then
              if( (.not. initial) .and. (abs( r3mdet( smat)) .gt. 1.d-6)) then
                call r3minv( smat, sinv)
                do i1 = 1, 3
                  wf_rguide( i1, ist) = 0.d0
                  do i2 = 1, 3
                    wf_rguide( i1, ist) = wf_rguide( i1, ist) + sinv( i1, i2)*svec( i2)
                  end do
                end do
              end if
            end if
          end do
        end do
      end do

      wf_sheet = 0.d0
      do ik = 1, wf_kset%nkpt
        do i = 1, wf_n_nshells 
          s = wf_n_usedshells( i)
          do n = 1, wf_n_n( s)
            idxn = wf_n_ns2n( n, s)
            do ist = 1, wf_nwf
              do i1 = 1, 3
                wf_sheet( ist, ik, idxn) = wf_sheet( ist, ik, idxn) + wf_n_vc( i1, n, s)*wf_rguide( i1, ist)
              end do
            end do
          end do
        end do
      end do

      return
    end subroutine wannier_phases

    ! print out projection local-orbitals
    subroutine wannier_showproj
      ! local variables
      integer :: iproj

      write(*,*)
      write(*,*) 'Local-orbitals for generating Wannier functions via projection.'
      write(*,*) '     nr   species   atom   l    m    used'
      do iproj = 1, wf_nprojtot
        write( *, '(5x,i3,3x)', advance='no') iproj
        write( *, '(2x,a5,3x)', advance='no') spsymb( wf_projst( iproj, 1))
        write( *, '(i4,3x)', advance='no') wf_projst( iproj, 2)
        write( *, '(i2,3x)', advance='no') wf_projst( iproj, 4)
        write( *, '(i2,3x)', advance='no') wf_projst( iproj, 5)
        if( wf_projused( iproj)) write( *, '(1x,a)', advance='no') '*'
        write(*,*)
      end do
      write(*,*)
      return
    end subroutine wannier_showproj

!*********************************!
!          FILE HANDLING          !
!*********************************!
    
    ! reads transformation matrices from file
    subroutine wannier_readfile( success)
      logical, intent( out) :: success

      ! local variables
      integer :: ik, ix, iy, iz, un
      integer :: fst_, lst_, nst_, nwf_, nkpt_
      real(8) :: vkl_( 3), vkl_tmp( 3, wf_kset%nkpt)
      logical :: disentangle_

      call getunit( un)

      success = .true.
      inquire( file=trim( wf_filename)//"_TRANSFORM"//trim( filext), exist=success)
      if( .not. success) then
        write(*,*) 'ERROR (wannier_readfile): File '//trim( wf_filename)//"_TRANSFORM"//trim( filext)//' does not exist.'
        return
      end if
      open( un, file=trim( wf_filename)//"_TRANSFORM"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
      read( un) fst_, lst_, nst_, nwf_, nkpt_, disentangle_
      wf_disentangle = disentangle_
      if( (fst_ .ne. wf_fst) .or. (lst_ .ne. wf_lst)) then
        write( *, '(" ERROR (wannier_readfile): different band-ranges in input (",I4,":",I4,") and file (",I4,":",I4,").")'), wf_fst, wf_lst, fst_, lst_
        call terminate
      end if
      if( nwf_ .ne. wf_nwf) then
        write( *, '(" ERROR (wannier_readfile): different number of Wannier functions in input (",I4,") and file (",I4,").")') wf_nwf, nwf_
        call terminate
      end if
      if( nkpt_ .ne. wf_kset%nkpt) then
        write( *, '(" ERROR (wannier_readfile): different number of k-points in input (",I4,") and file (",I4,").")'), wf_kset%nkpt, nkpt_
        call terminate
      end if
      if( allocated( wf_omega)) deallocate( wf_omega)
      allocate( wf_omega( wf_nwf))
      if( allocated( wf_omega_i)) deallocate( wf_omega_i)
      allocate( wf_omega_i( wf_nwf))
      if( allocated( wf_omega_d)) deallocate( wf_omega_d)
      allocate( wf_omega_d( wf_nwf))
      if( allocated( wf_omega_od)) deallocate( wf_omega_od)
      allocate( wf_omega_od( wf_nwf))
      if( allocated( wf_omega)) deallocate( wf_omega)
      allocate( wf_omega( wf_nwf))
      if( allocated( wf_centers)) deallocate( wf_centers)
      allocate( wf_centers( 3, wf_nwf))
      if( allocated( wf_subspace)) deallocate( wf_subspace)
      allocate( wf_transform( wf_nwf, wf_nwf, wf_kset%nkpt))
      do ik = 1, wf_kset%nkpt
        read( un) vkl_
        vkl_tmp( 1, :) = wf_kset%vkl( 1, :) - vkl_( 1)
        vkl_tmp( 2, :) = wf_kset%vkl( 2, :) - vkl_( 2)
        vkl_tmp( 3, :) = wf_kset%vkl( 3, :) - vkl_( 3)
        iz = minloc( norm2( vkl_tmp( :, :), 1), 1)
        if( norm2( vkl_tmp( :, iz)) .gt. input%structure%epslat) then
          write( *, '(" ERROR (wannier_readfile): k-point in file not in k-point-set.")')
          write( *, '(3F23.6)') vkl_
          call terminate
        end if
        do iy = 1, wf_nwf
          do ix = 1, wf_nwf
            read( un) wf_transform( ix, iy, iz)
          end do
        end do
      end do
      do ix = 1, wf_nprojtot
        read( un) wf_projused( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega_i( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega_d( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega_od( ix)
      end do
      do iy = 1, wf_nwf
        do ix = 1, 3
          read( un) wf_centers( ix, iy)
        end do
      end do
      close( un)
      if( success) write(*,*) 'Transformation matrices successfully read.'
      return
    end subroutine wannier_readfile
    
    ! writes transformation matrices to file
    subroutine wannier_writefile
      ! local variables
      integer :: ik, ix, iy, un
      
      call getunit( un)

      open( un, file=trim( wf_filename)//"_TRANSFORM"//trim( filext), action='WRITE', form='UNFORMATTED')
      write( un) wf_fst, wf_lst, wf_nst, wf_nwf, wf_kset%nkpt, wf_disentangle
      do ik = 1, wf_kset%nkpt
        write( un) wf_kset%vkl( :, ik)
        do iy = 1, wf_nwf
          do ix = 1, wf_nwf
            write( un) wf_transform( ix, iy, ik)
          end do
        end do
      end do
      do ix = 1, wf_nprojtot
        write( un) wf_projused( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega_i( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega_d( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega_od( ix)
      end do
      do iy = 1, wf_nwf
        do ix = 1, 3
          write( un) wf_centers( ix, iy)
        end do
      end do
      close( un)
      write( *, '(a,a)') ' Transformation matrices written to file ', trim( wf_filename)//"_TRANSFORM"//trim( filext)
      return
    end subroutine wannier_writefile  

    subroutine wannier_writefun( nshell)
      !!!!!
      ! correct wf_nst and wf_nwf here
      !!!!!
      use m_wsweight
      integer, intent( in) :: nshell

      integer :: ik, ist, jst, is, ia, ias, l, m, o, lm, lmo, ilo, ig, igk, ifg, ir
      integer :: lmaxapw, nlmomax, ngknr, nrpt
      integer :: un, recl, offset
      real(8) :: x

      integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:)
      complex(8), allocatable :: wanfmt(:,:,:,:), wanfir(:,:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), match_combined(:,:,:), wfmt(:,:), wfir(:,:)
      complex(8), allocatable :: auxmat(:,:), auxmat2(:,:), wswgt(:)
      real(8), allocatable :: rptl(:,:)

      nrpt = (1 + 2*nshell)**3

      allocate( rptl( 3, nrpt))
      ir = 0
      do l = -nshell, nshell
        do m = -nshell, nshell
          do o = -nshell, nshell
            ir = ir + 1
            rptl( :, ir) = dble( (/o, m, l/))
          end do
        end do
      end do

      lmaxapw = input%groundstate%lmaxapw
      ! count combined (l,m,o) indices and build index maps
      allocate( nlmo( nspecies))
      allocate( lmo2l( (lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2m( (lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2o( (lmaxapw + 1)**2*apwordmax, nspecies))
      nlmomax = 0
      do is = 1, nspecies
        nlmo( is) = 0
        do l = 0, lmaxapw
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
      if( nlmomax .ne. lmmaxapw*apwordmax) write(*,*) "ERROR (wannier_writefun): wrong nlmomax"

      call readstate
      call readfermi
      call linengy
      call genapwfr
      call genlofr
      call olprad
      call genidxlo

      allocate( wanfmt( nlmomax+nlotot, wf_fst:wf_lst, natmtot, nrpt))
      allocate( wanfir( wf_Gset%ngrtot, wf_fst:wf_lst, nrpt))

      wanfmt = zzero
      wanfir = zzero

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      allocate( match_combined( nlmomax, ngkmax, natmtot))
      allocate( wfmt( nlmomax+nlotot, wf_fst:wf_lst))
      allocate( wfir( wf_Gset%ngrtot, wf_fst:wf_lst))
      allocate( auxmat( nlmomax+nlotot, wf_fst:wf_lst))
      allocate( auxmat2( wf_Gset%ngrtot, wf_fst:wf_lst))
      allocate( wswgt( nrpt))

      do ik = 1, wf_kset%nkpt
        ngknr = wf_Gkset%ngk( 1, ik)

        do ir = 1, nrpt
          call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), wswgt( ir), kgrid=.true.)
        end do

        ! get matching coefficients
        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
          
        ! read eigenvector      
        if( input%properties%wannier%input .eq. "gs") then
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
            do lmo = 1, nlmo( is)
              l = lmo2l( lmo, is)
              m = lmo2m( lmo, is)
              o = lmo2o( lmo, is)
              lm = idxlm( l, m)
              match_combined( lmo, 1:ngknr, ias) = apwalm( 1:ngknr, o, lm, ias, 1)
            end do
          end do
        end do

        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            wfmt = zzero
            call zgemm( 'N', 'N', nlmo( is), wf_nst, ngknr, zone, &
                 match_combined( 1:nlmo( is), 1:ngknr, ias), nlmo( is), &
                 evecfv( 1:ngknr, wf_fst:wf_lst, 1), ngknr, zone, &
                 wfmt( 1:nlmo( is), :), nlmo( is))
            do ilo = 1, nlorb( is)
              l = lorbl( ilo, is)
              do m = -l, l
                lm = idxlm( l, m)
                wfmt( nlmo( is)+idxlo( lm, ilo, ias), :) = wfmt( nlmo( is)+idxlo( lm, ilo, ias), :) + evecfv( ngknr+idxlo( lm, ilo, ias), wf_fst:wf_lst, 1)
              end do
            end do

            call zgemm( 'N', 'N', nlmomax+nlotot, wf_nst, wf_nst, zone, &
                 wfmt, nlmomax+nlotot, &
                 wf_transform( :, :, ik), wf_nst, zzero, &
                 auxmat, nlmomax+nlotot)

            do ir = 1, nrpt
              wanfmt( :, :, ias, ir) = wanfmt( :, :, ias, ir) + conjg( wswgt( ir))/wf_kset%nkpt*auxmat(:,:)
            end do

          end do
        end do
      
        wfir = zzero
        do igk = 1, wf_Gkset%ngk( 1, ik)
          ifg = igfft( wf_Gkset%igkig( igk, 1, ik))
          wfir( ifg, :) = evecfv( igk, wf_fst:wf_lst, 1)
        end do

        do ist = wf_fst, wf_lst
          call zfftifc( 3, ngrid, 1, wfir( :, ist))
        end do

        do ig = 1, wf_Gset%ngrtot
          x = wf_kset%vkl( 1, ik)*wf_Gset%ivg( 1, ig)/ngrid(1) + &
              wf_kset%vkl( 2, ik)*wf_Gset%ivg( 2, ig)/ngrid(2) + &
              wf_kset%vkl( 3, ik)*wf_Gset%ivg( 3, ig)/ngrid(3)
          ifg = igfft( ig)
          wfir( ifg, :) = wfir( ifg, :)*cmplx( cos( twopi*x), sin( twopi*x), 8)
        end do

        call zgemm( 'N', 'N', wf_Gset%ngrtot, wf_nst, wf_nst, conjg( wswgt)/wf_kset%nkpt, &
             wfir, wf_Gset%ngrtot, &
             wf_transform( :, :, ik), wf_nst, zzero, &
             auxmat2, wf_Gset%ngrtot)

        do ir = 1, nrpt
          wanfir( :, :, ir) = wanfir( :, :, ir) + wswgt( ir)/wf_kset%nkpt*auxmat2(:,:)
        end do

      end do

      deallocate( apwalm, match_combined, evecfv, wfir, wswgt)

      call getunit( un)
      open( un, file=trim( wf_filename)//"_FUN"//trim( filext), action='write', form='unformatted')
      ! constants
      write( un) wf_fst, wf_lst, lmaxapw, apwordmax, nlmomax, nlotot, lolmmax, nlomax, natmtot, nrmtmax, wf_Gset%ngrtot, ngrid
      ! maps
      write( un) nlmo
      write( un) lmo2l
      write( un) lmo2m
      write( un) lmo2o
      write( un) idxlo
      ! Wannier functions
      write( un) wanfmt
      write( un) wanfir
          
      close( un)
      write(*,*) "WANNIER FUNCTIONS WRITTEN"
      
      deallocate( wanfmt, wanfir, nlmo, lmo2l, lmo2m, lmo2o)

      return
    end subroutine wannier_writefun

    subroutine wannier_readfun( nshell, wanfmt, wanfir, wnlmo, wlmo2l, wlmo2m, wlmo2o)
      integer, intent( in) :: nshell
      complex(8), intent( out) :: wanfmt( lmmaxapw*apwordmax+nlotot, wf_fst:wf_lst, natmtot, (1+2*nshell)**3)
      complex(8), intent( out) :: wanfir( wf_Gset%ngrtot, wf_fst:wf_lst, (1+2*nshell)**3)
      integer, intent( out) :: wnlmo( nspecies)
      integer, intent( out) :: wlmo2l( lmmaxapw*apwordmax, nspecies)
      integer, intent( out) :: wlmo2m( lmmaxapw*apwordmax, nspecies)
      integer, intent( out) :: wlmo2o( lmmaxapw*apwordmax, nspecies)

      integer :: i, is, ia, ias, l, o, ilo, un, wf_fst_, wf_lst_, lmaxapw_, apwordmax_, nlmomax_, nlotot_, lolmmax_, nlomax_, natmtot_, nrmtmax_, ngrtot_, ngrid_(3), idxlo_( lolmmax, nlomax, natmtot), nrpt
      logical :: exist
    
      nrpt = (1 + 2*nshell)**3

      call getunit( un)
    
      do i = 1, 100
        inquire( file=trim( wf_filename)//"_FUN"//trim( filext), exist=exist)
        if( exist) then
          open( un, file=trim( wf_filename)//"_FUN"//trim( filext), action='read', form='unformatted')
          exit
        else
          call system( 'sync')
          write(*,*) "Waiting for other process to write"
          call sleep( 1)
        end if
      end do
    
      read( un) wf_fst_, wf_lst_, lmaxapw_, apwordmax_, nlmomax_, nlotot_, lolmmax_, nlomax_, natmtot_, nrmtmax_, ngrtot_, ngrid_
      if( (wf_fst_ .ne. wf_fst) .or. (wf_lst_ .ne. wf_lst)) then
        write(*, '("Error (wannier_readfun): invalid band ranges")')
        write (*, '(" current	   : ", 2I8)') wf_fst, wf_lst
        write (*, '(" in file      : ", 2I8)') wf_fst_, wf_lst_
        stop
      end if
      if( nrmtmax_ .ne. nrmtmax) then
        write(*, '("Error (wannier_readfun): invalid number of radial points")')
        write (*, '(" current	   : ", I8)') nrmtmax
        write (*, '(" in file      : ", I8)') nrmtmax_
        stop
      end if
      if( lmaxapw_ .ne. input%groundstate%lmaxapw) then
        write(*, '("Error (wannier_readfun): invalid maximum l")')
        write (*, '(" current	   : ", I8)') input%groundstate%lmaxapw
        write (*, '(" in file      : ", I8)') lmaxapw_
        stop
      end if
      if( natmtot_ .ne. natmtot) then
        write(*, '("Error (wannier_readfun): invalid number of atoms")')
        write (*, '(" current	   : ", I8)') natmtot
        write (*, '(" in file      : ", I8)') natmtot_
        stop
      end if
      if( nlotot_ .ne. nlotot) then
        write(*, '("Error (wannier_readfun): invalid number of local-orbitals")')
        write (*, '(" current	   : ", I8)') nlotot
        write (*, '(" in file      : ", I8)') nlotot_
        stop
      end if

      read( un) wnlmo
      read( un) wlmo2l
      read( un) wlmo2m
      read( un) wlmo2o
      read( un) idxlo_

      read( un) wanfmt
      read( un) wanfir

      close( un)
      write(*,*) "WANNIER FUNCTIONS READ"
    
      return
    end subroutine wannier_readfun
    
    subroutine wannier_delfun
    
      integer :: un
      logical :: exist
    
      inquire( file=trim( wf_filename)//"_FUN"//trim( filext), exist=exist)
      if( exist) then
        call getunit( un)
        open( un, file=trim( wf_filename)//"_FUN"//trim( filext))
        close( un, status='delete')
        write(*,*) "WANNIER FUNCTIONS DELETED"
      end if
    
      return
    end subroutine wannier_delfun

    subroutine wannier_reademat( success)
      logical, intent( out) :: success

      ! local variables
      integer :: i, j, is, n, idxn, ik, ix, iy, iz, un
      integer :: fst_, lst_, nst_, ntot_, nkpt_
      real(8) :: vln(3), vkl_( 3), vkl_tmp( 3, wf_kset%nkpt)
      complex(8) :: ztmp

      call getunit( un)

      success = .true.
      inquire( file=trim( wf_filename)//"_EMAT"//trim( filext), exist=success)
      if( .not. success) then
        write(*,*) 'ERROR (wannier_reademat): File '//trim( wf_filename)//"_EMAT"//trim( filext)//' does not exist.'
        return
      end if
      open( un, file=trim( wf_filename)//"_EMAT"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
      read( un) fst_, lst_, nst_, ntot_, nkpt_
      if( (fst_ .gt. wf_fst) .or. (lst_ .lt. wf_lst)) then
        write( *, '(" ERROR (wannier_reademat): bands in input (",I4,":",I4,") out of file band range (",I4,":",I4,").")'), wf_fst, wf_lst, fst_, lst_
        success = .false.
        return
      end if
      if( ntot_ .ne. wf_n_ntot) then
        write( *, '(" ERROR (wannier_reademat): different number of BZ-neighbors in input (",I4,") and file (",I4,").")'), wf_n_ntot, ntot_
        success = .false.
        return
      end if
      if( nkpt_ .ne. wf_kset%nkpt) then
        write( *, '(" ERROR (wannier_reademat): different number of k-points in input (",I4,") and file (",I4,").")'), wf_kset%nkpt, nkpt_
        success = .false.
        return
      end if
      if( allocated( wf_pwmat)) deallocate( wf_pwmat)
      allocate( wf_pwmat( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt, wf_n_ntot))
      do i = 1, wf_n_ntot
        read( un) vln
        ! find index of neighbor
        idxn = 0
        do j = 1, wf_n_nshells
          is = wf_n_usedshells( j)
          do n = 1, wf_n_n( is)
            if( norm2( wf_n_vl( :, n, is) - vln) .lt. input%structure%epslat) then
              idxn = wf_n_ns2n( n, is)
              exit
            end if
          end do
        end do
        if( idxn .gt. 0) then
          do ik = 1, wf_kset%nkpt
            read( un) vkl_
            vkl_tmp( 1, :) = wf_kset%vkl( 1, :) - vkl_( 1)
            vkl_tmp( 2, :) = wf_kset%vkl( 2, :) - vkl_( 2)
            vkl_tmp( 3, :) = wf_kset%vkl( 3, :) - vkl_( 3)
            iz = minloc( norm2( vkl_tmp( :, :), 1), 1)
            if( norm2( wf_kset%vkl( :, iz) - vkl_(:)) .gt. input%structure%epslat) then
              write( *, '(" ERROR (wannier_reademat): k-point in file not in k-point-set.")')
              write( *, '(3F23.6)') vkl_
              success = .false.
              call terminate
            end if
            do iy = fst_, lst_
              do ix = fst_, lst_
                read( un) ztmp
                if( (ix .ge. wf_fst) .and. (ix .le. wf_lst) .and. (iy .ge. wf_fst) .and. (iy .le. wf_lst)) wf_pwmat( ix, iy, ik, idxn) = ztmp
              end do
            end do
          end do
        else
          write( *, '(" ERROR (wannier_reademat): neighboring vector in file not consistent with input.")')
          write( *, '(3F23.6)') vln
        end if
      end do
      close( un)
      if( success) write(*,*) 'Plane-wave matrix-elements successfully read.'
      return
    end subroutine wannier_reademat
    
    ! writes transformation matrices to file
    subroutine wannier_writeemat
      ! local variables
      integer :: i, ik, is, n, idxn, ix, iy, un
      
      call getunit( un)

      open( un, file=trim( wf_filename)//"_EMAT"//trim( filext), action='WRITE', form='UNFORMATTED')
      write( un) wf_fst, wf_lst, wf_nst, wf_n_ntot, wf_kset%nkpt
      do i = 1, wf_n_nshells
        is = wf_n_usedshells( i)
        do n = 1, wf_n_n( is)
          idxn = wf_n_ns2n( n, is)
          write( un) wf_n_vl( :, n, is)
          do ik = 1, wf_kset%nkpt
            write( un) wf_kset%vkl( :, ik)
            do iy = wf_fst, wf_lst
              do ix = wf_fst, wf_lst
                write( un) wf_pwmat( ix, iy, ik, idxn)
              end do
            end do
          end do
        end do
      end do
      close( un)
      write( *, '(a,a)') ' Plane-wave matrix-elements written to file ', trim( wf_filename)//"_EMAT"//trim( filext)
      return
    end subroutine wannier_writeemat  

    subroutine wannier_delemat
      integer :: un
      logical :: exist

      inquire( file=trim( wf_filename)//"_EMAT"//trim( filext), exist=exist)
      if( exist) then
        call getunit( un)
        open( un, file=trim( wf_filename)//"_EMAT"//trim( filext))
        close( un, status='delete')
      end if

      return
    end subroutine wannier_delemat

    !BOP
    ! !ROUTINE: wannier_geometry
    ! !INTERFACE:
    !
    subroutine wannier_geometry
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      ! !DESCRIPTION:
      !   Sets geometry environment for Wannier functions. Finds near neighbors
      !   and corresponding geometric weights as well as the number of shells
      !   needed for gradient calculation.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      integer :: ix, iy, iz, i, j, nshell, nkmax, nposvec, iknr, d, f
      integer :: vi(3)
      real(8) :: vl(3), vc(3), dist
      real(8) :: nvlt( 3, wf_kset%nkpt, wf_kset%nkpt), nvct( 3, wf_kset%nkpt, wf_kset%nkpt)
      real(8) :: coeff(6,6), coeffcpy(6,6), right(6), sval(6), lsvec(6,6), rsvec(6,6), coeffinv(6,6)
      real(8) :: mat1(3,3), mat2(3,3), mat3(3,3)
      logical :: stopshell

      integer :: lwork, info
      real(8), allocatable :: work(:), dt(:)
      integer, allocatable :: iwork(:), mt(:), equik(:), posvec(:,:)
      
      ! find possible distances (shells)
      nposvec = (wf_kset%ngridk(3)-1)*(2*wf_kset%ngridk(2)-1)*(2*wf_kset%ngridk(1)-1)+(wf_kset%ngridk(2)-1)*(2*wf_kset%ngridk(1)-1)+wf_kset%ngridk(1)-1
      allocate( dt( nposvec))
      allocate( mt( nposvec))
      allocate( posvec( 3, nposvec))

      i = 0
      posvec(:,:) = 0
      do iz = -wf_kset%ngridk(3)+1, wf_kset%ngridk(3)-1
        do iy = -wf_kset%ngridk(2)+1, wf_kset%ngridk(2)-1
          do ix = -wf_kset%ngridk(1)+1, wf_kset%ngridk(2)-1
            if( (abs(ix)+abs(iy)+abs(iz)) .ne. 0) then
              stopshell = .true.
              do j = 1, i
                if( (posvec(1,j) .eq. -ix) .and. (posvec(2,j) .eq. -iy) .and. (posvec(3,j) .eq. -iz)) then
                  stopshell = .false.
                  exit
                end if
              end do
              if( stopshell) then
                i = i + 1
                posvec(:,i) = (/ix, iy, iz/)
              end if
            end if
          end do
        end do
      end do

      dt = 0.d0
      mt = 0
      nvlt = 0.d0
      nvct = 0.d0
      i = 0
      do ix = 1, nposvec
        vl = dble( posvec(:,ix))/wf_kset%ngridk
        call r3mv( bvec, vl, vc)
        dist = norm2( vc)
        if( minval( abs( dt(:) - dist)) .gt. input%structure%epslat) then
          i = i + 1
          dt( i) = dist
        end if
        j = minloc( abs( dt(:) - dist), 1)
        if( abs( dt( j) - dist) .lt. input%structure%epslat) then
          mt( j) = mt( j) + 1
        end if
      end do
      iz = i
      nshell = min( iz, 6)

      ! allocating geometry arrays
      if( .not. allocated( wf_n_dist))          allocate( wf_n_dist( nshell))
      if( .not. allocated( wf_n_n))             allocate( wf_n_n( nshell))

      ! sort shells by distance
      wf_n_dist = 0.d0
      wf_n_n = 0
      dist = minval( dt)
      do i = 1, nshell
        j = minloc( abs( dt( 1:iz) - dist), 1)
        wf_n_dist( i) = dt( j)
        wf_n_n( i) = mt( j)
        dt( j) = 1.d13 
        dist = minval( dt)
      end do
      nkmax = maxval( wf_n_n)
      
      ! allocating geometry arrays
      if( .not. allocated( wf_n_ik))            allocate( wf_n_ik( nkmax, nshell, wf_kset%nkptnr))
      if( .not. allocated( wf_n_ik2))           allocate( wf_n_ik2( nkmax, nshell, wf_kset%nkptnr))
      if( .not. allocated( wf_n_vl))            allocate( wf_n_vl( 3, nkmax, nshell))
      if( .not. allocated( wf_n_vc))            allocate( wf_n_vc( 3, nkmax, nshell))
      if( .not. allocated( wf_n_wgt))           allocate( wf_n_wgt( nshell))
      if( .not. allocated( wf_n_usedshells))    allocate( wf_n_usedshells( nshell))
      if( .not. allocated( wf_n_ns2n))          allocate( wf_n_ns2n( nkmax, nshell))
      if( .not. allocated( wf_n_n2neg))         allocate( wf_n_n2neg( nkmax, nshell))

      ! find all possible neighbors
      wf_n_vl = 0.d0
      wf_n_vc = 0.d0
      mt = 0
      do ix = 1, nposvec
        vl = dble( posvec(:,ix))/wf_kset%ngridk
        call r3mv( bvec, vl, vc)
        dist = norm2( vc)
        j = minloc( abs( wf_n_dist(:) - dist), 1)
        if( abs( wf_n_dist( j) - dist) .lt. input%structure%epslat) then
          mt( j) = mt( j) + 1
          wf_n_vl( :, mt( j), j) = vl
          wf_n_vc( :, mt( j), j) = vc
        end if
      end do

      !! find index of negative neighbor
      !do i = 1, nshell
      !  do ix = 1, wf_n_n( i)
      !    dist = 1.d13
      !    iz = 0
      !    do iy = 1, wf_n_n( i)
      !      if( norm2( wf_n_vl( :, ix, i) + wf_n_vl( :, iy, i)) .lt. dist) then
      !        iz = iy
      !        dist = norm2( wf_n_vl( :, ix, i) + wf_n_vl( :, iy, i))
      !      end if
      !    end do
      !    if( (iz .gt. 0) .and. (dist .lt. input%structure%epslat)) then
      !      wf_n_n2neg( ix, i) = iz
      !    else
      !      write( *, '(" ERROR (wannier_geomentry): negative neighbor not in shell.")')
      !      call terminate
      !    end if
      !  end do
      !end do
      !stop
      
      ! find k-point indices of neighbors
      allocate( equik( nsymcrys))
      do iknr = 1, wf_kset%nkpt
        do j = 1, nshell
          do i = 1, wf_n_n( j)
            vl = wf_kset%vkl( :, iknr) + wf_n_vl( :, i, j)
            call r3frac( input%structure%epslat, vl, vi) 
            call myfindkpt( vl, wf_kset, ix, iy)
            if( ix .gt. 0) then
              wf_n_ik( i, j, iknr) = iy
            else
              write( *, '(" ERROR (wannier_geometry): wrong neighboring vector.")')
              write( *, '(3F13.6)') wf_n_vl( :, i, j)
              call terminate
            end if

            vl = wf_kset%vkl( :, iknr) - wf_n_vl( :, i, j)
            call r3frac( input%structure%epslat, vl, vi) 
            call myfindkpt( vl, wf_kset, ix, iy)
            if( ix .gt. 0) then
              wf_n_ik2( i, j, iknr) = iy
            else
              write( *, '(" ERROR (wannier_geometry): wrong neighboring vector.")')
              write( *, '(3F13.6)') wf_n_vl( :, i, j)
              call terminate
            end if
          end do
        end do
      end do
      
      ! find number of shells needed for gradient calculation and geometric weights
      j = 1
      stopshell = .false.
      coeff = 0.d0
      wf_n_wgt = 0.d0
      allocate( work( 1), iwork( 8*j))

      if( minval( wf_kset%ngridk) .eq. 1) then
        if( sum( wf_kset%ngridk) .eq. maxval( wf_kset%ngridk) + 2) then
          d = 1
        else
          d = 3
        end if
      else
        d = 6
      end if
      
      call r3minv( bvec, mat1)
      call r3minv( transpose( bvec), mat2)
      call r3mm( mat1, mat2, mat3)

      if( d .gt. 1) then
        if( d .eq. 3) then
          ix = minloc( wf_kset%ngridk, 1)
          iy = mod( ix+1, 3)
          iz = mod( ix+2, 3)
          if( iy .eq. 0) iy = 3
          if( iz .eq. 0) iz = 3
          right(:) = (/mat3( iy, iy), mat3( iz, iz), mat3( iy, iz), 0.d0, 0.d0, 0.d0/)
        else
          right(:) = (/mat3(1,1), mat3(2,2), mat3(3,3), mat3(1,2), mat3(1,3), mat3(2,3)/)
        end if
        do while( (j .le. d) .and. (.not. stopshell))
          if( d .eq. 3) then
            do i = 1, wf_n_n( j)
              coeff( 1, j) = coeff( 1, j) + wf_n_vl( iy, i, j)**2
              coeff( 2, j) = coeff( 2, j) + wf_n_vl( iz, i, j)**2
              coeff( 3, j) = coeff( 3, j) + wf_n_vl( iy, i, j)*wf_n_vl( iz, i, j)
            end do
          else
            do i = 1, wf_n_n( j)
              coeff( 1, j) = coeff( 1, j) + wf_n_vl( 1, i, j)**2
              coeff( 2, j) = coeff( 2, j) + wf_n_vl( 2, i, j)**2
              coeff( 3, j) = coeff( 3, j) + wf_n_vl( 3, i, j)**2
              coeff( 4, j) = coeff( 4, j) + wf_n_vl( 1, i, j)*wf_n_vl( 2, i, j)
              coeff( 5, j) = coeff( 5, j) + wf_n_vl( 1, i, j)*wf_n_vl( 3, i, j)
              coeff( 6, j) = coeff( 6, j) + wf_n_vl( 2, i, j)*wf_n_vl( 3, i, j)
            end do
          end if
          coeffcpy = coeff
          ! find pseudo-inverse of coeff
          call dgesdd( 'A', d, j, coeffcpy( 1:d, 1:j), d, sval( 1:j), lsvec( 1:d, 1:d), d, rsvec( 1:j, 1:j), j, work, -1, iwork, info)
          lwork = work(1)
          if( allocated( work)) deallocate( work)
          allocate( work( lwork))
          sval = 0.d0
          rsvec = 0.d0
          call dgesdd( 'A', d, j, coeffcpy( 1:d, 1:j), d, sval( 1:j), lsvec( 1:d, 1:d), d, rsvec( 1:j, 1:j), j, work, lwork, iwork, info)
          do i = 1, j
            if( sval( i) .gt. input%structure%epslat) then
              rsvec( i, 1:d) = rsvec( i, 1:j)/sval( i)
            else
              rsvec( i, 1:d) = 0.d0
            end if
          end do
          do i = j+1, d
            rsvec( i, :) = 0.d0
          end do
          coeffinv( 1:j, 1:d) = matmul( transpose( rsvec( 1:d, 1:j)), transpose( lsvec( 1:d, 1:d)))
          sval( 1:d) = matmul( matmul( coeff( 1:d, 1:j), coeffinv( 1:j, 1:d)), right( 1:d))
          if( (sum( abs( sval( 1:d) - right( 1:d))) .lt. input%structure%epslat) .or. (j .eq. d)) then
            stopshell = .true.
            wf_n_wgt( 1:j) = matmul( coeffinv( 1:j, 1:d), right( 1:d))
            wf_n_nshells = j
          end if
          j = j+1
        end do
      else
        wf_n_wgt( 1) = 1.0d0/norm2( wf_n_vc)**2
        wf_n_nshells = 1
      end if
      
      j = 0
      nkmax = 0
      wf_n_usedshells = 0
      wf_n_ntot = 0
      wf_n_ns2n = 0
      do i = 1, wf_n_nshells
        if( abs( wf_n_wgt( i)) .gt. input%structure%epslat) then
          j = j + 1
          wf_n_usedshells( j) = i
          do d = 1, wf_n_n( i)
            wf_n_ntot = wf_n_ntot + 1
            wf_n_ns2n( d, i) = wf_n_ntot
          end do
        end if
      end do
      wf_n_nshells = j
      wf_n_wgt = 0.5d0*wf_n_wgt

      deallocate( dt, mt, work, iwork, posvec)

      return
      !EOC
    end subroutine wannier_geometry
    !EOP
    
    subroutine wannier_readinput
      integer :: iproj, is, ia, ilo, l, m, ik, fst, lst, ist
      real(8) :: e
      real(8), allocatable :: evalfv(:,:)

      if( input%properties%wannier%method .eq. "disentangle") then
        wf_win_i(1) = minval( input%properties%wannier%innerwindow)
        wf_win_i(2) = maxval( input%properties%wannier%innerwindow)
        wf_win_o(1) = minval( input%properties%wannier%outerwindow)
        wf_win_o(2) = maxval( input%properties%wannier%outerwindow)

        if( norm2( wf_win_i) .gt. input%structure%epslat) then
          if( (wf_win_i(1) .lt. wf_win_o(1)) .or. (wf_win_i(2) .gt. wf_win_o(2))) then
            write( *, '(" ERROR (wannier_init): The inner window must be fully contained within the outer window.")')
            call terminate
          end if
        end if

        call wannier_geteval( evalfv, fst, lst)
        allocate( wf_win_ii( fst:lst, wf_kset%nkpt), wf_win_io( fst:lst, wf_kset%nkpt))
        allocate( wf_win_ni( wf_kset%nkpt), wf_win_no( wf_kset%nkpt))
        wf_win_ii = 0
        wf_win_io = 0
        wf_win_ni = 0
        wf_win_no = 0

        wf_fst = 100000000
        wf_lst = 0
        do ik = 1, wf_kset%nkpt
          do ist = fst, lst
            e = evalfv( ist, ik)
            if( (e .ge. wf_win_i(1)) .and. (e .le. wf_win_i(2))) then
              wf_win_ni( ik) = wf_win_ni( ik) + 1
              wf_win_ii( wf_win_ni( ik), ik) = ist
            else if( (e .ge. wf_win_o(1)) .and. (e .le. wf_win_o(2))) then
              wf_win_no( ik) = wf_win_no( ik) + 1
              wf_win_io( wf_win_no( ik), ik) = ist
            end if
          end do
          wf_fst = min( wf_fst, minval( wf_win_ii( 1:wf_win_ni( ik), ik)))
          wf_fst = min( wf_fst, minval( wf_win_io( 1:wf_win_no( ik), ik)))
          wf_lst = max( wf_lst, maxval( wf_win_ii( 1:wf_win_ni( ik), ik)))
          wf_lst = max( wf_lst, maxval( wf_win_io( 1:wf_win_no( ik), ik)))
        end do
        wf_nwf = input%properties%wannier%nwf
        if( wf_nwf .lt. 1) then
          wf_nwf = nint( 0.5d0*(maxval( wf_win_ni) + minval( wf_win_no)))
          write( *, '(" WARNING (wannier_init): Number of Wannier functions (nwf) not set. I chose nwf = ",I3)') wf_nwf
        end if

        do ik = 1, wf_kset%nkpt
          if( wf_win_no( ik) + wf_win_ni( ik) .lt. wf_nwf) then
            write( *, '(" ERROR (wannier_init): Outer window contains less than nwf (",I3,") bands for k-point ",3F13.6,".")') wf_nwf, wf_kset%vkl( :, ik)
            call terminate
          end if
          if( wf_win_ni( ik) .gt. wf_nwf) then
            write( *, '(" ERROR (wannier_init): Inner window contains more than nwf (",I3,") bands for k-point ",3F13.6,".")') wf_nwf, wf_kset%vkl( :, ik)
            call terminate
          end if
        end do

        wf_disentangle = .true.
      else
        wf_fst = min( input%properties%wannier%fst, input%properties%wannier%lst)       
        wf_lst = max( input%properties%wannier%fst, input%properties%wannier%lst)       
        wf_nwf = wf_lst - wf_fst + 1
      end if  

      if( wf_fst .lt. 1) then
        write( *, '(" ERROR (wannier_init): The lowest band (fst) is smaller than 1")')
        call terminate
      end if

      if( wf_lst .gt. nstfv) then
        write( *, '(" ERROR (wannier_init): The highest band (lst = ",I4,") is greater than total number of states (nstfv = ",I4,").")') wf_lst, nstfv
        call terminate
      end if

      wf_nst = wf_lst - wf_fst + 1

      if( input%properties%wannier%input .eq. "gw") then
        if( (wf_fst .lt. input%gw%ibgw) .or. (wf_fst .gt. input%gw%nbgw)) then
          write( *, '(" ERROR (wannier_init): lower band-index (",I4,") out of range (",I4,":",I4,").")'), &
              wf_fst, input%gw%ibgw, input%gw%nbgw-1
          call terminate
        end if
      else
        if( (wf_fst .lt. 1) .or. (wf_fst .gt. nstfv)) then
          write( *, '(" ERROR (wannier_init): lower band-index (",I4,") out of range (",I4,":",I4,").")'), &
              wf_fst, 1, nstfv-1
          call terminate
        end if
      end if

      if( input%properties%wannier%input .eq. "gw") then
        if( wf_lst .gt. input%gw%nbgw) then
          write( *, '(" ERROR (wannier_init): upper band-index (",I4,") out of range (",I4,":",I4,").")'), &
              wf_lst, wf_fst+1, input%gw%nbgw
          call terminate
        end if
      else
        if( wf_lst .gt. nstfv) then
          write( *, '(" ERROR (wannier_init): upper band-index (",I4,") out of range (",I4,":",I4,").")'), &
              wf_lst, wf_fst+1, nstfv
          call terminate
        end if
      end if

      if( allocated( wf_projst)) deallocate( wf_projst)
      allocate( wf_projst( nlotot, 5))
      wf_projst(:,:) = 0
      wf_nprojtot = 0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          do ilo = 1, nlorb( is)
            l = lorbl( ilo, is)
            do m = -l, l
              wf_nprojtot = wf_nprojtot + 1
              wf_projst( wf_nprojtot, 1) = is
              wf_projst( wf_nprojtot, 2) = ia
              wf_projst( wf_nprojtot, 3) = ilo
              wf_projst( wf_nprojtot, 4) = l
              wf_projst( wf_nprojtot, 5) = m 
            end do
          end do
        end do
      end do
      if( wf_nprojtot .eq. 0) then
        write(*,*) 'ERROR (wannier_init): No local-orbitals found for projection.'
        stop
      end if
      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( wf_nprojtot))
      wf_projused = 0
      wf_nprojused = size( input%properties%wannier%projectorarray, 1)

      if( (input%properties%wannier%method .eq. "pro") .or. (input%properties%wannier%method .eq. "promax")) then
        if( wf_nprojused .ne. wf_nst) then
          write( *, '(" ERROR (wannier_init): The number of projectors must be equal to the size of the band-range given.")')
          call terminate
        end if
        do iproj = 1, wf_nst
          if( (input%properties%wannier%projectorarray( iproj)%projector%nr .lt. 1) .or. (input%properties%wannier%projectorarray( iproj)%projector%nr .gt. wf_nprojtot)) then
            write( *, '(" ERROR (wannier_init): ",I4," is not a valid index for projection local-orbitals.")'), input%properties%wannier%projectorarray( iproj)%projector%nr
            write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
            call wannier_showproj
            call terminate
          end if
          wf_projused( input%properties%wannier%projectorarray( iproj)%projector%nr) = 1
        end do
      else if( (input%properties%wannier%method .eq. "opf") .or. (input%properties%wannier%method .eq. "opfmax") .or. (input%properties%wannier%method .eq. "disentangle")) then
        if( wf_nprojused .eq. 0) then
          wf_nprojused = wf_nprojtot
          wf_projused = 1
        else
          if( wf_nprojused .lt. wf_nst) then
            write( *, '(" ERROR (wannier_init): The number of projectors must be greater or equal to the size of the band-range given.")')
            call terminate
          end if
          do iproj = 1, wf_nprojused
            if( (input%properties%wannier%projectorarray( iproj)%projector%nr .lt. 1) .or. (input%properties%wannier%projectorarray( iproj)%projector%nr .gt. wf_nprojtot)) then
              write( *, '(" ERROR (wannier_init): ",I4," is not a valid index for projection local-orbitals.")'), input%properties%wannier%projectorarray( iproj)%projector%nr
              write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
              call wannier_showproj
              call terminate
            end if
            wf_projused( input%properties%wannier%projectorarray( iproj)%projector%nr) = 1
          end do
        end if
      else
      end if

      return
    end subroutine wannier_readinput

    subroutine wannier_setkpts
      !write(*,*) "k set"
      select case (input%properties%wannier%input)
        case( "gs")
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
          wf_kset%vkl = vklnr
          wf_kset%vkc = vkcnr
          call generate_k_vectors( wf_kset_red, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .true.)
        case( "gw")
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
          vkl = wf_kset%vkl
          vkc = wf_kset%vkc
          nkpt = wf_kset%nkpt
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
          vklnr = wf_kset%vkl
          vkcnr = wf_kset%vkc
          nkptnr = wf_kset%nkpt
          call generate_k_vectors( wf_kset_red, bvec, input%gw%ngridq, input%gw%vqloff, .true.)
        case( "qsgw")
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
          vkl = wf_kset%vkl
          vkc = wf_kset%vkc
          nkpt = wf_kset%nkpt
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
          vklnr = wf_kset%vkl
          vkcnr = wf_kset%vkc
          nkptnr = wf_kset%nkpt
          call generate_k_vectors( wf_kset_red, bvec, input%gw%ngridq, input%gw%vqloff, .true.)
        case( "hybrid")
          call readkpts
          !call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)
          !vkl = wf_kset%vkl
          !vkc = wf_kset%vkc
          !nkpt = wf_kset%nkpt
          !! recalculate G+k-vectors
          !do ik = 1, nkpt
          !  do is = 1, nspnfv
          !    vl (:) = vkl (:, ik)
          !    vc (:) = vkc (:, ik)
          !    call gengpvec( vl, vc, ngk( is, ik), igkig( :, is, ik), vgkl( :, :, is, ik), vgkc( :, :, is, ik), gkc( :, is, ik), tpgkc( :, :, is, ik))
          !    call gensfacgp( ngk(is, ik), vgkc( :, :, is, ik), ngkmax, sfacgk( :, :, is, ik))
          !  end do
          !end do
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
          call generate_k_vectors( wf_kset_red, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .true., .false.)
        case default
          write(*, '(" ERROR (wannier_init): ",a," is not a valid input.")') input%properties%wannier%input
          call terminate
      end select
      !write(*,*) "G set"
      !write(*,*) intgv
      call generate_G_vectors( wf_Gset, bvec, intgv, input%groundstate%gmaxvr)
      !write(*,*) "G+k set"
      !write(*,*) gkmax
      call generate_Gk_vectors( wf_Gkset, wf_kset, wf_Gset, gkmax)
      return
    end subroutine wannier_setkpts

    subroutine wannier_genradfun
      character(256) :: fname

      fname = filext

      if( input%properties%wannier%input .eq. 'hybrid') then
        filext = '_PBE.OUT'
      else
        filext = '.OUT'
      end if
      
      call readstate
      call linengy
      call genapwfr     ! APW radial functions
      call genlofr      ! LO radial functions
      call olprad

      filext = fname

      return
    end subroutine wannier_genradfun

    subroutine wannier_projection
      integer :: ik, iproj, is, ias, io, ilo, l, m, lm, ig
      real(8) :: t0, t1

      real(8), allocatable :: rolpi(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), auxmat(:,:)

      call timesec( t0)

      write( wf_info, '(" calculate projection overlap matrices...")')
      allocate( rolpi( wf_nprojtot, apwordmax+nlomax))
      rolpi(:,:) = 0.d0
      do iproj = 1, wf_nprojtot
        is = wf_projst( iproj, 1)
        ias = idxas( wf_projst( iproj, 2), is)
        lm = idxlm( wf_projst( iproj, 4), wf_projst( iproj, 5))
        do io = 1, apword( wf_projst( iproj, 4), is)
          rolpi( iproj, io) = oalo( io, wf_projst( iproj, 3), ias)
        end do
        do ilo = 1, nlorb( is)
          l = lorbl( ilo, is)
          if( l .eq. wf_projst( iproj, 4)) then
            rolpi( iproj, apwordmax+ilo) = ololo( ilo, wf_projst( iproj, 3), ias)
          end if
        end do
      end do

      if( allocated( wf_projection)) deallocate( wf_projection)
      allocate( wf_projection( nstfv, wf_nprojtot, wf_kset%nkpt))

      allocate( evecfv( nmatmax, nstfv, nspinor))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
      allocate( auxmat( nmatmax, wf_nprojtot))

      do ik = 1, wf_kset%nkpt   
        call wannier_getevec( ik, evecfv)
        call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm(:, :, :, :, 1))

        ! projection matrices 
        auxmat(:,:) = zzero
        do iproj = 1, wf_nprojtot
          is = wf_projst( iproj, 1)
          ias = idxas( wf_projst( iproj, 2), is)
          lm = idxlm( wf_projst( iproj, 4), wf_projst( iproj, 5))
          do io = 1, apword( wf_projst( iproj, 4), is)
            auxmat( 1:wf_Gkset%ngk( 1, ik), iproj) = auxmat( 1:wf_Gkset%ngk( 1, ik), iproj) + &
                conjg( apwalm( 1:wf_Gkset%ngk( 1, ik), io, lm, ias, 1))*cmplx( rolpi( iproj, io), 0, 8)
          end do
          do ilo = 1, nlorb( is)
            l = lorbl( ilo, is)
            if( l .eq. wf_projst( iproj, 4)) then
              do m = -l, l
                if( m .eq. wf_projst( iproj, 5)) then
                  ig = wf_Gkset%ngk( 1, ik) + idxlo( lm, ilo, ias)
                  auxmat( ig, iproj) = cmplx( rolpi( iproj, apwordmax+ilo), 0, 8)
                end if
              end do
            end if
          end do
        end do

        call zgemm( 'c', 'n', nstfv, wf_nprojtot, wf_Gkset%ngk( 1, ik)+nlotot, zone, &
             evecfv( :, :, 1), nmatmax, &
             auxmat, nmatmax, zzero, &
             wf_projection( :, :, ik), nstfv)

        !write(*,'(I,3F13.6)') ik, wf_kset%vkl( :, ik)
        !call plotmat( wf_projection( :, :, ik))
        !write(*,*)
      end do 

      deallocate( evecfv, apwalm, auxmat)
      deallocate( rolpi)

      !call writematlab( wf_projection(:,:,12), 'pmfull')
      call timesec( t1)
      write( wf_info, '(" ...projection overlap matrices calculated.")')
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
      write( wf_info, '(5x,"#states: ",T40,7x,I6)') nstfv
      write( wf_info, '(5x,"#projection functions: ",T40,7x,I6)') wf_nprojtot
      write( wf_info, *)
      call flushifc( wf_info)

      return
    end subroutine wannier_projection

    subroutine wannier_removeldprojf( eps)
      real(8), intent( in), optional :: eps

      integer :: i, is, ias, ilo, ilo2, l, lm, j, is2, ias2, l2, lm2, iproj, iproj2, lapack_info, lapack_lwork
      real(8) :: tmp, eps_

      real(8), allocatable :: projolp(:,:)
      real(8), allocatable :: lapack_tau(:), lapack_work(:)

      eps_ = 1.d-3
      if( present( eps)) eps_ = eps

      allocate( projolp( wf_nprojused, wf_nprojused))
      projolp(:,:) = 0.d0
      i = 0
      do iproj = 1, wf_nprojtot
        if( wf_projused( iproj) .eq. 1) then
          i = i + 1
          is = wf_projst( iproj, 1)
          ias = idxas( wf_projst( iproj, 2), is)
          ilo = wf_projst( iproj, 3)
          l = wf_projst( iproj, 4)
          lm = idxlm( l, wf_projst( iproj, 5))
          j = 0
          do iproj2 = 1, wf_nprojtot
            if( wf_projused( iproj2) .eq. 1) then
              j = j + 1
              is2 = wf_projst( iproj2, 1)
              ias2 = idxas( wf_projst( iproj2, 2), is2)
              ilo2 = wf_projst( iproj2, 3)
              l2 = wf_projst( iproj2, 4)
              lm2 = idxlm( l2, wf_projst( iproj2, 5))
              if( (lm .eq. lm2) .and. (ias .eq. ias2)) then
                projolp( i, j) = ololo( ilo, ilo2, ias)
              end if
            end if
          end do
        end if
      end do
      allocate( lapack_tau( wf_nprojused))
      allocate( lapack_work( 1))
      call dgeqrf( wf_nprojused, wf_nprojused, projolp, wf_nprojused, lapack_tau, lapack_work, -1, lapack_info)
      lapack_lwork = nint( lapack_work( 1))
      deallocate( lapack_work)
      allocate( lapack_work( lapack_lwork))
      !call writematlab( cmplx( projolp, 0, 8), 'olp')
      call dgeqrf( wf_nprojused, wf_nprojused, projolp, wf_nprojused, lapack_tau, lapack_work, lapack_lwork, lapack_info)
      !call writematlab( cmplx( projolp, 0, 8), 'qr')
      tmp = 0.d0
      do iproj = 1, wf_nprojused
        tmp = max( tmp, abs( projolp( iproj, iproj)))
      end do
      
      wf_nprojused = 0
      i = 0
      do iproj = 1, wf_nprojtot
        if( wf_projused( iproj) .eq. 1) then
          wf_projused( iproj) = 0
          i = i + 1
          !write(*,'(I4,3x,F13.6)') iproj, abs( projolp( i, i))/tmp
          if( abs( projolp( i, i))/tmp .gt. eps_) then
            wf_nprojused = wf_nprojused + 1
            wf_projused( iproj) = 1
          else
            wf_projused( iproj) = -1
          end if
        end if
      end do
      deallocate( projolp, lapack_tau, lapack_work)
      write(*,*) wf_nprojused, wf_nprojtot
      if( wf_nprojused .lt. wf_nst) then
        write(*,'(" ERROR (wannier_init): Number of linear independent projection functions (",I4,") is smaller than number of bands (",I4,"). Please add more projection functions.")') wf_nprojused, wf_nst
        !call terminate
      end if

      return
    end subroutine wannier_removeldprojf

    subroutine wannier_getevec( ik, evec)
      integer, intent( in) :: ik
      complex(8), intent( out) :: evec( nmatmax, nstfv, nspinor)

      if( (input%properties%wannier%input .eq. "gs") .or. (input%properties%wannier%input .eq. "qsgw")) then
        call getevecfv( wf_kset%vkl(:, ik), wf_Gkset%vgkl( :, :, :, ik), evec)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl(:, ik), wf_GKset%vgkl( :, :, :, ik), evec)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspinor, evec)
      else
        call terminate
      end if

      return
    end subroutine wannier_getevec

    subroutine wannier_geteval( eval, fst, lst)
      real(8), allocatable, intent( out) :: eval(:,:)
      integer, intent( out) :: fst, lst

      integer :: ik, ist, un, recl, nkpqp, fstqp, lstqp, nk, isym, iq, nkequi, isymequi( wf_kset%nkpt), ikequi( wf_kset%nkpt)
      real(8) :: vl(3), efermiqp, efermiks
      character(256) :: fname
      logical :: exist

      real(8), allocatable :: evalqp(:), evalks(:), evalfv(:,:)

      if( input%properties%wannier%input .eq. "gw") then
        call getunit( un)
        write( fname, '("EVALQP.OUT")')
        inquire( file=trim( fname), exist=exist)
        if( .not. exist) then
          write( *, '("Error (wfint_init): File EVALQP.OUT does not exist!")')
          call terminate
        end if
        inquire( iolength=recl) nkpqp, fstqp, lstqp
        open( un, file=trim( fname), action='read', form='unformatted', access='direct', recl=recl)
        read( un, rec=1), nkpqp, fstqp, lstqp
        close( un)
        allocate( evalqp( fstqp:lstqp))
        allocate( evalks( fstqp:lstqp))
        allocate( eval( fstqp:lstqp, wf_kset%nkpt))
        fst = fstqp
        lst = lstqp
        inquire( iolength=recl) nkpqp, fstqp, lstqp, vl, evalqp, evalks, efermiqp, efermiks
        open( un, file=trim( fname), action='read', form='unformatted', access='direct', recl=recl)
        do ik = 1, nkpqp
          read( un, rec=ik) nkpqp, fstqp, lstqp, vl, evalqp, evalks, efermiqp, efermiks
          call findequivkpt( vl, wf_kset, nkequi, isymequi, ikequi)
          do iq = 1, nkequi
            eval( :, ikequi( iq)) = evalqp(:)
          end do
        end do
        close( un)
        deallocate( evalqp, evalks)
      else
        allocate( evalfv( nstfv, nspnfv))
        allocate( eval( nstfv, wf_kset%nkpt))
        fst = 1
        lst = nstfv
        do ik = 1, wf_kset%nkpt
          call getevalfv( wf_kset%vkl( :, ik), evalfv)
          eval( :, ik) = evalfv( :, 1)
        end do
        deallocate( evalfv)
      end if

      return
    end subroutine wannier_geteval

    subroutine wannier_writeinfo_lo
      integer :: i

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call printbox( wf_info, '*', "Local-orbitals for projection")
      write( wf_info, *)
      write( wf_info, '(12x,"#",6x,"species",9x,"atom",12x,"l",12x,"m",9x,"used")')
      write( wf_info, '(80("-"))')
      
      do i = 1, wf_nprojtot
        write( wf_info, '(9x,I4,11x,a2,10x,I3,11x,I2,11x,I2,10x)', advance='no') &
            i, &
            spsymb( wf_projst( i, 1)), &
            wf_projst( i, 2), &
            wf_projst( i, 4), &
            wf_projst( i, 5)
        if( wf_projused( i) .eq. 1) then
          write( wf_info, '("[X]")', advance='no')
        else if( wf_projused( i) .eq. -1) then
          write( wf_info, '("[-]")', advance='no')
        else
          write( wf_info, '("[ ]")', advance='no')
        end if
        write( wf_info, *)
      end do

      write( wf_info, '(80("-"))')
      write( wf_info, '(36x,"local-orbitals used in total:",4x,I4,"/",I4)') wf_nprojused, wf_nprojtot
      write( wf_info, *)
      call flushifc( wf_info)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
    end subroutine wannier_writeinfo_lo
    
    subroutine wannier_writeinfo_geometry
      integer :: i, j
      real(8) :: v(3,1), m(3,3)

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call printbox( wf_info, '*', "Brillouin zone neighbors for k-gradient")
      write( wf_info, *)
      write( wf_info, '(11x,"shell",6x,"#neighbors",8x,"distance",10x,"weight",10x,"used")')
      write( wf_info, '(80("-"))')

      do i = 1, min( size( wf_n_dist), 6)
        write( wf_info, '(15x,I1,13x,I3,F16.10,F16.10,11x)', advance='no') &
            i, &
            wf_n_n( i), &
            wf_n_dist( i), &
            wf_n_wgt( i)
        if( any( wf_n_usedshells .eq. i, 1)) then
          write( wf_info, '("[X]")', advance='no')
        else
          write( wf_info, '("[ ]")', advance='no')
        end if
        write( wf_info, *)
      end do

      write( wf_info, '(80("-"))')
      write( wf_info, '(40x,"neighbors used in total:",11x,I3)') wf_n_ntot
      write( wf_info, *)
      
      write( wf_info, '(" vectors to neighboring k-points (cartesian)")')
      write( wf_info, *)
      
      m = 0.d0
      do i = 1, wf_n_nshells
        write( wf_info, '(5x,"shell ",I1)') wf_n_usedshells( i)
        do j = 1, wf_n_n( wf_n_usedshells( i))
          write( wf_info, '(12x,2I2,3(F22.10))') j, wf_n_ns2n( j, wf_n_usedshells( i)), wf_n_vc( :, j, wf_n_usedshells( i))
          v(:,1) = wf_n_vc( :, j, wf_n_usedshells( i))
          m = m + 2.0d0*wf_n_wgt( wf_n_usedshells( i))*matmul( v, transpose( v)) 
        end do
      end do

      write( wf_info, *)
      write( wf_info, '(" vectors to neighboring k-points (lattice)")')
      write( wf_info, *)
      
      do i = 1, wf_n_nshells
        write( wf_info, '(5x,"shell ",I1)') wf_n_usedshells( i)
        do j = 1, wf_n_n( wf_n_usedshells( i))
          write( wf_info, '(12x,2I2,3(F22.10))') j, wf_n_ns2n( j, wf_n_usedshells( i)), wf_n_vl( :, j, wf_n_usedshells( i))
        end do
      end do

      write( wf_info, *)
      write( wf_info, '(" consistency check (fullfilled if unity, consider grid-dimension)")')
      write( wf_info, *)
      write( wf_info, '(14x,3(F22.10))') m(1,:)
      write( wf_info, '(14x,3(F22.10))') m(2,:)
      write( wf_info, '(14x,3(F22.10))') m(3,:)
      

      write( wf_info, *)
      call flushifc( wf_info)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
    end subroutine wannier_writeinfo_geometry
    
    subroutine wannier_writeinfo_task( task)
      character(*), intent( in) :: task
     
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call printbox( wf_info, '*', "starting Wannierization...")
      write( wf_info, *)
      write( wf_info, '(" lowest band:",T30,13x,I3)') wf_fst
      write( wf_info, '(" highest band:",T30,13x,I3)') wf_lst
      write( wf_info, '(" #Wannier functions:",T30,13x,I3)') wf_nwf
      write( wf_info, '(" method:",T30,A16)') trim( task)
      write( wf_info, *)
      call flushifc( wf_info)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
    end subroutine wannier_writeinfo_task
    
    subroutine wannier_writeinfo_finish
      real(8) :: t
      
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call timesec( t)
      write( wf_info, '(" total duration (seconds):",T30,F16.1)') t-wf_t0
      call printbox( wf_info, '*', "...Wannierization finished")
      write( wf_info, *)
      call flushifc( wf_info)
      call wannier_writeinfo_results
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
    end subroutine wannier_writeinfo_finish
    
    subroutine wannier_writeinfo_results
      use mod_lattice, only: ainv 
      integer :: i
      real(8) :: vl(3)

      call printbox( wf_info, '*', "Wannier functions")
      write( wf_info, *)
      write( wf_info, '(6x,"band",17x,"localization center (lattice)",5x,"Omega (bohr^2)",3x,"Omega_I (bohr^2)",3x,"Omega_D (bohr^2)",2x,"Omega_OD (bohr^2)")')
      write( wf_info, '(80("-"))')

      do i = 1, wf_nwf
        call r3mv( ainv, wf_centers( :, i), vl)
        write( wf_info, '(6x,I4,7x,3F13.6,4(3x,F16.6))') i, vl, wf_omega( i), wf_omega_i( i), wf_omega_d( i), wf_omega_od( i)
      end do

      write( wf_info, '(80("-"))')
      write( wf_info, '(50x,"total:",4(3x,F16.6))') sum( wf_omega), sum( wf_omega_i), sum( wf_omega_d), sum( wf_omega_od)
      write( wf_info, '(48x,"average:",4(3x,F16.6))') sum( wf_omega)/wf_nwf, sum( wf_omega_i)/wf_nwf, sum( wf_omega_d)/wf_nwf, sum( wf_omega_od)/wf_nwf

      write( wf_info, *)
      call flushifc( wf_info)
      !close( wf_info)
    end subroutine wannier_writeinfo_results

    subroutine wannier_destroy
      if( allocated( wf_projst)) deallocate( wf_projst)
      if( allocated( wf_projused)) deallocate( wf_projused)
      if( allocated( wf_transform)) deallocate( wf_transform)
      if( allocated( wf_projection)) deallocate( wf_projection)
      if( allocated( wf_opf)) deallocate( wf_opf)
      if( allocated( wf_pwmat)) deallocate( wf_pwmat)
      if( allocated( wf_m0)) deallocate( wf_m0)
      if( allocated( wf_m)) deallocate( wf_m)
      if( allocated( wf_subspace)) deallocate( wf_subspace)
      if( allocated( wf_centers)) deallocate( wf_centers)
      if( allocated( wf_omega)) deallocate( wf_omega)
      if( allocated( wf_omega_i)) deallocate( wf_omega_i)
      if( allocated( wf_omega_d)) deallocate( wf_omega_d)
      if( allocated( wf_omega_od)) deallocate( wf_omega_od)
      if( allocated( wf_n_dist)) deallocate( wf_n_dist)
      if( allocated( wf_n_n)) deallocate( wf_n_n)
      if( allocated( wf_n_ns2n)) deallocate( wf_n_ns2n)
      if( allocated( wf_n_n2neg)) deallocate( wf_n_n2neg)
      if( allocated( wf_n_ik)) deallocate( wf_n_ik)
      if( allocated( wf_n_ik2)) deallocate( wf_n_ik2)
      if( allocated( wf_n_vl)) deallocate( wf_n_vl)
      if( allocated( wf_n_vc)) deallocate( wf_n_vc)
      if( allocated( wf_n_wgt)) deallocate( wf_n_wgt)
      if( allocated( wf_n_usedshells)) deallocate( wf_n_usedshells)
      wf_initialized = .false.
    end subroutine wannier_destroy
end module mod_wannier
