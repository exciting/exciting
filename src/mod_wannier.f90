module mod_wannier
  use mod_APW_LO
  use mod_atoms 
  use mod_kpoint
  use mod_constants
  use mod_muffin_tin
  use mod_Gkvector
  use mod_kpointset
  use mod_eigensystem, only: nmatmax, idxlo
  use mod_spin, only: nspinor, nspnfv
  use mod_lattice, only: bvec
  use mod_eigenvalue_occupancy, only: nstfv
  use modgw, only: kset
  use m_getunit
  use m_ematqk
  use m_plotmat
  use modmpi
  use mod_misc

  implicit none

! module variables
  integer :: wf_nprojtot                                ! total number of projection functions
  integer :: wf_fst                                     ! lower band-boundary for generation
  integer :: wf_lst                                     ! upper band-boundary for generation
  integer :: wf_nst                                     ! number of states within boundaries
  integer :: wf_info                                    ! unit to write WANNIERINFO.OUT
  integer :: wf_n_nshells                               ! number of shells used for gradient calculation
  integer :: wf_n_ntot                                  ! total numbers of neighbors in used shells
  logical :: wf_initialized = .false.
  logical :: wf_disentangle = .false.
  real(8) :: wf_t0
  real(8) :: wf_win_i(2), wf_win_o(2)
  type( k_set) :: wf_kset
  
  integer, allocatable :: wf_projst(:,:)                ! l, m, species, atom of local-orbitals for projection
  integer, allocatable :: wf_projused(:)                ! flag wether a local-obital was used for projection
  integer, allocatable :: wf_win_ni(:), wf_win_no(:), wf_win_idxi(:,:), wf_win_idxo(:,:)
  complex(8), allocatable :: wf_transform(:,:,:)        ! unitary transformation matrices
  complex(8), allocatable :: wf_projection(:,:,:)       ! full overlap of wafefunctions and local-orbitals
  complex(8), allocatable :: wf_opf(:,:)                ! expansion coefficients for optimal projection functions
  complex(8), allocatable :: wf_emat(:,:,:,:)           ! plane-wave matrix elements for neighboring k-points
  complex(8), allocatable :: wf_subspace(:,:,:)         ! selected subspace for band-disentanglement
  real(8), allocatable :: wf_centers(:,:)               ! centers of Wannier functions
  real(8), allocatable :: wf_omega(:)                   ! localization functional for each Wannier function
  real(8), allocatable :: wf_n_dist(:)                  ! distance of neighbors in given shell
  integer, allocatable :: wf_n_n(:)                     ! number of neighbors in given shell
  integer, allocatable :: wf_n_ik(:,:,:)                ! k-point index of given neighbor in given shell for given k-point
  integer, allocatable :: wf_n_usedshells(:)            ! index of used shells
  integer, allocatable :: wf_n_ns2n(:,:)                ! map from shell and neighbor in shell to total neighbor index
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
      ! !DESCRIPTION:
      !   Reads local-orbitals from species files which are indicated to be used
      !   as projection functions for the generation of Wannier functions to the
      !   module {\tt mod\_wannier} and constructs necessary geometry (nearest 
      !   neighbors and geometric weights) as well as the overlap of the
      !   wavefunctions and the local-orbitals.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      ! local variables
      integer :: iproj, is, ia, js, ja, ias, nr, l, m, lm, io, ilo, ir, ig
      integer :: ik, iknr, ngknr
      real(8) :: ecenter, t0, t1, vl(3), vc(3)
      integer, allocatable :: igkignr(:), cnt(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:), sval(:)
      complex(8), allocatable :: sfacgknr(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), auxmat(:,:), lsvec(:,:), rsvec(:,:) 
      real(8), allocatable :: rolpi(:,:), uf(:), gf(:), cf(:,:), evalfv(:,:)
      
      write(*,*) "wannier_init"
      call init0
      call init1

      !********************************************************************
      ! open WANNIERINFO.OUT file
      !********************************************************************
      call getunit( wf_info)
      open( wf_info, file='WANNIERINFO'//trim(filext), action='WRITE', form='FORMATTED')
      call timesec( wf_t0)

      !********************************************************************
      ! read projection functions
      !********************************************************************
      allocate( wf_projst( nlotot, 5))
      wf_projst(:,:) = 0
      wf_nprojtot = 0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          do ilo = 1, nlorb( is)
            !if( lorbwfproj( ilo, is)) then
              l = lorbl( ilo, is)
              do m = -l, l
                wf_nprojtot = wf_nprojtot + 1
                wf_projst( wf_nprojtot, 1) = is
                wf_projst( wf_nprojtot, 2) = ia
                wf_projst( wf_nprojtot, 3) = ilo
                wf_projst( wf_nprojtot, 4) = l
                wf_projst( wf_nprojtot, 5) = m 
              end do
            !end if
          end do
        end do
      end do
      if( wf_nprojtot .eq. 0) then
        write(*,*) 'ERROR (wannier_init): No local-orbitals found for projection.'
        stop
      end if
      allocate( wf_projused( wf_nprojtot))
      wf_projused = 0

      !********************************************************************
      ! check input parameters
      !********************************************************************
      if( (input%properties%wannier%fst .lt. 1) .or. (input%properties%wannier%fst .gt. nstfv)) then
        write( *, '(" ERROR (wannier_init): lower band-index (",I4,") out of range (",I4,":",I4,").")'), input%properties%wannier%fst, 1, nstfv
        call terminate
      else
        wf_fst = input%properties%wannier%fst
      end if
      if( (input%properties%wannier%lst .lt. wf_fst) .or. (input%properties%wannier%lst .gt. nstfv)) then
        write( *, '(" ERROR (wannier_init): upper band-index (",I4,") out of range (",I4,":",I4,").")'), input%properties%wannier%lst, wf_fst, nstfv
        call terminate
      else
        wf_lst = input%properties%wannier%lst
      end if
      wf_nst = wf_lst-wf_fst+1

      if( (input%properties%wannier%method .eq. "pro") .or. (input%properties%wannier%method .eq. "promax")) then
        if( size( input%properties%wannier%projectors, 1) .ne. wf_nst) then
          write( *, '(" ERROR (wannier_init): The number of projectors must be equal to the size of the band-range given.")')
          call terminate
        end if
        do iproj = 1, wf_nst
          if( (input%properties%wannier%projectors( iproj) .lt. 1) .or. (input%properties%wannier%projectors( iproj) .gt. wf_nprojtot)) then
            write( *, '(" ERROR (wannier_init): ",I4," is not a valid index for projection local-orbitals.")'), input%properties%wannier%projectors( iproj)
            write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
            call wannier_showproj
            call terminate
          end if
          wf_projused( input%properties%wannier%projectors( iproj)) = 1
        end do
      else if( (input%properties%wannier%method .eq. "opf") .or. (input%properties%wannier%method .eq. "opfmax")) then
        wf_projused = 1
      else
      end if
      call wannier_writeinfo_lo

      !********************************************************************
      ! load used k-grid
      !********************************************************************
      select case (input%properties%wannier%input)
        case( "groundstate")
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
          wf_kset%vkl = vklnr
          wf_kset%vkc = vkcnr
          write(*,*) wf_kset%nkpt, nkpt, nkptnr
        case( "gw")
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
          vkl = wf_kset%vkl
          vkc = wf_kset%vkc
          nkpt = wf_kset%nkpt
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
        case( "hybrid")
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)
          vkl = wf_kset%vkl
          vkc = wf_kset%vkc
          nkpt = wf_kset%nkpt
          ! recalculate G+k-vectors
          do ik = 1, nkpt
            do is = 1, nspnfv
              vl (:) = vkl (:, ik)
              vc (:) = vkc (:, ik)
              call gengpvec( vl, vc, ngk( is, ik), igkig( :, is, ik), vgkl( :, :, is, ik), vgkc( :, :, is, ik), gkc( :, is, ik), tpgkc( :, :, is, ik))
              call gensfacgp( ngk(is, ik), vgkc( :, :, is, ik), ngkmax, sfacgk( :, :, is, ik))
            end do
          end do
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
        case default
          write(*, '(" ERROR (wannier_init): ",a," is not a valid input.")') input%properties%wannier%input
          call terminate
      end select

      !********************************************************************
      ! find geometry
      !********************************************************************
      call wannier_geometry
      call wannier_writeinfo_geometry
      
      call wannier_writeinfo_task( input%properties%wannier%method)
      
      ! analyze energy windows for disentanglement
      !if( norm2( input%properties%wannier%outerwindow) .gt. 0.d0) then
      !  wf_win_o(1) = minval( input%properties%wannier%outerwindow)
      !  wf_win_o(2) = maxval( input%properties%wannier%outerwindow)
      !  if( wf_win_o(2) - wf_win_o(1) .lt. input%structure%epslat) then
      !    write( *, '(" ERROR (wannier_init): The given outer energy-window is too narrow. Width: ",F13.6," Hartree")') wf_win_o(2) - wf_win_o(1)
      !    call terminate
      !  end if
      !  if( norm2( input%properties%wannier%innerwindow) .gt. 0.d0) then
      !    wf_win_i(1) = minval( input%properties%wannier%innerwindow)
      !    wf_win_i(2) = maxval( input%properties%wannier%innerwindow)
      !    if( wf_win_i(2) - wf_win_i(1) .lt. input%structure%epslat) then
      !      write( *, '(" ERROR (wannier_init): The given inner energy-window is too narrow. Width: ",F13.6," Hartree")') wf_win_i(2) - wf_win_i(1)
      !      call terminate
      !    end if
      !    if( (wf_win_i(1) .lt. wf_win_o(1)) .or. (wf_win_i(2) .gt. wf_win_o(2))) then
      !      write( *, '(" ERROR (wannier_init): The inner energy-window is not fully contained inside the outer energy-window.")')
      !      call terminate
      !    end if
      !  else
      !    wf_win_i = wf_win_o(1)-0.1d0
      !  end if
      !  allocate( wf_win_ni( wf_kset%nkpt), wf_win_no( wf_kset%nkpt))
      !  allocate( wf_win_idxi( nstfv, wf_kset%nkpt), wf_win_idxo( nstfv, wf_kset%nkpt))
      !  allocate( evalfv( nstfv, nspinor), cnt( nstfv))
      !  wf_win_ni = 0
      !  wf_win_no = 0
      !  wf_win_idxi = 0
      !  wf_win_idxo = 0
      !  cnt = 0
      !  if( abs( wf_win_i(2) - wf_win_i(1)) .lt. input%structure%epslat) then
      !    ecenter = 0.5d0*sum( wf_win_o)
      !  else
      !    ecenter = 0.5d0*sum( wf_win_i)
      !  end if
      !  write(*,*) ecenter
      !  do iknr = 1, wf_kset%nkpt
      !    call findkpt( wf_kset%vkl( :, iknr), ig, ik)
      !    call getevalfv( vkl( :, ik), evalfv)
      !    do is = 1, nstfv
      !      if( (evalfv( is, 1) .ge. wf_win_o(1)) .and. (evalfv( is, 1) .le. wf_win_o(2))) then
      !        if( (evalfv( is, 1) .ge. wf_win_i(1)) .and. (evalfv( is, 1) .le. wf_win_i(2))) then
      !          wf_win_ni( iknr) = wf_win_ni( iknr) + 1
      !          wf_win_idxi( wf_win_ni( iknr), iknr) = is
      !        else
      !          wf_win_no( iknr) = wf_win_no( iknr) + 1
      !          wf_win_idxo( wf_win_no( iknr), iknr) = is
      !        end if
      !      end if
      !    end do
      !    l = minloc( abs( evalfv( :, 1) - ecenter), 1)
      !    cnt( l) = cnt( l) + 1
      !    if( wf_win_no( iknr)+wf_win_ni( iknr) .lt. input%properties%wannier%nst) then
      !      write( *, '(" ERROR (wannier_init): The outer energy-window must contain at least as many bands as the number of states given.")')
      !      write( *, '(" states: ",T20,I4)') input%properties%wannier%nst
      !      write( *, '(" bands in window: ",T20,I4)') wf_win_no( iknr)+wf_win_ni( iknr)
      !      write( *, '(" k-point: ",T20,I4)') ik
      !      call terminate
      !    end if
      !    if( wf_win_ni( iknr) .gt. input%properties%wannier%nst) then
      !      write( *, '(" ERROR (wannier_init): The inner energy-window must contain maximally as many bands as the number of states given.")')
      !      write( *, '(" states: ",T20,I4)') input%properties%wannier%nst
      !      write( *, '(" bands in window: ",T20,I4)') wf_win_ni( iknr)
      !      write( *, '(" k-point: ",T20,I4)') ik
      !      call terminate
      !    end if
      !    !write(*,'(I3)') iknr
      !    !write(*,'(100I3)') wf_win_idxo( 1:wf_win_no( iknr), iknr)
      !    !write(*,'(100I3)') wf_win_idxi( 1:wf_win_ni( iknr), iknr)
      !    !write(*,'(2I3)') wf_win_ni( iknr), wf_win_no( iknr)
      !    !write(*,'("-------")')
      !  end do
      !  write(*,*) cnt
      !  wf_disentangle = .true.
      !end if
      
      write( wf_info, '(" calculate projection overlap matrices...")')
      call timesec( t0)
      !********************************************************************
      ! radial overlap integrals
      !********************************************************************
      call readstate
      call linengy
      call genapwfr     ! APW radial functions
      call genlofr      ! LO radial functions

      allocate( rolpi( wf_nprojtot, apwordmax+maxlorb))
      allocate( uf( nrmtmax), gf( nrmtmax), cf( 3, nrmtmax))
      rolpi(:,:) = zzero
      do iproj = 1, wf_nprojtot
        uf(:) = 0d0
        is = wf_projst( iproj, 1)
        nr = nrmt( is)
        ias = idxas( wf_projst( iproj, 2), is)
        lm = idxlm( wf_projst( iproj, 4), wf_projst( iproj, 5))
        do io = 1, apword( wf_projst( iproj, 4), is)
          do ir = 1, nr
            uf( ir) = apwfr( ir, 1, io, wf_projst( iproj, 4), ias)*lofr( ir, 1, wf_projst( iproj, 3), ias)*spr( ir, is)**2
          end do
          call fderiv( -1, nr, spr(:, is), uf, gf, cf)
          rolpi( iproj, io) = gf( nr)
        end do
        do ilo = 1, nlorb( is)
          l = lorbl( ilo, is)
          if( l .eq. wf_projst( iproj, 4)) then
            do ir = 1, nr
              uf( ir) = lofr( ir, 1, ilo, ias)*lofr( ir, 1, wf_projst( iproj, 3), ias)*spr( ir, is)**2
            end do
            call fderiv( -1, nr, spr( :, is), uf, gf, cf)
            rolpi( iproj, apwordmax+ilo) = gf( nr)
          end if
        end do
      end do

      !********************************************************************
      ! build full projection overlap matrices
      !********************************************************************
      allocate( evecfv( nmatmax, nstfv, nspinor))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
      allocate( wf_projection( nstfv, wf_nprojtot, wf_kset%nkpt))
      allocate( auxmat( ngkmax+nlotot, wf_nprojtot))
      allocate( sval( min( nstfv, wf_nprojtot)), &
                lsvec( nstfv, nstfv), &
                rsvec( wf_nprojtot, wf_nprojtot))
      allocate( igkignr( ngkmax))
      allocate( vgklnr( 3, ngkmax, nspinor), vgkcnr( 3, ngkmax, nspinor), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
      allocate( sfacgknr( ngkmax, natmtot))

      do iknr = 1, wf_kset%nkpt
        vgklnr = zzero
        ! find G+k-vectors for non-reduced k-point
        call gengpvec( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
        ! find structure factors
        call gensfacgp( ngknr, vgkcnr, ngkmax, sfacgknr)
        ! get basis function coefficients and matching coefficients
        call match( ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm(:, :, :, :, 1))
        if( input%properties%wannier%input .eq. "groundstate") then
          call getevecfv( wf_kset%vkl(:, iknr), vgklnr, evecfv)
        else if( input%properties%wannier%input .eq. "hybrid") then
          call getevecfv( wf_kset%vkl(:, iknr), vgklnr, evecfv)
        else if( input%properties%wannier%input .eq. "gw") then
          call getevecsvgw_new( "GW_EVECSV.OUT", iknr, wf_kset%vkl( :, iknr), nmatmax, nstfv, nspinor, evecfv)
        else
          call terminate
        end if

        ! projection matrices 
        do iproj = 1, wf_nprojtot
          is = wf_projst( iproj, 1)
          ias = idxas( wf_projst( iproj, 2), is)
          lm = idxlm( wf_projst( iproj, 4), wf_projst( iproj, 5))
          do ig = 1, ngknr
            auxmat( ig, iproj) = zzero
            do io = 1, apword( wf_projst( iproj, 4), is)
              auxmat( ig, iproj) = auxmat( ig, iproj) + conjg( apwalm( ig, io, lm, ias, 1))*rolpi( iproj, io)
            end do
          end do
          do js = 1, nspecies
            do ja = 1, natoms( js)
              do ilo = 1, nlorb( js)
                l = lorbl( ilo, js)
                do m = -l, l
                  ig = ngknr + idxlo( idxlm( l, m), ilo, idxas( ja, js))
                  auxmat( ig, iproj) = zzero
                  if( (idxas( ja, js) .eq. ias) .and. (idxlm( l, m) .eq. lm)) then
                    auxmat( ig, iproj) = cmplx( rolpi( iproj, apwordmax+ilo), 0, 8)
                  end if
                end do
              end do
            end do
          end do
        end do

        call ZGEMM( 'C', 'N', nstfv, wf_nprojtot, ngknr+nlotot, zone, &
             evecfv( 1:(ngknr+nlotot), :, 1), ngknr+nlotot, &
             auxmat( 1:(ngknr+nlotot), :), ngknr+nlotot, zzero, &
             wf_projection( :, :, iknr), nstfv)
      end do 
      call timesec( t1)
      write( wf_info, '(" ...projection overlap matrices calculated.")')
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
      write( wf_info, '(5x,"#states: ",T40,7x,I6)') nstfv
      write( wf_info, '(5x,"#projection functions: ",T40,7x,I6)') wf_nprojtot
      write( wf_info, *)
      call flushifc( wf_info)

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
        case default
      end select
      
      wf_initialized = .true.
      !call wannier_showproj
      deallocate( auxmat, rolpi, apwalm, evecfv, uf, gf, cf, sfacgknr, vgklnr, vgkcnr, gkcnr, tpgkcnr, igkignr, sval, lsvec, rsvec)
      if( wf_disentangle) call wannier_subspace( maxloc( cnt, 1))
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
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      ! local variables
      integer :: i, is, n, k1, k2, iknr, ngknr, maxn, idxn, n2
      real(8) :: t0, t1
      logical :: succes

      ! allocatable arrays
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:), dist(:)
      !complex(8), allocatable :: evecfv_tmp(:,:,:,:)
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)
      
      write(*,*) "wannier_emat..."
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(" calculate plane-wave matrix-elements...")')
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      call timesec( t0)

      ! check for existing file
      succes = .true.
      inquire( file="WANNIER_EMAT"//trim( filext), exist=succes)
      if( succes) then
        call wannier_reademat( "WANNIER", succes)
        if( succes) then
          call timesec( t1)
          write( wf_info, '(" ...plane-wave matrix-elements read from file.")')
          write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
          write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
          write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
        end if
      end if

      if( .not. succes) then          
        maxn = 0
        do i = 1, wf_n_nshells
          maxn = max( maxn, wf_n_n( wf_n_usedshells( i)))
        end do
        write(*,*) maxn

        if( allocated( wf_emat)) deallocate( wf_emat)
        !allocate( wf_emat( nstfv, nstfv, maxval( wf_n_n( 1:wf_n_nshells)), wf_n_nshells, wf_kset%nkpt))
        allocate( wf_emat( wf_fst:wf_lst, wf_fst:wf_lst, wf_n_ntot, wf_kset%nkpt))
        wf_emat = zzero
        allocate( dist( wf_n_ntot))
        
        ! read eigenvectors
        allocate( igkignr( ngkmax))
        allocate( vgklnr( 3, ngkmax, nspinor), vgkcnr( 3, ngkmax, nspinor), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
        !allocate( evecfv_tmp( nmatmax, wf_fst:wf_lst, nspinor, wf_kset%nkpt))
        allocate( evecfv1( nmatmax, nstfv, nspinor))
        allocate( evecfv2( nmatmax, nstfv, nspinor))
        !do iknr = 1, wf_kset%nkpt 
        !  if( input%properties%wannier%input .eq. "groundstate") then
        !    ! find G+k-vectors and eigenvectors for non-reduced k-point k
        !    call gengpvec( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
        !    call getevecfv( wf_kset%vkl( :, iknr), vgklnr, evecfv1)
        !    evecfv_tmp( :, :, :, iknr) = evecfv1( :, wf_fst:wf_nst, :)
        !  else if( input%properties%wannier%input .eq. "gw") then
        !    call getevecsvgw_new( "GW_EVECSV.OUT", iknr, wf_kset%vkl( :, iknr), nmatmax, nstfv, nspinor, evecfv1)
        !    evecfv_tmp( :, :, :, iknr) = evecfv1( :, wf_fst:wf_nst, :)
        !  else
        !    call terminate
        !  end if
        !end do
        !deallocate( igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
        !write( *, '("size of evecs: ",I)') sizeof( evecfv_tmp)

        do i = 1, wf_n_nshells
          is = wf_n_usedshells( i) 
          write(*,*) is 
          do n = 1, wf_n_n( is)/2
            idxn = wf_n_ns2n( n, is)
            call emat_init( wf_n_vl( :, n, is), (/0, 0, 0/), input%groundstate%lmaxapw, 8)
            k1 = 1
            k2 = wf_kset%nkpt
!#ifdef MPI
!          k1 = firstofset( rank, wf_kset%nkpt)
!          k2 = lastofset( rank, wf_kset%nkpt)
!#endif
#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr)
!!$OMP DO  
#endif
            do iknr = k1, k2   
              write(*,*) iknr
              if( input%properties%wannier%input .eq. "groundstate") then
                ! find G+k-vectors and eigenvectors for non-reduced k-point k
                call gengpvec( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
                call getevecfv( wf_kset%vkl( :, iknr), vgklnr, evecfv1)
                call gengpvec( wf_kset%vkl( :, wf_n_ik( n, is, iknr)), wf_kset%vkc( :, wf_n_ik( n, is, iknr)), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
                call getevecfv( wf_kset%vkl( :, wf_n_ik( n, is, iknr)), vgklnr, evecfv2)
              else if( input%properties%wannier%input .eq. "hybrid") then
                ! find G+k-vectors and eigenvectors for non-reduced k-point k
                call gengpvec( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
                call getevecfv( wf_kset%vkl( :, iknr), vgklnr, evecfv1)
                call gengpvec( wf_kset%vkl( :, wf_n_ik( n, is, iknr)), wf_kset%vkc( :, wf_n_ik( n, is, iknr)), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
                call getevecfv( wf_kset%vkl( :, wf_n_ik( n, is, iknr)), vgklnr, evecfv2)
              else if( input%properties%wannier%input .eq. "gw") then
                call getevecsvgw_new( "GW_EVECSV.OUT", iknr, wf_kset%vkl( :, iknr), nmatmax, nstfv, nspinor, evecfv1)
                call getevecsvgw_new( "GW_EVECSV.OUT", wf_n_ik( n, is, iknr), wf_kset%vkl( :, wf_n_ik( n, is, iknr)), nmatmax, nstfv, nspinor, evecfv2)
              else
                call terminate
              end if
              ! generate plane-wave matrix elements
              call emat_genemat( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), wf_fst, wf_lst, wf_fst, wf_lst, &
                   evecfv1( :, wf_fst:wf_lst, 1), &
                   evecfv2( :, wf_fst:wf_lst, 1), &
                   wf_emat( :, :, idxn, iknr))
              !call emat_genemat( wf_kset%vkl( :, iknr), wf_kset%vkc( :, iknr), wf_fst, wf_lst, wf_fst, wf_lst, &
              !     evecfv_tmp(:,:,1, iknr), &
              !     evecfv_tmp(:,:,1, wf_n_ik( n, is, iknr)), &
              !     wf_emat( :, :, n, is, iknr))
            end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
!#ifdef MPI
!          call mpi_allgatherv_ifc( wf_kset%nkpt, nstfv**2, zbuf=wf_emat( :, :, idxn, :))
!#endif
            ! make use of symmetry: M(k,-b) = M(k-b,b)**H
            dist = 1.d13
            do n2 = 1, wf_n_n( is)
              dist( n2) = norm2( wf_n_vl( :, n, is) + wf_n_vl( :, n2, is))
            end do
            n2 = minloc( dist( 1:wf_n_n( is)), 1)
            if( dist( n2) .lt. input%structure%epslat) then
              do iknr = 1, wf_kset%nkpt
                wf_emat( :, :, wf_n_ns2n( n2, is), iknr) = conjg( transpose( wf_emat( :, :, idxn, wf_n_ik( n2, is, iknr)))) 
              end do
            else
              write( *, '(" ERROR (wannier_emat): Negative neighboring vector not found for ",3F13.6)') wf_n_vl( :, n, is)
              call terminate
            end if
          end do
        end do
        deallocate( igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr)
        call emat_destroy
        call wannier_writeemat( "WANNIER")
        !deallocate( evecfv_tmp)
        call timesec( t1)
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
        write( wf_info, '(" ...plane-wave matrix-elements calculated.")')
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
        write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
      end if

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
    end subroutine wannier_emat
    !EOP

    !BOP
    ! !ROUTINE: wannier_gen_pro
    ! !INTERFACE:
    !
    subroutine wannier_gen_pro()
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
      !   Created September 2016 (Tillack)
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

      allocate( projm( wf_nst, wf_nst), sval( wf_nst), lsvec( wf_nst, wf_nst), rsvec( wf_nst, wf_nst))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))
         
#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, i, projm, sval, lsvec, rsvec)
!!$OMP DO  
#endif
      do iknr = 1, wf_kset%nkpt
        j = 0
        do i = 1, wf_nprojtot
          if( wf_projused( i) .eq. 1) then
            j = j+1
            projm( :, j) = wf_projection( wf_fst:wf_lst, j, iknr)
          end if
        end do
        call zgesdd_wrapper( projm, wf_nst, wf_nst, sval, lsvec, rsvec)
        if( minval( sval) .lt. input%structure%epslat) then
          write(*,*) ' WARNING (wannier_gen_pro): The local-orbitals selected for projection may not describe the angular character of the selected bands sufficiently. Add local-orbitals with different angular momentum.'
        end if
        call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
             lsvec, wf_nst, &
             rsvec, wf_nst, zzero, &
             wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, iknr), wf_nst)
      end do 
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
      call wannier_loc
      call wannier_writefile( 'WANNIER')
      deallocate( projm, sval, lsvec, rsvec)
      call timesec( t1)
!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      write( wf_info, '(" ...simple projection step performed.")')
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"Omega): ",T40,F13.6)') sum( wf_omega)
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
    subroutine wannier_gen_opf()
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
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      ! local variables
      integer :: iknr, is, n, nx, iproj, i, j, l, ix, minit, ntot, idxn
      real(8) :: lambda, phi, phimean, uncertainty, theta, evalmin, omegamin, limit
      real(8) :: opf_min_q(3,3), opf_min_p(3,1), opf_min_c, opf_min_sol(3), eval(3), revec(3,3), opf_min_big(6,6), opf_min_aux(3,3), tp(2)
      real(8) :: t0, t1
      complex(8) :: opf_min_z(3,1), evec(3,3), evec_big(6,6), eval_big(6), r1, r2

      ! allocatable arrays
      complex(8), allocatable :: opf_transform(:,:,:), opf_projm(:,:,:), lsvec(:,:), rsvec(:,:), opf_mixing(:,:), opf_x(:,:,:), opf_x0(:,:,:)
      complex(8), allocatable :: auxmat(:,:), projm(:,:), lsvec2(:,:), rsvec2(:,:)
      real(8), allocatable :: sval(:), opf_t(:), sval2(:), phihist(:)

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
      allocate( opf_transform( wf_nst, wf_nprojtot, wf_kset%nkpt))
      allocate( opf_projm( wf_nst, wf_nprojtot, wf_kset%nkpt))
      allocate( lsvec( wf_nst, wf_nst), rsvec( wf_nprojtot, wf_nprojtot), sval( wf_nst))
      do iknr = 1, wf_kset%nkpt
        if( wf_disentangle) then
          i = min( minval( wf_win_idxi( 1:wf_win_ni( iknr), iknr)), minval( wf_win_idxo( 1:wf_win_no( iknr), iknr)))
          j = max( maxval( wf_win_idxi( 1:wf_win_ni( iknr), iknr)), maxval( wf_win_idxo( 1:wf_win_no( iknr), iknr)))
          write(*,*) iknr, i, j, j-i+1, wf_win_no( iknr)+wf_win_ni( iknr)
          call zgemm( 'C', 'N', wf_nst, wf_nprojtot, wf_win_no( iknr)+wf_win_ni( iknr), zone, &
               wf_subspace( 1:(wf_win_no( iknr)+wf_win_ni( iknr)), :, iknr), wf_win_no( iknr)+wf_win_ni( iknr), &
               wf_projection( i:j, :, iknr), wf_win_no( iknr)+wf_win_ni( iknr), zzero, &
               opf_projm( :, :, iknr), wf_nst)
        else
          opf_projm( :, :, iknr) = wf_projection( wf_fst:wf_lst, :, iknr)
        end if
        call zgesdd_wrapper( opf_projm( :, :, iknr), wf_nst, wf_nprojtot, sval, lsvec, rsvec)
        call zgemm( 'N', 'N', wf_nst, wf_nprojtot, wf_nst, zone, &
             lsvec, wf_nst, &
             rsvec( 1:wf_nst, :), wf_nst, zzero, &
             opf_transform( :, :, iknr), wf_nst)
      end do

      lambda = 0.d0
      do i = 1, wf_n_nshells
        lambda = lambda + 0.5d0*wf_n_n( wf_n_usedshells( i))*wf_n_wgt( wf_n_usedshells( i))
      end do
      minit = nint( wf_nst*(wf_nprojtot - 0.5d0*(wf_nst + 1)))
      limit = 1.d-3

      !********************************************************************
      ! build enlarged overlap and constraint matrices
      !********************************************************************
      allocate( opf_x( wf_nprojtot, wf_nprojtot, wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_x0( wf_nprojtot, wf_nprojtot, wf_kset%nkpt*(1+wf_n_ntot)))
      allocate( opf_t( wf_kset%nkpt*(1+wf_n_ntot)))
      ! build enlarged overlap matrices
      allocate( auxmat( wf_nst, wf_nprojtot))
      i = 0
      do iknr = 1, wf_kset%nkpt
        do j = 1, wf_n_nshells
          is = wf_n_usedshells( j)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            i = i + 1
            call zgemm( 'N', 'N', wf_nst, wf_nprojtot, wf_nst, zone, &
                 wf_emat( wf_fst:wf_lst, wf_fst:wf_lst, idxn, iknr), wf_nst, &
                 opf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nst, zzero, &
                 auxmat, wf_nst)
            call zgemm( 'C', 'N', wf_nprojtot, wf_nprojtot, wf_nst, zone, &
                 opf_transform( :, :, iknr), wf_nst, &
                 auxmat, wf_nst, zzero, &
                 opf_x0( :, :, i), wf_nprojtot)
            opf_t( i) = -wf_n_wgt( is)
          end do
        end do
      end do
      ! build constraint matrices
      do iknr = 1, wf_kset%nkpt
        i = i + 1
        call zgemm( 'C', 'N', wf_nprojtot, wf_nprojtot, wf_nst, zone, &
             opf_projm( :, :, iknr), wf_nst, &
             opf_projm( :, :, iknr), wf_nst, zzero, &
             opf_x0( :, :, i), wf_nprojtot)
        do is = 1, wf_nprojtot
          opf_x0( is, is, i) = opf_x0( is, is, i) - zone
        end do
        opf_t( i) = lambda
      end do
      nx = i

      !********************************************************************
      ! minimize Lagrangian
      !********************************************************************
      if( allocated( wf_opf)) deallocate( wf_opf)
      allocate( wf_opf( wf_nprojtot, wf_fst:wf_lst))
      allocate( opf_mixing( wf_nprojtot, wf_nprojtot))
      allocate( projm( wf_nst, wf_nst), sval2( wf_nst), lsvec2( wf_nst, wf_nst), rsvec2( wf_nst, wf_nst))
      allocate( phihist( minit))
      deallocate( lsvec, rsvec)
      allocate( lsvec(3,3), rsvec(3,3))
      omegamin = 1.d13
      opf_mixing = zzero
      opf_x = opf_x0
      do i = 1, wf_nprojtot
        opf_mixing( i, i) = zone
      end do
      l = 1
      i = 1
      n = 0
      phihist = 0.d0
      uncertainty = 1.d0
      ! start minimization
      do while( (n .lt. minit) .or. ((n .lt. 500000) .and. (uncertainty .gt. limit)))
        i = 1
        do while( (i .le. wf_nst) .and. ((n .lt. minit) .or. (uncertainty .gt. limit)))
          j = i + 1
          do while( (j .le. wf_nprojtot) .and. ((n .lt. minit) .or. (uncertainty .gt. limit)))
            n = n + 1
            ! calculate Lagrangian
            phi = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, iknr) reduction(+:phi)
!$OMP DO  
#endif
            do is = 1, nx
              do iknr = 1, wf_nst
                phi = phi + opf_t( is)*opf_x( iknr, iknr, is)*conjg( opf_x( iknr, iknr, is))
              end do
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            phi = phi/wf_kset%nkpt
            do is = 1, wf_n_nshells
              phi = phi + wf_nst*wf_n_n( wf_n_usedshells( is))*wf_n_wgt( wf_n_usedshells( is))
            end do

            ! convergence analysis
            if( n .eq. 1) phihist(:) = phi
            phihist = cshift( phihist, -1)
            phihist(1) = phi
            phimean = sum( phihist(:))/minit
            if( n .gt. 1) then
              uncertainty = sqrt( sum( (phihist(:)-phimean)**2)/(min( n, minit)-1))/abs( phi)
            else
              uncertainty = 1.d0
            end if

            write(*,'(4(I7,3x),2(F23.16,3x))') n, l, i, j, phi, uncertainty

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
            if( j .le. wf_nst) then
              opf_min_sol(:) = revec( :, 1)
              if( opf_min_sol(1) .lt. 0.d0) opf_min_sol = -opf_min_sol
              !write(*,*) "Case 1"
            else
              opf_min_big = 0.d0
              opf_min_big( 1:3, 1:3) = 2*opf_min_q
              opf_min_big( 1:3, 4:6) = 0.25d0*matmul( opf_min_p, transpose( opf_min_p)) - matmul( opf_min_q, opf_min_q)
              opf_min_big( 4:6, :) = 0.d0
              opf_min_big( 4, 1) = 1.d0
              opf_min_big( 5, 2) = 1.d0
              opf_min_big( 6, 3) = 1.d0
              call diaggenmat( 6, cmplx( opf_min_big, 0, 8), eval_big, evec_big)
              evalmin = maxval( abs( eval_big))
              do is = 1, 6
                if( abs( aimag( eval_big( is))) .lt. input%structure%epslat) evalmin = min( evalmin, real( eval_big( is)))
              end do
              do is = 1, 3
                opf_min_q( is, is) = opf_min_q( is, is) - evalmin
              end do
              if( minval( eval(:) - evalmin) .lt. 1.d-6) then
                call zgesdd_wrapper( cmplx( opf_min_q, 0.d0, 8), 3, 3, eval, lsvec, rsvec)
                do is = 1, 3
                  if( eval( is) .gt. 1.d-6) then
                    rsvec( is, :) = rsvec( is, :)/eval( is)
                  else
                    rsvec( is, :) = 0.d0
                  end if
                end do
                opf_min_aux = real( matmul( transpose( rsvec), transpose( lsvec)))
                call r3mv( opf_min_aux, -0.5d0*opf_min_p(:,1), opf_min_sol)
                !write(*,*) "Case 2"
                call r3mv( opf_min_q, opf_min_sol, opf_min_aux(:,1))
                if( (norm2( opf_min_aux(:,1) + 0.5d0*opf_min_p(:,1)) .ge. 1.d-6) .or. (norm2( opf_min_sol) .gt. 1.d0)) then
                  opf_min_sol = 10.d0
                else
                  opf_min_sol = opf_min_sol + revec(:,1)/norm2( revec(:,1))*sqrt( 1.d0 - norm2( opf_min_sol))
                end if
              else
                !write(*,*) "Case 3"
                call r3minv( opf_min_q, opf_min_aux)
                call r3mv( opf_min_aux, -0.5d0*opf_min_p(:,1), opf_min_sol)
              end if
            end if
            if( abs( 1.d0 - norm2( opf_min_sol)) .le. input%structure%epslat) then
              call sphcrd( (/opf_min_sol(2), opf_min_sol(3), opf_min_sol(1)/), r1, tp)
              phi = tp(2)
              theta = tp(1)/2
              r1 = cmplx( cos( theta), 0, 8)
              r2 = exp( zi*phi)*sin( theta)
              ! update mixing matrix
              auxmat(1,:) = opf_mixing(:,i)
              auxmat(2,:) = opf_mixing(:,j)
              opf_mixing(:,i) = auxmat(1,:)*r1 - auxmat(2,:)*conjg( r2)
              opf_mixing(:,j) = auxmat(2,:)*r1 + auxmat(1,:)*r2
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, auxmat)
!$OMP DO  
#endif
              ! update X-matrices
              do is = 1, nx
                auxmat(1,:) = opf_x( :, i, is)
                auxmat(2,:) = opf_x( :, j, is)
                opf_x( :, i, is) = auxmat(1,:)*r1 - auxmat(2,:)*conjg( r2)
                opf_x( :, j, is) = auxmat(2,:)*r1 + auxmat(1,:)*r2
                auxmat(1,:) = opf_x( i, :, is)
                auxmat(2,:) = opf_x( j, :, is)
                opf_x( i, :, is) = auxmat(1,:)*r1 - auxmat(2,:)*r2
                opf_x( j, :, is) = auxmat(2,:)*r1 + auxmat(1,:)*conjg( r2)
              end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            else
              !write(*,*) "Doh!"
            end if
            j = j + 1
          end do
          i = i + 1
        end do
        l = l + 1
      end do
      wf_opf( :, wf_fst:wf_lst) = opf_mixing( :, 1:wf_nst)

      !********************************************************************
      ! generate transformation matrices
      !********************************************************************
      if( allocated( wf_transform)) deallocate( wf_transform)
      if( wf_disentangle) then
        allocate( wf_transform( maxval( wf_win_no+wf_win_ni), wf_nst, iknr))
      else
        allocate( wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, iknr))
      end if
      deallocate( auxmat)
      allocate( auxmat( wf_nst, wf_nst))
      do iknr = 1, wf_kset%nkpt
        call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nprojtot, zone, &
             opf_projm( :, :, iknr), wf_nst, &
             wf_opf( :, wf_fst:wf_lst), wf_nprojtot, zzero, &
             projm, wf_nst)
        call zgesdd_wrapper( projm, wf_nst, wf_nst, sval2, lsvec2, rsvec2)
        call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
             lsvec2, wf_nst, &
             rsvec2, wf_nst, zzero, &
             auxmat( :, :), wf_nst)
        if( wf_disentangle) then
          call zgemm( 'N', 'N', wf_win_no( iknr)+wf_win_ni( iknr), wf_nst, wf_nst, zone, &
               wf_subspace( 1:(wf_win_no( iknr)+wf_win_ni( iknr)), :, iknr), wf_win_no( iknr)+wf_win_ni( iknr), &
               auxmat, wf_nst, zzero, &
               wf_transform( 1:(wf_win_no( iknr)+wf_win_ni( iknr)), :, iknr), wf_win_no( iknr)+wf_win_ni( iknr))
        else
          wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, iknr) = auxmat(:,:)
        end if
      end do
      call wannier_loc
      call wannier_writefile( 'WANNIER')
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

    !BOP
    ! !ROUTINE: wannier_gen_max
    ! !INTERFACE:
    !
    subroutine wannier_maxloc()
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
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      ! local variables
      integer :: ix, iy, iz, i, n, is, iknr, ngknr, k1, k2, z0, z1, idxn
      real(8) :: mixing, alpha, grad, dg1, dg2, nwgt
      real(8) :: omegastart, omega, omegamean, uncertainty, omega2
      real(8) :: t0, t1, v1(3), v2(3)
      real(8), allocatable :: omegai(:), omegad(:), omegaod(:)
      integer :: minit, maxit
      logical :: succes

      ! allocatable arrays
      integer, allocatable :: nn(:), integerlist(:)
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      real(8), allocatable :: ndist(:), nvl(:,:,:), nvc(:,:,:), bwgt(:), ravg(:,:)
      real(8), allocatable :: omegahist(:), gradhist(:), eval(:)
      complex(8), allocatable :: auxmat(:,:), mlwf_m0(:,:,:,:), mlwf_m(:,:,:,:), evec(:,:), mlwf_r(:,:), mlwf_t(:,:), mlwf_dw(:,:), mlwf_transform(:,:,:)
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)
      complex(8), allocatable :: auxmatread(:,:)

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
        allocate( auxmat( wf_nst, wf_nst))
        allocate( mlwf_m0( wf_nst, wf_nst, wf_n_ntot, wf_kset%nkpt))
        allocate( omegai( wf_nst))

        write( wf_info, '(" minimize localization functional Omega...")')
        call timesec( t0)
        omegai = 0.d0

      !********************************************************************
      ! generate M0 matrices
      !********************************************************************

        do i = 1, wf_n_nshells
          is = wf_n_usedshells( i)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            write(*,*) idxn
            k1 = 1
            k2 = wf_kset%nkpt
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ix, auxmat) reduction(+:omegai)
!$OMP DO  
#endif
            do iknr = k1, k2   
                               
              call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
                   wf_emat( wf_fst:wf_lst, wf_fst:wf_lst, idxn, iknr), wf_nst, &
                   wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, wf_n_ik( n, is, iknr)), wf_nst, zzero, &
                   auxmat, wf_nst)
              call ZGEMM( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
                   wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, iknr), wf_nst, &
                   auxmat, wf_nst, zzero, &
                   mlwf_m0( 1:wf_nst, 1:wf_nst, idxn, iknr), wf_nst)
  
              ! independent part of localization functional
              do ix = 1, wf_nst
                omegai( ix) = omegai( ix) + wf_n_wgt( is)*(1.d0 - sum( abs( mlwf_m0( :, ix, idxn, iknr))**2))
              end do

            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
          end do
        end do
        
        omegai = omegai/wf_kset%nkpt
        
      !********************************************************************
      ! minimize localization functional
      !********************************************************************
        ! mixing parameter for self consistent minimization
        mixing = dble( 0.5)
        alpha = mixing
        ! minimum/maximum number of iterations
        minit = 200
        maxit = 20000

        ! initialize transformation matrices
        allocate( mlwf_transform( wf_nst, wf_nst, wf_kset%nkpt))
        mlwf_transform = zzero
        do ix = 1, wf_nst
          mlwf_transform( ix, ix, :) = zone
        end do

        ! start minimization loop
        allocate( evec( wf_nst, wf_nst), eval( wf_nst))
        allocate( mlwf_r( wf_nst, wf_nst), &
                  mlwf_t( wf_nst, wf_nst), &
                  mlwf_dw( wf_nst, wf_nst))
        allocate( mlwf_m( wf_nst, wf_nst, wf_n_ntot, wf_kset%nkpt))
        allocate( ravg( 3, wf_nst))
        allocate( omegahist( minit), gradhist( minit), omegaod( wf_nst), omegad( wf_nst))
        allocate( integerlist( minit))
        write( *, '("shape of M: ",100I)') shape( mlwf_m)
        do iz = 1, minit
          integerlist( iz) = iz
        end do
        mlwf_m = mlwf_m0
        iz = 0
        omegastart = sum( wf_omega)
        write(*,*) omegastart
        omegahist = 0.d0
        gradhist = 1.d0
        grad = 1.d0
        succes = .false.

        nwgt = 0.d0
        do i = 1, wf_n_nshells
          is = wf_n_usedshells( i)
          nwgt = nwgt + wf_n_wgt( is)*wf_n_n( is)
        end do
        write(*, '("nwgt: ",F13.6)') nwgt
      !************************************************************
      !* start of minimization
      !************************************************************
        do while( .not. succes)
          iz = iz + 1
          ! centers of Wannier functions
          ravg = 0.d0
          do ix = 1, wf_nst
            do i = 1, wf_n_nshells
              is = wf_n_usedshells( i)
              do n = 1, wf_n_n( is)
                idxn = wf_n_ns2n( n, is)
                ravg( :, ix) = ravg( :, ix) - wf_n_wgt( is)/wf_kset%nkpt*sum( real( aimag( log( mlwf_m( ix, ix, idxn, :)))))*wf_n_vc( :, n, is)
              end do
            end do
          end do
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, i, is, n, idxn, mlwf_r, mlwf_t, mlwf_dw, eval, evec, ix, auxmat)
!$OMP DO
#endif
          do iknr = 1, wf_kset%nkpt
            mlwf_dw(:,:) = zzero
            do i = 1, wf_n_nshells
              is = wf_n_usedshells( i)
              do n = 1, wf_n_n( is)
                idxn = wf_n_ns2n( n, is)
                ! calculating R and T
                do ix = 1, wf_nst
                  mlwf_r( :, ix) = mlwf_m( :, ix, idxn, iknr)*conjg( mlwf_m( ix, ix, idxn, iknr))
                  mlwf_t( :, ix) = mlwf_m( :, ix, idxn, iknr)/mlwf_m( ix, ix, idxn, iknr)*(real( aimag( log( mlwf_m( ix, ix, idxn, iknr)))) + dot_product( wf_n_vc( :, n, is), ravg( :, ix)))
                end do
        
                ! calculating dW
                mlwf_r = mlwf_r - conjg( transpose( mlwf_r))
                mlwf_t = mlwf_t + conjg( transpose( mlwf_t))
                mlwf_dw(:,:) = mlwf_dw(:,:) + wf_n_wgt( is)*( 0.5*mlwf_r(:,:) + 0.5*zi*mlwf_t(:,:))
              end do
            end do
            mlwf_dw(:,:) = mlwf_dw(:,:)*alpha/nwgt

            ! updating transformation matrices
            call diaghermat( wf_nst, zi*mlwf_dw, eval, evec)
            do ix = 1, wf_nst
              mlwf_dw( :, ix) = exp( -zi*eval( ix))*evec( :, ix) 
            end do
            call ZGEMM( 'N', 'C', wf_nst, wf_nst, wf_nst, zone, &
                 mlwf_dw, wf_nst, &
                 evec, wf_nst, zzero, &
                 auxmat, wf_nst)
            call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
                 mlwf_transform( :, :, iknr), wf_nst, &
                 auxmat, wf_nst, zzero, &
                 mlwf_dw, wf_nst)
            mlwf_transform( :, :, iknr) = mlwf_dw(:,:)
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
    
          ! gauge-dependent part of localization functional
          omegaod = 0.d0
          omegad = 0.d0
          omega2 = 0.d0

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, i, n, is, idxn, ix, auxmat) reduction(+:omegaod, omegad, omega2)
!$OMP DO  
#endif
          ! updating M
          do iknr = 1, wf_kset%nkpt
            do i = 1, wf_n_nshells 
              is = wf_n_usedshells( i)
              do n = 1, wf_n_n( is)
                idxn = wf_n_ns2n( n, is)
                call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
                     mlwf_m0( :, :, idxn, iknr), wf_nst, &
                     mlwf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nst, zzero, &
                     auxmat, wf_nst)
                call ZGEMM( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
                     mlwf_transform( :, :, iknr), wf_nst, &
                     auxmat, wf_nst, zzero, &
                     mlwf_m( :, :, idxn, iknr), wf_nst)
                do ix = 1, wf_nst
                  omegad( ix) = omegad( ix) + wf_n_wgt( is)*(real( aimag( log( mlwf_m( ix, ix, idxn, iknr)))) &
                                  + dot_product( wf_n_vc( :, n, is), ravg( :, ix)))**2
                  do iy = 1, wf_nst
                    if( iy .ne. ix) then
                      omegaod( ix) = omegaod( ix) + wf_n_wgt( is)*abs( mlwf_m( iy, ix, idxn, iknr))**2
                    end if
                  end do
                  omega2 = omega2 + wf_n_wgt( is)*(1.d0 - abs( mlwf_m( ix, ix, idxn, iknr))**2 + &
                      & ( real( aimag( log( mlwf_m( ix, ix, idxn, iknr))))&
                      &   + dot_product( wf_n_vc( :, n, is), ravg( :, ix)))**2 )
                end do
              end do
            end do
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
          ! convergence analysis
          omegaod = omegaod/wf_kset%nkpt
          omegad = omegad/wf_kset%nkpt
          omega2 = omega2/wf_kset%nkpt

          omega = sum( omegai) + sum( omegad) + sum( omegaod)
          !omega = omega2

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

          if( (minval( v1)*maxval( v1) .ge. 0.d0) .and. &
              (minval( v2)*maxval( v2) .ge. 0.d0) .and. &
              (minval( v1)+maxval( v1) .ge. 0.d0) .and. &
              (minval( v2)+maxval( v2) .le. 0.d0)) succes = .true.
          !if( omega .gt. omegastart) succes = .true.
          if( uncertainty .le. input%properties%wannier%uncertainty) succes = .true.
          if( maxval( gradhist) .le. input%properties%wannier%uncertainty) succes = .true.
          if( iz .lt. minit) succes = .false.
          if( iz .ge. maxit) succes = .true.
          write(*,'(I4,10F23.16)') iz, omega, uncertainty, grad, maxval( gradhist)!, omega2 !omegai+omegad+omegaod
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
          write(*,'(" SUCCES: Convergence reached after ",I4," cycles.")') iz 
          write(*,'(" Localization gain: ",I3,"%")') nint( 100d0*(omegastart-omega)/omega)
        end if
        
      !********************************************************************
      ! generate transformation matrices
      !********************************************************************
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, auxmat)
!$OMP DO  
#endif
        do iknr = 1, wf_kset%nkpt
          call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
               wf_transform( wf_fst:wf_lst, :, iknr), wf_nst, &
               mlwf_transform( :, :, iknr), wf_nst, zzero, &
               auxmat, wf_nst)
          wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, iknr) = auxmat(:,:)
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        call wannier_loc
        call wannier_writefile( 'WANNIER')
  
        deallocate( auxmat, eval, evec, mlwf_m0, mlwf_m, mlwf_r, mlwf_t, mlwf_dw, mlwf_transform, omegahist, gradhist)
        call timesec( t1)
        write( wf_info, '(" localization functional Omega minimized.")')
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"minimum/maximum iterations: ",T40,I6,"/",I6)') minit, maxit
        write( wf_info, '(5x,"iterations: ",T40,7x,I6)') iz
        write( wf_info, '(5x,"aimed uncertainty: ",T40,E13.6)') input%properties%wannier%uncertainty
        write( wf_info, '(5x,"achieved uncertainty: ",T40,E13.6)') uncertainty
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
    
    subroutine wannier_gen_fromfile()
      logical :: succes
      real(8), allocatable :: ravg(:,:), omega(:)

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call wannier_readfile( 'WANNIER', succes)
      if( .not. succes) then
        write( wf_info, '(" Failed to read Wannier functions from file. Aborted.")')
        write( wf_info, *)
        call terminate
      else
        write( wf_info, '(" file successfully read")')
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
    
    subroutine wannier_subspace( icenter)
      integer, intent( in) :: icenter

!      integer :: iknr, i, j, s, n, nk, bk, nkb, bkb
!      real(8) :: mixing, alpha
!      logical :: converged
!      complex(8), allocatable :: c0(:,:,:), c1(:,:,:), mdiagi(:,:,:), mdiago(:,:,:), emat(:,:,:,:,:), auxmat(:,:), auxmat2(:,:), lsvec(:,:), rsvec(:,:)
!      real(8), allocatable :: eval(:,:), sval(:)
!      
!      write(*,*) icenter
!      ! write reduced matrix-elements
!      if( .not. allocated( wf_emat)) call wannier_emat
!      allocate( emat( maxval( wf_win_no), maxval( wf_win_no), maxval( wf_n_n( 1:wf_n_nshells)), wf_n_nshells, wf_kset%nkpt))
!      allocate( lsvec( maxval( wf_win_no), maxval( wf_win_no)), rsvec( wf_nst, wf_nst), sval( wf_nst))
!      allocate( c0( maxval( wf_win_no), maxval( wf_win_no), wf_kset%nkpt), &
!                c1( maxval( wf_win_no), maxval( wf_win_no), wf_kset%nkpt))
!      allocate( mdiagi( maxval( wf_win_no), maxval( wf_win_no), wf_kset%nkpt))
!      allocate( mdiago( maxval( wf_win_no), maxval( wf_win_no), wf_kset%nkpt))
!      allocate( eval( maxval( wf_win_no), wf_kset%nkpt))
!      allocate( auxmat2( maxval( wf_win_no), wf_nprojtot))
!
!      if( allocated( wf_subspace)) deallocate( wf_subspace)
!      s = maxval( wf_win_no + wf_win_ni)
!      allocate( wf_subspace( s, wf_nst, wf_kset%nkpt))
!
!      wf_disentangle = .false.  ! find OPF without disentanglement first
!      call wannier_gen_opf( nint( icenter - 0.5d0*wf_nst), wf_nst)
!      emat = zzero
!      c0 = zzero
!      c1 = zzero
!      do iknr = 1, wf_kset%nkpt
!        do i = 1, wf_win_no( iknr)
!          do s = 1, wf_n_nshells
!            do n = 1, wf_n_n( s)
!              do j = 1, wf_win_no( wf_n_ik( n, s, iknr))
!                emat( i, j, n, s, iknr) = wf_emat( wf_win_idxo( i, iknr), wf_win_idxo( j, wf_n_ik( n, s, iknr)), n, s, iknr)
!              end do
!            end do
!          end do
!          auxmat2( i, :) = wf_projection( wf_win_idxo( i, iknr), :, iknr)
!        end do
!        call zgemm( 'N', 'N', wf_win_no( iknr), wf_nst, wf_nprojtot, zone, &
!             auxmat2( 1:wf_win_no( iknr), :), wf_win_no( iknr), &
!             wf_opf, wf_nprojtot, zzero, &
!             c0( 1:wf_win_no( iknr), 1:wf_nst, iknr), wf_win_no( iknr))
!        call zgesdd_wrapper( c0( 1:wf_win_no( iknr), 1:wf_nst, iknr), wf_win_no( iknr), wf_nst, sval, lsvec( 1:wf_win_no( iknr), 1:wf_win_no( iknr)), rsvec)
!        call ZGEMM( 'N', 'N', wf_win_no( iknr), wf_nst, wf_nst, zone, &
!             lsvec( 1:wf_win_no( iknr), 1:wf_nst), wf_win_no( iknr), &
!             rsvec, wf_nst, zzero, &
!             c1( 1:wf_win_no( iknr), (wf_win_no( iknr)-(wf_nst-wf_win_ni( iknr))+1):wf_win_no( iknr), iknr), wf_win_no( iknr))
!        !write(*,*) iknr
!        !call plotmat( c1( :, :, iknr))
!      end do
!
!      deallocate( lsvec, rsvec, sval)
!      allocate( auxmat( maxval( wf_win_no), maxval( wf_nst - wf_win_ni)))
!      allocate( sval( maxval( wf_win_no)))
!      converged = .false.
!      j = 0
!      mixing = 0.5d0
!      alpha = 1.d0
!      mdiagi = zzero
!      do while( (.not. converged) .and. (j .lt. 100000))
!        j = j + 1
!        converged = .true.
!        !write(*,*) j
!        c0 = zzero
!        do iknr = 1, wf_kset%nkpt
!          mdiagi(:,:,iknr) = (1.d0-alpha)*mdiagi(:,:,iknr)
!          nk = wf_win_no( iknr)
!          bk = nk - (wf_nst - wf_win_ni( iknr)) + 1
!          !write(*,*) nk, bk
!          do s = 1, wf_n_nshells
!            do n = 1, wf_n_n( s)
!              nkb = wf_win_no( wf_n_ik( n, s, iknr))
!              bkb = nkb - (wf_nst - wf_win_ni( wf_n_ik( n, s, iknr))) + 1
!              call zgemm( 'N', 'N', nk, nkb-bkb+1, nkb, zone, &
!                   emat( 1:nk, 1:nkb, n, s, iknr), nk, &
!                   c1( 1:nkb, bkb:nkb, wf_n_ik( n, s, iknr)), nkb, zzero, &
!                   auxmat( 1:nk, 1:(nkb-bkb+1)), nk) 
!              call zgemm( 'N', 'C', nk, nk, nkb-bkb+1, cmplx( alpha*wf_n_wgt( s), 0, 8), &
!                   auxmat( 1:nk, 1:(nkb-bkb+1)), nk, &
!                   auxmat( 1:nk, 1:(nkb-bkb+1)), nk, zone, &
!                   mdiagi( 1:nk, 1:nk, iknr), nk)
!            end do
!          end do
!          !mdiagi(:,:,iknr) = mdiagi(:,:,iknr) + (1.d0-alpha)*mdiago(:,:,iknr)
!          !c0( :, :, iknr) = zzero
!          sval( 1:nk) = eval( 1:nk, iknr)
!          call diaghermat( nk, mdiagi( 1:nk, 1:nk, iknr), eval( 1:nk, iknr), c0( 1:nk, 1:nk, iknr))
!          if( iknr .eq. 42) then
!            !write(*,*) bk, nk
!            !write(*,'(100F23.16)') eval( 1:nk, iknr)
!            !call plotmat( c0( 1:i, 1:(wf_nst - wf_win_ni( iknr)), iknr))
!            write(*,'(I6,F23.16)') j, norm2( eval( 1:nk, iknr) - sval( 1:nk))
!          end if
!          if( norm2( eval( 1:nk, iknr) - sval( 1:nk)) .gt. 1.d-6) converged = .false.
!          !mdiago(:,:,iknr) = mdiagi(:,:,iknr)
!        end do
!        c1 = c0
!        alpha = mixing
!      end do
!
!      ! write subspaces
!      wf_subspace = zzero
!      do iknr = 1, wf_kset%nkpt
!        bk = minval( wf_win_idxo( 1:wf_win_no( iknr), iknr))
!        if( wf_win_ni( iknr) .gt. 0) then
!          nk = minval( wf_win_idxi( 1:wf_win_ni( iknr), iknr))
!          nkb = maxval( wf_win_idxi( 1:wf_win_ni( iknr), iknr)) + 1
!        else
!          nk = maxval( wf_win_idxo( 1:wf_win_no( iknr), iknr)) + 1
!          nkb = nk
!        end if
!        write(*,*) iknr, wf_win_no( iknr), wf_win_ni( iknr)!, bk, nk, nkb
!        s = 0
!        do j = 1, wf_nst
!          if( any( wf_win_idxi( 1:wf_win_ni( iknr), iknr) .eq. j+bk-1, 1)) then
!            !write(*,'("in inner: ",3I4)') j, j+bk-1, iknr
!            wf_subspace( j, j, iknr) = zone
!          else
!            !write(*,'("in outer: ",3I4)') j, j+bk-1, iknr
!            s = s + 1
!            do i = 1, nk-bk
!              !write(*,*) i, wf_win_no( iknr)-(wf_nst-wf_win_ni( iknr))+s
!              wf_subspace( i, j, iknr) = c1( i, wf_win_no( iknr)-(wf_nst-wf_win_ni( iknr))+s, iknr)
!            end do
!            do i = nkb-bk+1, wf_win_no( iknr)+wf_win_ni( iknr)
!              !write(*,*) i, wf_win_no( iknr)-(wf_nst-wf_win_ni( iknr))+s
!              wf_subspace( i, j, iknr) = c1( i-(nkb-nk), wf_win_no( iknr)-(wf_nst-wf_win_ni( iknr))+s, iknr)
!            end do
!          end if
!        end do 
!        call plotmat( wf_subspace( 1:(wf_win_no( iknr)+wf_win_ni( iknr)), :, iknr))
!        write(*,*)
!      end do
!
!      deallocate( c0, c1, mdiagi, mdiago, eval, emat, auxmat, auxmat2)
!      wf_disentangle = .true.
    end subroutine wannier_subspace

    subroutine wannier_loc()
      integer :: iknr, is, n, i, j, k1, k2, idxn
      complex(8), allocatable :: loc_m(:,:,:,:), auxmat(:,:)
      
      allocate( auxmat( wf_fst:wf_lst, wf_fst:wf_lst))
      allocate( loc_m( wf_fst:wf_lst, wf_fst:wf_lst, wf_n_ntot, wf_kset%nkpt))

      if( .not. wf_initialized) call wannier_init
      if( .not. allocated( wf_emat)) call wannier_emat
      if( .not. allocated( wf_centers)) allocate( wf_centers( 3, wf_fst:wf_lst))
      if( .not. allocated( wf_omega)) allocate( wf_omega( wf_fst:wf_lst))

#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, is, i, n, idxn, auxmat)
!!$OMP DO
#endif
      do iknr = 1, wf_kset%nkpt
        do i = 1, wf_n_nshells 
          is = wf_n_usedshells( i)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            call ZGEMM( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
                 wf_emat( wf_fst:wf_lst, wf_fst:wf_lst, idxn, iknr), wf_nst, &
                 wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, wf_n_ik( n, is, iknr)), wf_nst, zzero, &
                 auxmat( wf_fst:wf_lst, wf_fst:wf_lst), wf_nst)
            call ZGEMM( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
                 wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, iknr), wf_nst, &
                 auxmat( 1:wf_nst, :), wf_nst, zzero, &
                 loc_m( wf_fst:wf_lst, wf_fst:wf_lst, idxn, iknr), wf_nst)
          end do
        end do
      end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
      wf_centers = 0.d0
      do j = wf_fst, wf_lst
        do i = 1, wf_n_nshells
          is = wf_n_usedshells( i)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            wf_centers( :, j) = wf_centers( :, j) - wf_n_wgt( is)/wf_kset%nkpt*sum( real( aimag( log( loc_m( j, j, idxn, :)))))*wf_n_vc( :, n, is)
          end do
        end do
      end do
      wf_omega = 0.d0
#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, is, n, i, idxn, j) reduction(+:wf_omega)
!!$OMP DO
#endif
      do iknr = 1, wf_kset%nkpt
        do i = 1, wf_n_nshells 
          is = wf_n_usedshells( i)
          do n = 1, wf_n_n( is)
            idxn = wf_n_ns2n( n, is)
            do j = wf_fst, wf_lst
              !write( *, '(3(3F13.6,3x))') wf_n_vc( :, n, is), wf_centers( :, j), dot_product( wf_n_vc( :, n, is), wf_centers( :, j)), real( aimag( log( loc_m( j, j, idxn, iknr)))), abs( loc_m( j, j, idxn, iknr))
              wf_omega( j) = wf_omega( j) + wf_n_wgt( is)*(1.d0 - abs( loc_m( j, j, idxn, iknr))**2 + &
                  & ( real( aimag( log( loc_m( j, j, idxn, iknr))))&
                  &   + dot_product( wf_n_vc( :, n, is), wf_centers( :, j)))**2 )
            end do
          end do
        end do
      end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
      wf_omega = wf_omega/wf_kset%nkpt

      deallocate( loc_m, auxmat)
      return
    end subroutine wannier_loc

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
    
    ! reads transformation matrices from file
    subroutine wannier_readfile( filename, succes)
      character(*), intent( in) :: filename
      logical, intent( out) :: succes

      ! local variables
      integer :: ik, ix, iy, iz, un
      integer :: fst_, lst_, nst_, nkpt_
      real(8) :: vkl_( 3), vkl_tmp( 3, wf_kset%nkpt)

      call getunit( un)

      succes = .true.
      inquire( file=trim( filename)//trim( filext), exist=succes)
      if( .not. succes) then
        write(*,*) 'ERROR (wannier_readfile): File '//trim( filename)//trim( filext)//' does not exist.'
        return
      end if
      open( un, file=trim( filename)//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
      read( un) fst_, lst_, nst_, nkpt_
      if( (fst_ .ne. wf_fst) .or. (lst_ .ne. wf_lst)) then
        write( *, '(" ERROR (wannier_readfile): different band-ranges in input (",I4,":",I4,") and file (",I4,":",I4,").")'), wf_fst, wf_lst, fst_, lst_
        call terminate
      end if
      if( nkpt_ .ne. wf_kset%nkpt) then
        write( *, '(" ERROR (wannier_readfile): different number of k-points in input (",I4,") and file (",I4,").")'), wf_kset%nkpt, nkpt_
        call terminate
      end if
      if( allocated( wf_omega)) deallocate( wf_omega)
      allocate( wf_omega( wf_fst:wf_lst))
      if( allocated( wf_centers)) deallocate( wf_centers)
      allocate( wf_centers( 3, wf_fst:wf_lst))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))
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
        do iy = wf_fst, wf_lst
          do ix = wf_fst, wf_lst
            read( un) wf_transform( ix, iy, iz)
          end do
        end do
      end do
      do ix = 1, wf_nprojtot
        read( un) wf_projused( ix)
      end do
      do ix = wf_fst, wf_lst
        read( un) wf_omega( ix)
      end do
      do iy = wf_fst, wf_lst
        do ix = 1, 3
          read( un) wf_centers( ix, iy)
        end do
      end do
      close( un)
      if( succes) write(*,*) 'Transformation matrices succesfully read.'
      return
    end subroutine wannier_readfile
    
    ! writes transformation matrices to file
    subroutine wannier_writefile( filename)
      character(*), intent( in) :: filename
      ! local variables
      integer :: ik, ix, iy, un
      
      call getunit( un)

      open( un, file=trim( filename)//trim( filext), action='WRITE', form='UNFORMATTED')
      write( un) wf_fst, wf_lst, wf_nst, wf_kset%nkpt
      do ik = 1, wf_kset%nkpt
        write( un) wf_kset%vkl( :, ik)
        do iy = wf_fst, wf_lst
          do ix = wf_fst, wf_lst
            write( un) wf_transform( ix, iy, ik)
          end do
        end do
      end do
      do ix = 1, wf_nprojtot
        write( un) wf_projused( ix)
      end do
      do ix = wf_fst, wf_lst
        write( un) wf_omega( ix)
      end do
      do iy = wf_fst, wf_lst
        do ix = 1, 3
          write( un) wf_centers( ix, iy)
        end do
      end do
      close( un)
      write( *, '(a,a)') ' Transformation matrices written to file ', trim( filename)//trim( filext)
      return
    end subroutine wannier_writefile  

    subroutine wannier_reademat( filename, succes)
      character(*), intent( in) :: filename
      logical, intent( out) :: succes

      ! local variables
      integer :: i, j, is, n, idxn, ik, ix, iy, iz, un
      integer :: fst_, lst_, nst_, ntot_, nkpt_
      real(8) :: vln(3), vkl_( 3), vkl_tmp( 3, wf_kset%nkpt)

      call getunit( un)

      succes = .true.
      inquire( file=trim( filename)//"_EMAT"//trim( filext), exist=succes)
      if( .not. succes) then
        write(*,*) 'ERROR (wannier_reademat): File '//trim( filename)//"_EMAT"//trim( filext)//' does not exist.'
        return
      end if
      open( un, file=trim( filename)//"_EMAT"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
      read( un) fst_, lst_, nst_, ntot_, nkpt_
      if( (fst_ .ne. wf_fst) .or. (lst_ .ne. wf_lst)) then
        write( *, '(" ERROR (wannier_reademat): different band-ranges in input (",I4,":",I4,") and file (",I4,":",I4,").")'), wf_fst, wf_lst, fst_, lst_
        call terminate
      end if
      if( ntot_ .ne. wf_n_ntot) then
        write( *, '(" ERROR (wannier_reademat): different number of BZ-neighbors in input (",I4,") and file (",I4,").")'), wf_n_ntot, ntot_
        call terminate
      end if
      if( nkpt_ .ne. wf_kset%nkpt) then
        write( *, '(" ERROR (wannier_reademat): different number of k-points in input (",I4,") and file (",I4,").")'), wf_kset%nkpt, nkpt_
        call terminate
      end if
      if( allocated( wf_emat)) deallocate( wf_emat)
      allocate( wf_emat( wf_fst:wf_lst, wf_fst:wf_lst, wf_n_ntot, wf_kset%nkpt))
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
          if( idxn .gt. 0) then
            do ik = 1, wf_kset%nkpt
              read( un) vkl_
              vkl_tmp( 1, :) = wf_kset%vkl( 1, :) - vkl_( 1)
              vkl_tmp( 2, :) = wf_kset%vkl( 2, :) - vkl_( 2)
              vkl_tmp( 3, :) = wf_kset%vkl( 3, :) - vkl_( 3)
              iz = minloc( norm2( vkl_tmp( :, :), 1), 1)
              if( norm2( vkl_tmp( :, iz)) .gt. input%structure%epslat) then
                write( *, '(" ERROR (wannier_reademat): k-point in file not in k-point-set.")')
                write( *, '(3F23.6)') vkl_
                call terminate
              end if
              do iy = wf_fst, wf_lst
                do ix = wf_fst, wf_lst
                  read( un) wf_emat( ix, iy, idxn, ik)
                end do
              end do
            end do
          else
            write( *, '(" ERROR (wannier_reademat): neighboring vector in file not consistens with input.")')
            write( *, '(3F23.6)') vln
          end if
        end do
      end do
      close( un)
      if( succes) write(*,*) 'Plane-wave matrix-elements succesfully read.'
      return
    end subroutine wannier_reademat
    
    ! writes transformation matrices to file
    subroutine wannier_writeemat( filename)
      character(*), intent( in) :: filename
      ! local variables
      integer :: i, ik, is, n, idxn, ix, iy, un
      
      call getunit( un)

      open( un, file=trim( filename)//"_EMAT"//trim( filext), action='WRITE', form='UNFORMATTED')
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
                write( un) wf_emat( ix, iy, idxn, ik)
              end do
            end do
          end do
        end do
      end do
      close( un)
      write( *, '(a,a)') ' Plane-wave matrix-elements written to file ', trim( filename)//"_EMAT"//trim( filext)
      return
    end subroutine wannier_writeemat  

    !BOP
    ! !ROUTINE: wannier_geometry
    ! !INTERFACE:
    !
    subroutine wannier_geometry
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      ! !DESCRIPTION:
      !   Sets geometry environment for Wannier functions. Finds near neighbors
      !   and corresponding geometric weights as well as ne number of shells
      !   needed for gradient calculation.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      integer :: ix, iy, iz, i, j, nshell, nkmax, iknr, d
      integer :: mt( wf_kset%nkpt), vi(3)
      real(8) :: vl(3), vc(3), dist
      real(8) :: nvlt( 3, wf_kset%nkpt, wf_kset%nkpt), nvct( 3, wf_kset%nkpt, wf_kset%nkpt), dt( wf_kset%nkpt)
      real(8) :: coeff(6,6), coeffcpy(6,6), right(6), sval(6), lsvec(6,6), rsvec(6,6), coeffinv(6,6)
      real(8) :: mat1(3,3), mat2(3,3), mat3(3,3)
      logical :: stopshell

      integer :: lwork, info
      real(8), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      
      ! find possible distances (shells)
      dt = 0.d0
      mt = 0
      nvlt = 0.d0
      nvct = 0.d0
      i = 0
      do iz = wf_kset%ngridk(3)-1, -wf_kset%ngridk(3)+1, -1
        do iy = wf_kset%ngridk(2)-1, -wf_kset%ngridk(2)+1, -1
          do ix = wf_kset%ngridk(1)-1, -wf_kset%ngridk(1)+1, -1
            if( (abs(ix)+abs(iy)+abs(iz)) .ne. 0) then
              vl = (/dble( ix)/wf_kset%ngridk(1), dble( iy)/wf_kset%ngridk(2), dble( iz)/wf_kset%ngridk(3)/)
              call r3mv( bvec, vl, vc)
              dist = norm2( vc)
              if( minval( abs( dt(:) - dist)) .gt. input%structure%epslat) then
                i = i + 1
                dt( i) = dist
              end if
            end if
          end do
        end do
      end do
      
      ! find all possible neighbors
      nshell = i
      nkmax = 0
      do iz = wf_kset%ngridk(3)-1, -wf_kset%ngridk(3)+1, -1
        do iy = wf_kset%ngridk(2)-1, -wf_kset%ngridk(2)+1, -1
          do ix = wf_kset%ngridk(1)-1, -wf_kset%ngridk(1)+1, -1
            if( (abs(ix)+abs(iy)+abs(iz)) .ne. 0) then
              vl = (/dble( ix)/wf_kset%ngridk(1), dble( iy)/wf_kset%ngridk(2), dble( iz)/wf_kset%ngridk(3)/)
              call r3mv( bvec, vl, vc)
              dist = norm2( vc)
              j = minloc( abs( dt(:) - dist), 1)
              if( abs( dt( j) - dist) .lt. input%structure%epslat) then
                mt( j) = mt( j) + 1
                nkmax = max( nkmax, mt(j))
                nvlt( :, mt( j), j) = vl
                nvct( :, mt( j), j) = vc
              end if
            end if
          end do
        end do
      end do
      
      ! allocating geometry arrays
      if( .not. allocated( wf_n_dist))          allocate( wf_n_dist( nshell))
      if( .not. allocated( wf_n_n))             allocate( wf_n_n( nshell))
      if( .not. allocated( wf_n_ik))            allocate( wf_n_ik( nkmax, nshell, wf_kset%nkpt))
      if( .not. allocated( wf_n_vl))            allocate( wf_n_vl( 3, nkmax, nshell))
      if( .not. allocated( wf_n_vc))            allocate( wf_n_vc( 3, nkmax, nshell))
      if( .not. allocated( wf_n_wgt))           allocate( wf_n_wgt( nshell))
      if( .not. allocated( wf_n_usedshells))    allocate( wf_n_usedshells( nshell))
      if( .not. allocated( wf_n_ns2n))          allocate( wf_n_ns2n( nkmax, nshell))

      ! sort everything by distance
      wf_n_n = 0
      wf_n_dist = 0.d0
      wf_n_vl = 0.d0
      wf_n_vc = 0.d0
      dist = minval( dt)
      do i = 1, nshell
        j = minloc( abs( dt( 1:nshell) - dist), 1)
        wf_n_dist( i) = dt( j)
        wf_n_n( i) = mt( j)
        wf_n_vl( :, :, i) = nvlt( :, :, j)
        wf_n_vc( :, :, i) = nvct( :, :, j)
        dt( j) = 1.d13 
        dist = minval( dt)
      end do
      
      ! find k-point indices of neighbors
      do iknr = 1, wf_kset%nkpt
        do j = 1, nshell
          do i = 1, wf_n_n( j)
            vl = wf_kset%vkl( :, iknr) + wf_n_vl( :, i, j)
            call r3frac( input%structure%epslat, vl, vi) 
            do iy = 1, wf_kset%nkpt
              dt( iy) = norm2( vl - wf_kset%vkl( :, iy))
            end do
            iz = minloc( dt, 1)
            if( norm2( vl - wf_kset%vkl( :, iz)) .lt. input%structure%epslat) then
              wf_n_ik( i, j, iknr) = iz
              !write(*,'(3(3F13.6,3x))') wf_kset%vkl( :, iknr), wf_n_vl( :, i, j), wf_kset%vkl( :, wf_n_ik( i, j, iknr))
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
        wf_n_wgt( 1) = 0.5d0/norm2( wf_n_vc( :, 1, 1))**2
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

      do i = 1, wf_n_nshells
        ix = wf_n_usedshells( i)
        do j = 1, wf_n_n( ix)
          write(*,*) i, ix, j, wf_n_ns2n( j, ix)
        end do
      end do

      return
      !EOC
    end subroutine wannier_geometry
    !EOP
    
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
        if( wf_projused( i)) then
          write( wf_info, '("[X]")', advance='no')
        else
          write( wf_info, '("[ ]")', advance='no')
        end if
        write( wf_info, *)
      end do

      write( wf_info, '(80("-"))')
      write( wf_info, '(36x,"local-orbitals used in total:",4x,I4,"/",I4)') sum( wf_projused), wf_nprojtot
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
          write( wf_info, '(12x,I2,3(F22.10))') j, wf_n_vc( :, j, wf_n_usedshells( i))
          v(:,1) = wf_n_vc( :, j, wf_n_usedshells( i))
          m = m + wf_n_wgt( wf_n_usedshells( i))*matmul( v, transpose( v)) 
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
      write( wf_info, '(" #Wannier functions:",T30,13x,I3)') wf_nst
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
      integer :: i

      call printbox( wf_info, '*', "Wannier functions")
      write( wf_info, *)
      write( wf_info, '(6x,"band",20x,"localization center (cartesian)",14x,"Omega")')
      write( wf_info, '(80("-"))')

      do i = wf_fst, wf_lst
        write( wf_info, '(6x,I4,3x,3F16.10,3x,F16.10)') i, wf_centers( :, i), wf_omega( i)
      end do

      write( wf_info, '(80("-"))')
      write( wf_info, '(55x,"total:",3x,F16.10)') sum( wf_omega)
      write( wf_info, '(53x,"average:",3x,F16.10)') sum( wf_omega)/wf_nst

      write( wf_info, *)
      call flushifc( wf_info)
      !close( wf_info)
    end subroutine wannier_writeinfo_results

    subroutine wannier_destroy
      if( allocated( wf_projst)) deallocate( wf_projst)
      if( allocated( wf_projused)) deallocate( wf_projused)
      if( allocated( wf_win_ni)) deallocate( wf_win_ni)
      if( allocated( wf_win_no)) deallocate( wf_win_no)
      if( allocated( wf_win_idxi)) deallocate( wf_win_idxi)
      if( allocated( wf_win_idxo)) deallocate( wf_win_idxo)
      if( allocated( wf_transform)) deallocate( wf_transform)
      if( allocated( wf_projection)) deallocate( wf_projection)
      if( allocated( wf_opf)) deallocate( wf_opf)
      if( allocated( wf_emat)) deallocate( wf_emat)
      if( allocated( wf_subspace)) deallocate( wf_subspace)
      if( allocated( wf_centers)) deallocate( wf_centers)
      if( allocated( wf_omega)) deallocate( wf_omega)
      if( allocated( wf_n_dist)) deallocate( wf_n_dist)
      if( allocated( wf_n_n)) deallocate( wf_n_n)
      if( allocated( wf_n_ik)) deallocate( wf_n_ik)
      if( allocated( wf_n_vl)) deallocate( wf_n_vl)
      if( allocated( wf_n_vc)) deallocate( wf_n_vc)
      if( allocated( wf_n_wgt)) deallocate( wf_n_wgt)
      wf_initialized = .false.
    end subroutine wannier_destroy
end module mod_wannier
