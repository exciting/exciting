module mod_wannier
  use mod_wannier_variables
  use mod_wannier_projection
  use mod_wannier_maxloc
  use mod_wannier_opf
  use mod_wannier_disentangle
  use mod_wannier_filehandling
  use mod_wannier_helper
  use mod_optkgrid
  use mod_kpointset
  use xlapack, only: svd_divide_conquer

  use mod_pwmat

  implicit none

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

      if( wf_initialized) call wannier_destroy
      call init0
      call init1

      !********************************************************************
      ! open INFO file
      !********************************************************************
      write( wf_filename, '("WANNIER")')
      call getunit( wf_info)
      if( mpiglobal%rank .eq. 0) then
        if( input%properties%wannier%do .ne. "fromfile") then
          open( wf_info, file=trim( wf_filename)//"_INFO"//trim(filext), action='write', form='formatted')
        end if
      end if
      call timesec( wf_t0)

      !********************************************************************
      ! load used k-grid
      !********************************************************************
      call wfhelp_setkpts

      !********************************************************************
      ! check input parameters
      !********************************************************************
      call wannier_readinput

      !********************************************************************
      ! find geometry
      !********************************************************************
      call wannier_geometry( wf_kset, wf_n_ntot, wf_n_vl, wf_n_vc, wf_n_ik, wf_n_ik2, wf_n_wgt, wf_n_dist, &
             nbzshell=input%properties%wannier%nbzshell, &
             nomult=.true.)

      if( input%properties%wannier%do .ne. "fromfile") call wffile_writeinfo_geometry
      
      !********************************************************************
      ! set up R-vector set
      !********************************************************************
      call wannier_rvectors( input%structure%crystal%basevect, wf_kset, wf_nrpt, wf_rvec, wf_rmul, wf_pkr)

      !********************************************************************
      ! generate radial functions and find phasefactors
      !********************************************************************
      call wfhelp_genradfun
      call wfhelp_fixphases

      !********************************************************************
      ! calculate plane-wave matrix elements
      !********************************************************************
      select case( input%properties%wannier%do)
        ! dont do
        case( "skip")
        case( "fromfile")
        ! do
        case default
          if( mpiglobal%rank .eq. 0) then
            call printbox( wf_info, '*', "Preparation")
            write( wf_info, *)
          end if
          call wannier_emat
      end select

      !********************************************************************
      ! initialize transformation matrices
      !********************************************************************
      call wannier_init_transform
        
      !********************************************************************
      ! build projection functions and overlap matrices
      !********************************************************************
      select case( input%properties%wannier%do)
        ! dont do
        case( "skip")
        case( "fromfile")
        case( "maxfromfile")
        ! do
        case default
          call wfpro_projection
      end select

      wf_initialized = .true.
      !call wannier_writesetup

      return
    end subroutine wannier_init
    !EOC

    subroutine wannier_init_transform
      integer :: igroup, ik, i, n, fst, lst
      real(8) :: fac, win(2)

      integer, allocatable :: idxo(:)
      real(8), allocatable :: eval(:,:), score(:,:), scoresum(:)

      fac = 4.d0

      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_fst:wf_lst, wf_nwf, wf_kset%nkpt))
      wf_transform = zzero

      if( allocated( wf_phases)) deallocate( wf_phases)
      allocate( wf_phases( wf_nwf, wf_kset%nkpt))
      wf_phases = 0.d0

      call wfhelp_geteval( eval, fst, lst)
      if( wf_fermizero) then
        eval = eval - wf_efermi
      end if

      do igroup = 1, wf_ngroups
        if( wf_groups( igroup)%method .eq. 'disSMV' .or. wf_groups( igroup)%method .eq. 'disFull') then
          allocate( idxo( wf_groups( igroup)%nwf))
          allocate( scoresum( wf_groups( igroup)%fst:wf_groups( igroup)%lst))
          allocate( score( wf_groups( igroup)%fst:wf_groups( igroup)%lst, wf_kset%nkpt))
          score = 0.d0
          idxo = 0
          win = wf_groups( igroup)%win_i
          if( sum( abs( win)) .lt. 1.d-23) win = wf_groups( igroup)%win_o
          do ik = 1, wf_kset%nkpt
            do i = 1, wf_groups( igroup)%win_ni( ik)
              n = wf_groups( igroup)%win_ii( i, ik)
              score( n, ik) = score( n, ik) + 1.d0
            end do
            do i = 1, wf_groups( igroup)%win_no( ik)
              n = wf_groups( igroup)%win_io( i, ik)
              score( n, ik) = score( n, ik) + &
                1.d0/( 1.d0 + exp( -fac*( eval( n, ik) - win(1))/( win(2) - win(1)))) * &
                1.d0/( 1.d0 + exp(  fac*( eval( n, ik) - win(2))/( win(2) - win(1))))
            end do
          end do
          scoresum = sum( score, 2)/wf_kset%nkpt
          i = 0
          do while( i .lt. wf_groups( igroup)%nwf)
            n = wf_groups( igroup)%fst + maxloc( scoresum, 1) - 1
            i = i + 1
            idxo( i) = n
            scoresum( n) = 0.d0
          end do
          do i = 1, wf_groups( igroup)%nwf
            wf_transform( idxo( i), wf_groups( igroup)%fwf+i-1, :) = zone
          end do
          deallocate( idxo, scoresum, score)
        else
          do i = 1, wf_groups( igroup)%nwf
            wf_transform( wf_groups( igroup)%fst+i-1, wf_groups( igroup)%fwf+i-1, :) = zone
          end do
        end if
      end do
      deallocate( eval)

      return
    end subroutine wannier_init_transform
        
    !BOP
    ! !ROUTINE: wannier_emat
    ! !INTERFACE:
    !
    subroutine wannier_emat
      ! !USES:
      use mod_symmat
      use mod_symmetry
      ! !DESCRIPTION:
      !   Generates plane-wave matrix elements for neighboring k-points
      !   needed in order to calculate well/maximally localized Wannier
      !   functions.
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      ! local variables
      integer :: k1, k2, iknr, idxn, cntk
      real(8) :: t0, t1
      logical :: success

      ! allocatable arrays
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)

      write( wf_info, '(" calculate plane-wave matrix-elements...")')
      call timesec( t0)

      ! check for existing file
      success = .true.
      inquire( file=trim( wf_filename)//"_EMAT"//trim( filext), exist=success)
      if( success) then
        call wffile_reademat( success)
        if( mpiglobal%rank .eq. 0) then
          if( success) then
            call timesec( t1)
            write( wf_info, '(5x,"plane-wave matrix-elements read from file")')
            write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
            write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
            write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
          else
            write(*, '( "Recalculate plane-wave matrix-elements.")')
          end if
        end if
      end if

      if( .not. success) then
        if( allocated( wf_m0)) deallocate( wf_m0)
        allocate( wf_m0( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt, wf_n_ntot))

        wf_m0(:,:,:,:) = zzero

        call pwmat_init( input%groundstate%lmaxapw, 8, wf_kset, wf_fst, wf_lst, wf_fst, wf_lst, &
               fft=.false.)

        allocate( evecfv1( nmatmax_ptr, nstfv, nspinor))
        allocate( evecfv2( nmatmax_ptr, nstfv, nspinor))

        k1 = firstofset( mpiglobal%rank, wf_kset%nkpt)
        k2 = lastofset( mpiglobal%rank, wf_kset%nkpt)

        do iknr = k1, k2
          call wfhelp_getevec( iknr, evecfv1)
          call pwmat_prepare( iknr, evecfv1( :, :, 1))
        end do
        call barrier
          
        do idxn = 1, wf_n_ntot
          call pwmat_init_qg( wf_n_vl( :, idxn), (/0, 0, 0/), 1)
          cntk = 0
          do iknr = k1, k2
            ! read eigenvectors
            call wfhelp_getevec( iknr, evecfv1)
            call wfhelp_getevec( wf_n_ik( idxn, iknr), evecfv2)
            ! generate plane-wave matrix elements
            call pwmat_genpwmat( iknr, &
                   evecfv1( :, wf_fst, 1), &
                   evecfv2( :, wf_fst, 1), &
                   wf_m0( :, :, iknr, idxn))
            cntk = cntk + 1
            if( mpiglobal%rank .eq. 0) then
              write(*,'(1a1,"Calculating plane-wave matrix elements (neighbor ",i2.2," of ",i2.2,"): ",f10.3,"%",$)') achar(13), idxn, wf_n_ntot, &
                  100.d0*(dble( idxn - 1)/wf_n_ntot + dble( cntk)/(k2 - k1 + 1)/wf_n_ntot)
            end if
          end do
          call mpi_allgatherv_ifc( set=wf_kset%nkpt, rlen=wf_nst*wf_nst, zbuf=wf_m0( :, :, :, idxn))
          call barrier
        end do
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
        end if
        deallocate( evecfv1, evecfv2)

        call pwmat_destroy
        call wffile_writeemat
        call timesec( t1)
        if( mpiglobal%rank .eq. 0) then
          write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
          write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
          write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
        end if
      end if

      if( mpiglobal%rank .eq. 0) then
        write( wf_info, *)
        call flushifc( wf_info)
      end if
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
      integer :: iknr
      real(8) :: t0, t1
      real(8), allocatable :: sval(:)
      complex(8), allocatable :: lsvec(:,:), rsvec(:,:)

      if( mpiglobal%rank .eq. 0) then
        write( wf_info, '(" perform simple projection step...")')
      end if
      call timesec( t0)

      !********************************************************************
      ! build transformation matrices
      !********************************************************************

      allocate( sval( wf_groups( wf_group)%nwf), &
                lsvec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst), &
                rsvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, sval, lsvec, rsvec)
!$OMP DO
#endif
      do iknr = 1, wf_kset%nkpt
        call svd_divide_conquer( wf_groups( wf_group)%projection(:,:,iknr), sval, lsvec, rsvec)
        ! for numerical stability
        !do i = 1, wf_nwf
        !  if( sval( i) .lt. 1.d-12) lsvec( :, i) = zzero
        !end do
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               lsvec, wf_groups( wf_group)%nwf, &
               rsvec, wf_groups( wf_group)%nwf, zzero, &
               wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, iknr), wf_groups( wf_group)%nst)
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      call wfomega_gen
      deallocate( sval, lsvec, rsvec)
      call timesec( t1)

      if( mpiglobal%rank .eq. 0) then
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
        write( wf_info, *)
        call flushifc( wf_info)
      end if

      return
    end subroutine wannier_gen_pro
    !EOC
    
    subroutine wannier_gen_fromfile
      logical :: success

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call wffile_readtransform( success)
      if( .not. success) then
        if( mpiglobal%rank .eq. 0) then
          write(*, '(" Failed to read Wannier functions from file. Recalculate them.")')
        end if
        input%properties%wannier%do = "fromscratch"
        call wannier_init
        do wf_group = 1, wf_ngroups
          if( mpiglobal%rank .eq. 0) call wfopf_gen
          if( mpiglobal%rank .eq. 0) call wfmax_gen
        end do
      end if
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      return
    end subroutine wannier_gen_fromfile

    subroutine wannier_geometry( kset, nn, nvl, nvc, nik, nik2, nwgt, ndist, nbzshell, nomult)
      ! !USES:
      use sorting, only: sort_index_1d
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
      type( k_set), intent( in)          :: kset
      integer, intent( out)              :: nn 
      real(8), allocatable, intent( out) :: nvl(:,:), nvc(:,:), nwgt(:), ndist(:)
      integer, allocatable, intent( out) :: nik(:,:), nik2(:,:)
      integer, optional, intent( in)     :: nbzshell(3)
      logical, optional, intent( in)     :: nomult

      integer :: ix, iy, iz, i, j, nposvec, iknr, d, nbzs(3), norder, o
      integer :: vi(3)
      real(8) :: vl(3), vc(3), wgord(3,3)
      real(8) :: rot(3,3), coeff(6,12), right(6), sval(6), coeffinv(12,6)
      real(8) :: mat1(3,3), mat2(3,3), mat3(3,3), umat(3,3), u(3), n(3), nvec(3,120), c, s
      logical :: stopshell, nomults

      real(8), allocatable :: dt(:)
      integer, allocatable :: idx(:), posvec(:,:)
      real(8), allocatable :: tmp_n_vl(:,:), tmp_n_vc(:,:), tmp_n_wgt(:)

      nbzs = 1
      if( present( nbzshell)) nbzs = nbzshell
      nomults = .false.
      if( present( nomult)) nomults = nomult
      
      ! find possible distances (shells)
      nposvec = (nbzs(3)*kset%ngridk(3)-1)*(2*nbzs(2)*kset%ngridk(2)-1)*(2*nbzs(1)*kset%ngridk(1)-1) + &
                (nbzs(2)*kset%ngridk(2)-1)*(2*nbzs(1)*kset%ngridk(1)-1) + &
                 nbzs(1)*kset%ngridk(1)-1
      allocate( dt( nposvec))
      allocate( posvec( 3, nposvec))

      i = 0
      posvec(:,:) = 0
      do iz = nbzs(3)*kset%ngridk(3)-1, -nbzs(3)*kset%ngridk(3)+1, -1
        do iy = nbzs(2)*kset%ngridk(2)-1, -nbzs(2)*kset%ngridk(2)+1, -1
          do ix = nbzs(1)*kset%ngridk(1)-1, -nbzs(1)*kset%ngridk(1)+1, -1
            vi = (/ix, iy, iz/)
            if( (abs(ix)+abs(iy)+abs(iz)) .ne. 0) then
              stopshell = .true.
              d = maxloc( abs( vi), 1)
              do j = 1, i
                if( nomults .and. (posvec( d, j) .ne. 0)) then
                  s = dble( vi(d))/dble( posvec( d, j))
                  if( norm2( s*dble( posvec( :, j)) - dble( vi)) .lt. input%structure%epslat) then
                    stopshell = .false.
                    if( norm2( dble( vi)) .lt. norm2( dble( posvec( :, j)))) posvec( :, j) = vi
                    exit
                  end if
                else
                  if( norm2( dble( posvec( :, j) + vi)) .lt. input%structure%epslat) then
                    stopshell = .false.
                    exit
                  end if
                end if
              end do
              if( stopshell) then
                i = i + 1
                posvec(:,i) = vi
              end if
            end if
          end do
        end do
      end do
      nposvec = i

      dt = 0.d0
      do ix = 1, nposvec
        vl = dble( posvec(:,ix))/kset%ngridk
        call r3mv( kset%bvec, vl, vc)
        dt( ix) = norm2( vc)
      end do

      ! sort neighbors by distance
      allocate( idx( nposvec))
      idx( 1:nposvec) =  sort_index_1d( nposvec, dt( 1:nposvec))
      posvec( :, 1:nposvec) = posvec( :, idx( 1:nposvec))

      ! find all possible neighbors
      allocate( tmp_n_vl( 3, nposvec))
      allocate( tmp_n_vc( 3, nposvec))

      tmp_n_vl = 0.d0
      tmp_n_vc = 0.d0
      do ix = 1, nposvec
        tmp_n_vl( :, ix) = dble( posvec(:,ix))/kset%ngridk
        call r3mv( kset%bvec, tmp_n_vl( :, ix), tmp_n_vc( :, ix))
      end do

      ! find number of shells needed for gradient calculation and geometric weights
      allocate( tmp_n_wgt( nposvec))

      j = 1
      stopshell = .false.
      coeff = 0.d0
      tmp_n_wgt = 0.d0

      if( sum( kset%ngridk) .eq. 1) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write(*,'("Error (wannier_geometry): Wannier Functions for single k-point only not yet implemented.")')
        end if
        stop
      else if( minval( kset%ngridk) .eq. 1) then
        if( sum( kset%ngridk) .eq. maxval( kset%ngridk) + 2) then
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
      rot = 0.d0
      rot(1,1) = 1.d0
      rot(2,2) = 1.d0
      rot(3,3) = 1.d0

      if( all( nbzs*kset%ngridk .gt. 1)) d = 6
      idx = 0
      nn = 0
      nvec = 0.d0

      if( d .gt. 1) then
        if( d .eq. 3) then
          ix = minloc( kset%ngridk, 1)
          iy = mod( ix+1, 3)
          iz = mod( ix+2, 3)
          if( iy .eq. 0) iy = 3
          if( iz .eq. 0) iz = 3

          ! rotate into xy-plane
          mat3 = 0.d0
          mat3(1,1) = 1.d0
          mat3(2,2) = 1.d0
          mat3(3,3) = 1.d0
          call r3cross( kset%bvec( :, iy), kset%bvec( :, iz), n)
          n = n/norm2( n)
          if( abs( abs( n(3)) - 1.d0) .gt. input%structure%epslat) then
            c = n(3)
            s = dsqrt( 1 - c*c)
            call r3cross( n, (/0.d0, 0.d0, 1.d0/), u)
            u = u/norm2( u)
            umat = 0.d0
            umat(1,2) = -u(3)
            umat(2,1) =  u(3)
            umat(1,3) =  u(2)
            umat(3,1) = -u(2)
            umat(2,3) = -u(1)
            umat(3,2) =  u(1)
            call r3mm( umat, umat, rot)
            rot = (1-c)*rot + s*umat + mat3
          end if
          mat3(3,3) = 0.d0
          call r3mm( mat3, rot, umat)
          call r3mm( transpose( rot), umat, mat3)
          call r3mm( mat1, mat3, umat)
          call r3mm( umat, mat2, mat3)
          right(:) = (/mat3(1,1), mat3(2,2), mat3(1,2), 0.d0, 0.d0, 0.d0/)

        else
          right(:) = (/mat3(1,1), mat3(2,2), mat3(3,3), mat3(1,2), mat3(1,3), mat3(2,3)/)
        end if
        do while( .not. stopshell)
          ! check if new shell is linear dependend
          idx( j) = 1
          !write(*,*) 'NEIGHBOR', j
          !write(*,'(3f13.6)') tmp_n_vl( :, j)
          !if( nn .gt. 0) then
          !  call rpinv( nvec( :, 1:nn), nveci( 1:nn, :))
          !  idx( j) = 0
          !  if( norm2( tmp_n_vl( :, j) - matmul( matmul( nvec( :, 1:nn), nveci( 1:nn, :)), tmp_n_vl( :, j))) &
          !  .gt. input%structure%epslat) idx( j) = 1
          !end if              
          !write(*,*) idx( j)

          if( idx( j) .ne. 0) then
            nn = nn + 1
            !write(*,*) 'nneighbors', nn
            if( d .eq. 3) then
              coeff( 1, nn) = tmp_n_vl( iy, j)**2
              coeff( 2, nn) = tmp_n_vl( iz, j)**2
              coeff( 3, nn) = tmp_n_vl( iy, j)*tmp_n_vl( iz, j)
              nvec( :, nn) = tmp_n_vl( :, j)
            else
              coeff( 1, nn) = tmp_n_vl( 1, j)**2
              coeff( 2, nn) = tmp_n_vl( 2, j)**2
              coeff( 3, nn) = tmp_n_vl( 3, j)**2
              coeff( 4, nn) = tmp_n_vl( 1, j)*tmp_n_vl( 2, j)
              coeff( 5, nn) = tmp_n_vl( 1, j)*tmp_n_vl( 3, j)
              coeff( 6, nn) = tmp_n_vl( 2, j)*tmp_n_vl( 3, j)
              nvec( :, nn) = tmp_n_vl( :, j)
            end if
            call rpinv( coeff( 1:d, 1:nn), coeffinv( 1:nn, 1:d))
            sval( 1:d) = matmul( matmul( coeff( 1:d, 1:nn), coeffinv( 1:nn, 1:d)), right( 1:d))
            tmp_n_wgt( 1:nn) = matmul( coeffinv( 1:nn, 1:d), right( 1:d))
            if( tmp_n_wgt( nn) .le. 1.d-6) then
              idx( j) = 0
              nn = nn -1
            end if
            if( norm2( sval( 1:d) - right( 1:d)) .lt. input%structure%epslat) then
              stopshell = .true.
            end if
          else
            write(*,'("neighbor ",i3," is linear dependend")') j
          end if
          j = j+1
        end do
      else
        tmp_n_wgt( 1) = 1.0d0/norm2( tmp_n_vc(:,1))**2
        nn = 1
      end if
      do j = 1, nn
        if( abs( tmp_n_wgt( j)) .lt. 1.d-6) then
          nn = nn - 1
          idx( j) = -1
        end if
      end do

      tmp_n_wgt = 0.5d0*tmp_n_wgt

      norder = 1
      wgord = 0.d0
      ! first oder gradient approximation
      wgord( 1, 1) =  1.d0
      ! second order gradient approximation
      wgord( 1, 2) =  4.d0/3.d0
      wgord( 2, 2) = -1.d0/12.d0

      ! find k-point indices of neighbors and copy results to global arrays
      if( allocated( nvl))            deallocate( nvl)
      allocate( nvl( 3, norder*nn))
      if( allocated( nvc))            deallocate( nvc)
      allocate( nvc( 3, norder*nn))
      if( allocated( nik))            deallocate( nik)
      allocate( nik( norder*nn, kset%nkptnr))
      if( allocated( nik2))           deallocate( nik2)
      allocate( nik2( norder*nn, kset%nkptnr))
      if( allocated( nwgt))           deallocate( nwgt)
      allocate( nwgt( norder*nn))
      if( allocated( ndist))          deallocate( ndist)
      allocate( ndist( norder*nn))

      d = 0
      i = 0
      do j = 1, nposvec
        if( d .eq. nn) exit
        if( idx( j) .ne. 0) then
          i = i + 1
          if( idx( j) .ne. -1) then
            d = d + 1
            do o = 1, norder
              nvl( :, (o-1)*nn+d) = dble( o)*tmp_n_vl( :, j)
              nvc( :, (o-1)*nn+d) = dble( o)*tmp_n_vc( :, j)
              nwgt( (o-1)*nn+d) = wgord( o, norder)*tmp_n_wgt( i)
              ndist( (o-1)*nn+d) = norm2( nvc( :, (o-1)*nn+d))
              do iknr = 1, kset%nkpt
                vl = kset%vkl( :, iknr) + nvl( :, (o-1)*nn+d)
                call r3frac( input%structure%epslat, vl, vi) 
                call findkptinset( vl, kset, ix, iy)
                if( ix .gt. 0) then
                  nik( (o-1)*nn+d, iknr) = iy
                else
                  if( mpiglobal%rank .eq. 0) then
                    write(*,*)
                    write( *, '("Error (wannier_geometry): wrong neighboring vector.")')
                    write( *, '(x,3F13.6)') nvl( :, (o-1)*nn+d)
                  end if
                  stop
                end if

                vl = kset%vkl( :, iknr) - nvl( :, (o-1)*nn+d)
                call r3frac( input%structure%epslat, vl, vi) 
                call findkptinset( vl, kset, ix, iy)
                if( ix .gt. 0) then
                  nik2( (o-1)*nn+d, iknr) = iy
                else
                  if( mpiglobal%rank .eq. 0) then
                    write(*,*)
                    write( *, '("Error (wannier_geometry): wrong neighboring vector.")')
                    write( *, '(x,3F13.6)') nvl( :, (o-1)*nn+d)
                  end if
                  stop
                end if
              end do
            end do
          end if
        end if
      end do

      nn = nn*norder
      
      deallocate( dt, idx, posvec)
      deallocate( tmp_n_vl, tmp_n_vc, tmp_n_wgt)

      return
      !EOC
    end subroutine wannier_geometry

    subroutine wannier_rvectors( latvec, kset, nrpt, rvec, rmul, pkr)
      use sorting, only: sort_index_1d
      real(8), intent( in)                         :: latvec(3,3)
      type( k_set), intent( in)                    :: kset

      integer, intent( out)                        :: nrpt
      integer, allocatable, intent( out)           :: rvec(:,:)
      integer, allocatable, intent( out)           :: rmul(:)
      complex(8), allocatable, intent( out)        :: pkr(:,:)

      integer :: i, j, k, is, js, ks, ir, ik, cnt
      integer :: tmpvec( 3, 4*kset%nkpt), tmpmul( 4*kset%nkpt), idx( 4*kset%nkpt)
      real(8) :: v1(3), v2(3), vs(3), vc(3), d, dist(125), length( 4*kset%nkpt)
      real(8) :: ilatvec(3,3)

      call r3minv( latvec, ilatvec)

      tmpvec = 0
      tmpmul = 0
      idx = 0
      length = 0.d0
      ! map R vectors from canonical SZ into its Wigner-Seitz cell
      ! and count equivalent vectors
      nrpt = 0
      do i = -kset%ngridk(1), kset%ngridk(1)
        v1 = i*latvec( :, 1)
        do j = -kset%ngridk(2), kset%ngridk(2)
          v2 = v1 + j*latvec( :, 2)
          do k = -kset%ngridk(3), kset%ngridk(3)
            vc = v2 + k*latvec( :, 3)
            cnt = 0
            do is = -2, 2               
              do js = -2, 2             
                do ks = -2, 2           
                  cnt = cnt + 1
                  call r3mv( latvec, dble( (/is, js, ks/)*kset%ngridk), vs)
                  dist( cnt) = norm2( vc - vs)
                end do
              end do
            end do
            d = minval( dist)
            if( abs( dist( 63) - d) .lt. input%structure%epslat) then
              nrpt = nrpt + 1
              do ir = 1, 125
                if( abs( dist( ir) - d) .lt. input%structure%epslat) tmpmul( nrpt) = tmpmul( nrpt) + 1
              end do
              tmpvec( :, nrpt) = (/i, j, k/)
              length( nrpt) = norm2( vc)
            end if
          end do
        end do
      end do
      ! sort R-vectors with increasing length
      idx( 1:nrpt) = sort_index_1d( nrpt, length( 1:nrpt))

      if( allocated( rvec)) deallocate( rvec)
      allocate( rvec( 3, nrpt))
      if( allocated( rmul)) deallocate( rmul)
      allocate( rmul( nrpt))

      do ir = 1, nrpt
        rvec( :, ir) = tmpvec( :, idx( ir))
        rmul( ir) = tmpmul( idx( ir))
      end do

      ! precompute phases for Fourier sums
      if( allocated( pkr)) deallocate( pkr)
      allocate( pkr( kset%nkpt, nrpt))

      do ir = 1, nrpt
        do ik = 1, kset%nkpt
          d = -twopi*dot_product( kset%vkl( :, ik), dble( rvec( :, ir)))
          pkr( ik, ir) = cmplx( cos( d), sin( d), 8)
        end do
      end do

      return
    end subroutine wannier_rvectors

    subroutine wannier_mindist( latvec, kset, nrpt, rvec, nst, centers, wdistvec, wdistmul)
      real(8), intent( in)               :: latvec(3,3)
      type( k_set), intent( in)          :: kset
      integer, intent( in)               :: nrpt
      integer, intent( in)               :: rvec( 3, nrpt)
      integer, intent( in)               :: nst
      real(8), intent( in)               :: centers( 3, nst)

      integer, allocatable, intent( out) :: wdistvec(:,:,:,:,:)
      integer, allocatable, intent( out) :: wdistmul(:,:,:)

      integer :: ir, i, j, is, js, ks
      real(8) :: ilatvec(3,3), vc(3), vs(3), vl(3), v1(3), d

      call r3minv( latvec, ilatvec)

      if( allocated( wdistvec)) deallocate( wdistvec)
      allocate( wdistvec( 3, 8, nst, 0:nst, nrpt))
      if( allocated( wdistmul)) deallocate( wdistmul)
      allocate( wdistmul( nst, 0:nst, nrpt))

      wdistvec = 0
      wdistmul = 0
#ifdef USEOMP
!$omp parallel default( shared) private( ir, j, i, vc, v1, d, is, js, ks, vs, vl)
!$omp do collapse( 3)
#endif
      do ir = 1, nrpt
        do j = 0, nst
          do i = 1, nst
            call r3mv( latvec, dble( rvec( :, ir)), vc)
            if( j .eq. 0) then
              vc = vc + centers( :, i)
            else
              vc = vc + centers( :, j) - centers( :, i)
            end if
            v1 = vc
            d = norm2( v1)
            ! map to super-WS cell
            do is = -2, 2    
              do js = -2, 2             
                do ks = -2, 2           
                  call r3mv( latvec, dble( (/is, js, ks/)*kset%ngridk), vs)
                  if( norm2( vc - vs) .lt. d) then
                    v1 = vc - vs
                    d = norm2( v1)
                  end if
                end do
              end do
            end do
            vc = v1
            ! find equivalent vectors
            do is = -2, 2    
              do js = -2, 2             
                do ks = -2, 2           
                  call r3mv( latvec, dble( (/is, js, ks/)*kset%ngridk), vs)
                  if( abs( norm2( vc - vs) - d) .lt. input%structure%epslat) then
                    wdistmul( i, j, ir) = wdistmul( i, j, ir) + 1
                    if( j .eq. 0) then
                      v1 = vs
                    else
                      v1 = vc - vs + centers( :, i) - centers( :, j)
                    end if
                    call r3mv( ilatvec, v1, vl)
                    wdistvec( :, wdistmul( i, j, ir), i, j, ir) = nint( vl)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      do ir = 1, nrpt
        do i = 1, nst
          wdistmul( i, i, ir) = 1
          wdistvec( :, 1, i, i, ir) = rvec( :, ir)
        end do
      end do

      return
    end subroutine wannier_mindist
    
    subroutine wannier_readinput
      integer :: igroup, ik, fst, lst, ist, un
      integer :: fst_, lst_, nst_, nwf_, nkpt_
      logical :: disentangle_, success
      real(8) :: e
      real(8), allocatable :: eval(:,:), occ(:,:)

      wf_fermizero = input%properties%wannier%fermizero
      wf_ngroups = size( input%properties%wannier%grouparray, 1)
      allocate( wf_groups( wf_ngroups))

      ! initialize Fermi energy
      call wfhelp_getefermi( wf_efermi)
      call wfhelp_geteval( eval, fst, lst)
      allocate( occ( fst:lst, wf_kset%nkpt))
      call wfhelp_occupy( wf_kset, eval, fst, lst, wf_efermi, occ)
      deallocate( occ)

      ! read input from input.xml
      if( input%properties%wannier%do .ne. "skip") then
        do igroup = 1, wf_ngroups
          wf_groups( igroup)%method = input%properties%wannier%grouparray( igroup)%group%method
          wf_groups( igroup)%neighcells = input%properties%wannier%grouparray( igroup)%group%neighcells
          wf_groups( igroup)%win_o(1) = minval( input%properties%wannier%grouparray( igroup)%group%outerwindow)
          wf_groups( igroup)%win_o(2) = maxval( input%properties%wannier%grouparray( igroup)%group%outerwindow)
          wf_groups( igroup)%nproj = 0

          ! set energy windows
          if( norm2( wf_groups( igroup)%win_o) .gt. 1.d-20) then
            wf_groups( igroup)%win_i(1) = minval( input%properties%wannier%grouparray( igroup)%group%innerwindow)
            wf_groups( igroup)%win_i(2) = maxval( input%properties%wannier%grouparray( igroup)%group%innerwindow)

            if( norm2( wf_groups( igroup)%win_i) .gt. 1.d-20) then
              if( (wf_groups( igroup)%win_i(1) .lt. wf_groups( igroup)%win_o(1)) .or. &
                  (wf_groups( igroup)%win_i(2) .gt. wf_groups( igroup)%win_o(2))) then
                if( mpiglobal%rank .eq. 0) then
                  write(*,*)
                  write( *, '("Error (wannier_readinput): The inner window must be fully contained within the outer window for group ",I2,".")') igroup
                end if
                stop
              end if
            end if

            allocate( wf_groups( igroup)%win_ii( lst-fst+1, wf_kset%nkpt), wf_groups( igroup)%win_io( lst-fst+1, wf_kset%nkpt))
            allocate( wf_groups( igroup)%win_ni( wf_kset%nkpt), wf_groups( igroup)%win_no( wf_kset%nkpt))
            ! assign states to windows
            wf_groups( igroup)%win_ii = 0
            wf_groups( igroup)%win_io = 0
            wf_groups( igroup)%win_ni = 0
            wf_groups( igroup)%win_no = 0

            wf_groups( igroup)%fst = 100000000
            wf_groups( igroup)%lst = 0
            do ik = 1, wf_kset%nkpt
              do ist = fst, lst
                e = eval( ist, ik)
                ! shift by Fermi energy
                if( wf_fermizero) e = e - wf_efermi
                ! state in inner window
                if( (e .ge. wf_groups( igroup)%win_i(1)) .and. (e .le. wf_groups( igroup)%win_i(2))) then
                  wf_groups( igroup)%win_ni( ik) = wf_groups( igroup)%win_ni( ik) + 1
                  wf_groups( igroup)%win_ii( wf_groups( igroup)%win_ni( ik), ik) = ist
                ! state in outer (but not in inner) window
                else if( (e .ge. wf_groups( igroup)%win_o(1)) .and. (e .le. wf_groups( igroup)%win_o(2))) then
                  wf_groups( igroup)%win_no( ik) = wf_groups( igroup)%win_no( ik) + 1
                  wf_groups( igroup)%win_io( wf_groups( igroup)%win_no( ik), ik) = ist
                end if
                !write(*,'(2i,2f13.6,2i)') ik, ist, eval( ist, ik), e, wf_groups( igroup)%win_ni( ik), wf_groups( igroup)%win_no( ik)
              end do
              ! set band ranges
              wf_groups( igroup)%fst = min( wf_groups( igroup)%fst, minval( wf_groups( igroup)%win_ii( 1:wf_groups( igroup)%win_ni( ik), ik)))
              wf_groups( igroup)%fst = min( wf_groups( igroup)%fst, minval( wf_groups( igroup)%win_io( 1:wf_groups( igroup)%win_no( ik), ik)))
              wf_groups( igroup)%lst = max( wf_groups( igroup)%lst, maxval( wf_groups( igroup)%win_ii( 1:wf_groups( igroup)%win_ni( ik), ik)))
              wf_groups( igroup)%lst = max( wf_groups( igroup)%lst, maxval( wf_groups( igroup)%win_io( 1:wf_groups( igroup)%win_no( ik), ik)))
            end do
            ! set number of Wannier functions
            wf_groups( igroup)%nwf = input%properties%wannier%grouparray( igroup)%group%nwf
            ! if not specified automatically select number
            if( wf_groups( igroup)%nwf .lt. 1) then
              wf_groups( igroup)%nwf = nint( 0.5d0*(maxval( wf_groups( igroup)%win_ni) + minval( wf_groups( igroup)%win_ni + wf_groups( igroup)%win_no)))
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Warning (wannier_readinput): Number of Wannier functions (nwf) not set for group ",I2,". I chose nwf = ",I3)') igroup, wf_groups( igroup)%nwf
              end if
            end if
            ! check compatibility of number of Wannier functions with windows
            do ik = 1, wf_kset%nkpt
              if( wf_groups( igroup)%win_no( ik) + wf_groups( igroup)%win_ni( ik) .lt. wf_groups( igroup)%nwf) then
                if( mpiglobal%rank .eq. 0) then
                  write(*,*)
                  write( *, '("Error (wannier_readinput): Outer window contains less than nwf (",I3,") bands for k-point ",3F13.6," in group ",I2,".")') &
                      wf_groups( igroup)%nwf, wf_kset%vkl( :, ik), igroup
                end if
                stop
              end if
              if( wf_groups( igroup)%win_ni( ik) .gt. wf_groups( igroup)%nwf) then
                if( mpiglobal%rank .eq. 0) then
                  write(*,*)
                  write( *, '("Error (wannier_readinput): Inner window contains more than nwf (",I3,") bands for k-point ",3F13.6," in group ",I2,".")') &
                      wf_groups( igroup)%nwf, wf_kset%vkl( :, ik), igroup
                end if
                stop
              end if
            end do
          ! if no windows are defined use band ranges from input
          else
            wf_groups( igroup)%fst = min( input%properties%wannier%grouparray( igroup)%group%fst, input%properties%wannier%grouparray( igroup)%group%lst)
            wf_groups( igroup)%lst = max( input%properties%wannier%grouparray( igroup)%group%fst, input%properties%wannier%grouparray( igroup)%group%lst)
            wf_groups( igroup)%nwf = wf_groups( igroup)%lst - wf_groups( igroup)%fst + 1
            allocate( wf_groups( igroup)%win_ii( wf_groups( igroup)%nwf, wf_kset%nkpt), &
                      wf_groups( igroup)%win_io( wf_groups( igroup)%nwf, wf_kset%nkpt))
            allocate( wf_groups( igroup)%win_ni( wf_kset%nkpt), wf_groups( igroup)%win_no( wf_kset%nkpt))
            wf_groups( igroup)%win_ni = wf_groups( igroup)%nwf
            wf_groups( igroup)%win_no = 0
            wf_groups( igroup)%win_io = 0
            do ist = 1, wf_groups( igroup)%nwf
              wf_groups( igroup)%win_ii(ist,:) = ist + wf_groups( igroup)%fst - 1
            end do
          end if 
          ! sanity checks for band ranges
          if( wf_groups( igroup)%fst .lt. 1) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wannier_readinput): The lowest band (fst) is smaller than 1 for group ",I2,".")') igroup
            end if
            stop
          end if
          if( wf_groups( igroup)%lst .gt. nstfv) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wannier_readinput): The highest band (lst = ",I4,") is greater than total number of states (nstfv = ",I4,") for group ",I2,".")') &
                  wf_groups( igroup)%lst, nstfv, igroup
            end if
            stop
          end if
          ! set number of envolved bands
          wf_groups( igroup)%nst = wf_groups( igroup)%lst - wf_groups( igroup)%fst + 1
          ! check if band ranges are compatible with available states
          if( input%properties%wannier%input .eq. "gw") then
            if( (wf_groups( igroup)%fst .lt. input%gw%ibgw) .or. (wf_groups( igroup)%fst .gt. input%gw%nbgw)) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wannier_readinput): lower band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")') &
                    wf_groups( igroup)%fst, input%gw%ibgw, input%gw%nbgw-1, igroup
              end if
              stop
            end if
          else
            if( (wf_groups( igroup)%fst .lt. 1) .or. (wf_groups( igroup)%fst .gt. nstfv)) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wannier_readinput): lower band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")') &
                    wf_groups( igroup)%fst, 1, nstfv-1, igroup
              end if
              stop
            end if
          end if
          if( input%properties%wannier%input .eq. "gw") then
            if( wf_groups( igroup)%lst .gt. input%gw%nbgw) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wannier_readinput): upper band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")') &
                    wf_groups( igroup)%lst, wf_groups( igroup)%fst+1, input%gw%nbgw, igroup
              end if
              stop
            end if
          else
            if( wf_groups( igroup)%lst .gt. nstfv) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wannier_readinput): upper band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")') &
                    wf_groups( igroup)%lst, wf_groups( igroup)%fst+1, nstfv, igroup
              end if
              stop
            end if
          end if
          ! set the appropriate method
          if( wf_groups( igroup)%method .eq. 'pro' .or. &
              wf_groups( igroup)%method .eq. 'promax' .or. &
              wf_groups( igroup)%method .eq. 'opf' .or. &
              wf_groups( igroup)%method .eq. 'opfmax') then
            ! input is not compatible with selected method, change to suitable method
            if( wf_groups( igroup)%nst .ne. wf_groups( igroup)%nwf) then
              wf_groups( igroup)%method = 'auto'
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Warning (wannier_readinput): The selected method is not compatible with the provided input on the energy range for group ",I2,". The method was reset to auto.")') &
                  igroup
              end if
            end if
          end if
          if( wf_groups( igroup)%method .eq. 'auto' .or. &
              wf_groups( igroup)%method .eq. 'disSMV' .or. &
              wf_groups( igroup)%method .eq. 'disFull') then
            ! in this case there is nothing to disentangle and the effort can be reduced
            if( sum( wf_groups( igroup)%win_ni) .eq. wf_kset%nkpt*wf_groups( igroup)%nwf) then
              wf_groups( igroup)%method = 'opfmax'
              wf_groups( igroup)%fst = minval( wf_groups( igroup)%win_ii( 1:wf_groups( igroup)%win_ni(1), 1))
              wf_groups( igroup)%lst = maxval( wf_groups( igroup)%win_ii( 1:wf_groups( igroup)%win_ni(1), 1))
              wf_groups( igroup)%nst = wf_groups( igroup)%lst - wf_groups( igroup)%fst + 1
            else if( wf_groups( igroup)%method .eq. 'auto') then
              wf_groups( igroup)%method = 'disFull'
            end if
          end if

        end do !igroup

        wf_fst = 10000000
        wf_lst = 0
        wf_nwf = 0
        do igroup = 1, wf_ngroups
          wf_groups( igroup)%fwf = wf_nwf + 1
          wf_fst = min( wf_fst, wf_groups( igroup)%fst)
          wf_lst = max( wf_lst, wf_groups( igroup)%lst)
          wf_nwf = wf_nwf + wf_groups( igroup)%nwf
          wf_groups( igroup)%lwf = wf_nwf
        end do
        wf_nst = wf_lst - wf_fst + 1

      ! read input from TRANSFORM file
      else if( (input%properties%wannier%do .eq. "fromfile") .or. (input%properties%wannier%do .eq. "maxfromfile")) then

        call getunit( un)

        inquire( file=trim( wf_filename)//"_TRANSFORM"//trim( filext), exist=success)
        if( .not. success) then
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write(*, '"(Error (wannier_readinput): File '//trim( wf_filename)//"_TRANSFORM"//trim( filext)//' does not exist.")')
          end if
          return
        end if
        open( un, file=trim( wf_filename)//"_TRANSFORM"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
        read( un) fst_, lst_, nst_, nwf_, nkpt_, disentangle_
        !wf_disentangle = disentangle_
        if( (fst_ .ne. wf_fst) .or. (lst_ .ne. wf_lst)) then
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *, '("Warning (wannier_readinput): different band-ranges in input (",I4,":",I4,") and file (",I4,":",I4,").")') wf_fst, wf_lst, fst_, lst_
            write( *, '(" Use data from file.")')
          end if
          wf_fst = fst_
          wf_lst = lst_
          wf_nst = nst_
        end if
        if( nwf_ .ne. wf_nwf) then
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *, '("Warning (wannier_readinput): different number of Wannier functions in input (",I4,") and file (",I4,").")') wf_nwf, nwf_
            write( *, '(" Use data from file.")')
          end if
          wf_nwf = nwf_
        end if
        close( un)
      end if

      if( allocated( eval)) deallocate( eval)
      return
    end subroutine wannier_readinput

end module mod_wannier
