module mod_wannier
  use mod_wannier_variables
  use mod_wannier_projection
  use mod_wannier_maxloc
  use mod_wannier_opf
  use mod_wannier_scdm
  use mod_wannier_subspace
  use mod_wannier_filehandling
  use mod_wannier_helper
  use mod_optkgrid
  use mod_kpointset

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
      integer :: ist, ik, igroup

      if( wf_initialized) call wannier_destroy
      call init0
      call init1

      !********************************************************************
      ! open INFO file
      !********************************************************************
      write( wf_filename, '("WANNIER")')
      call getunit( wf_info)
      if( input%properties%wannier%do .ne. "fromfile") then
        open( wf_info, file=trim( wf_filename)//"_INFO"//trim(filext), action='write', form='formatted')
      end if
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
      
      if( input%properties%wannier%do .ne. "fromfile") call wannier_writeinfo_geometry

      !********************************************************************
      ! find idices for phase correction
      !********************************************************************
      !select case (input%properties%wannier%method)
      !  case( "fromfile")
      !    call wannier_readsetup
      !  case( "maxfromfile")
      !    call wannier_readsetup
      !  case default
      !    call wannier_getevecphase
      !end select
      !do ik = 1, wf_kset%nkpt
      !  write(*,'(I5,5x,1000I4)') ik, wf_evecphase( :, ik)
      !end do

      call wannier_genradfun

      !********************************************************************
      ! calculate plane-wave matrix elements
      !********************************************************************
      select case( input%properties%wannier%do)
        ! dont do
        case( "skip")
        case( "fromfile")
        ! do
        case default
          call printbox( wf_info, '*', "Preparation")
          write( wf_info, *)
          call wannier_emat
      end select
        
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

      !********************************************************************
      ! initialize transformation matrices
      !********************************************************************
      allocate( wf_transform( wf_fst:wf_lst, wf_nwf, wf_kset%nkpt))
      wf_transform = zzero
      ik = 0
      do igroup = 1, wf_ngroups
        do ist = 1, wf_groups( igroup)%nwf
          wf_transform( wf_groups( igroup)%fst+ist-1, ik+ist, :) = zone
        end do
        ik = ik + wf_groups( igroup)%nwf
      end do

      wf_initialized = .true.
      !call wannier_writesetup

      return
    end subroutine wannier_init
    !EOC

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
        call wannier_reademat( success)
        if( success) then
          call timesec( t1)
          write( wf_info, '(5x,"plane-wave matrix-elements read from file")')
          write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
          write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
          write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
        end if
      end if

      if( .not. success) then          
        if( allocated( wf_m0)) deallocate( wf_m0)
        allocate( wf_m0( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt, wf_n_ntot))

        wf_m0(:,:,:,:) = zzero

        !write(*,'("shape of M0:     ",4I6)') shape( wf_pwmat)
        !write(*,'("shape of evec:   ",3I6)') nmatmax_ptr, nstfv, nspinor
        !write(*,'("size of evec:    ",F13.6," GB")') nmatmax_ptr*nstfv*nspinor*16*1.d-9
        !write(*,'("shape of apwalm: ",4I6)') ngkmax_ptr, apwordmax, lmmaxapw, natmtot
        !write(*,'("size of apwalm:  ",F13.6," GB")') ngkmax_ptr*apwordmax*lmmaxapw*natmtot*16*1.d-9
        
        call pwmat_init( input%groundstate%lmaxapw, 8)

        do idxn = 1, wf_n_ntot
          call pwmat_init_qg( wf_n_vl( :, idxn), (/0, 0, 0/))
          cntk = 0
          !write(*,'(1a1,"  neighbor ",i2," of ",i2,": ",i3,"%",$)') char(13), idxn, wf_n_ntot, nint( 100.d0*cntk/wf_kset%nkpt)
#ifdef USEOMP
!!$omp parallel default( shared) private( iknr, evecfv1, evecfv2)
#endif
          allocate( evecfv1( nmatmax_ptr, nstfv, nspinor))
          allocate( evecfv2( nmatmax_ptr, nstfv, nspinor))
#ifdef USEOMP
!!$omp do
#endif
          do iknr = 1, wf_kset%nkpt
            k1 = iknr
            k2 = wf_n_ik( idxn, k1)
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
                 wf_m0( :, :, iknr, idxn))
#ifdef USEOMP
!!$omp atomic update
#endif
            cntk = cntk + 1
#ifdef USEOMP
!!$omp end atomic
#endif
            !write(*,'(1a1,"o neighbor ",i2," of ",i2,": ",i3,"%",$)') char(13), idxn, wf_n_ntot, nint( 100.d0*cntk/wf_kset%nkpt)
          end do
#ifdef USEOMP
!!$omp end do
#endif
          deallocate( evecfv1, evecfv2)
#ifdef USEOMP
!!$omp end parallel
#endif
          !write(*,*)
        end do

        !write(*,*) "destroy"
        call pwmat_destroy
        !write(*,*) "write"
        call wannier_writeemat
        call timesec( t1)
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
        write( wf_info, '(5x,"#neighbors per k-point: ",T40,7x,I6)') wf_n_ntot
      end if

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
      integer :: iknr, iproj, n2
      real(8) :: t0, t1
      real(8), allocatable :: sval(:)
      complex(8), allocatable :: projection(:,:), lsvec(:,:), rsvec(:,:)

      write( wf_info, '(" perform simple projection step...")')
      call timesec( t0)

      !********************************************************************
      ! build projection matrices
      !********************************************************************
      allocate( projection( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%nprojused))

      !********************************************************************
      ! build transformation matrices
      !********************************************************************

      allocate( sval( wf_groups( wf_group)%nwf), &
                lsvec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst), &
                rsvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf))
         
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, n2, iproj, projection, sval, lsvec, rsvec)
!$OMP DO  
#endif
      do iknr = 1, wf_kset%nkpt
        n2 = 0
        do iproj = 1, wf_nprojtot
          if( wf_groups( wf_group)%projused( iproj) .eq. 1) then
            n2 = n2 + 1
            projection( :, n2) = wf_projection( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, iproj, iknr)
          end if
        end do
        call zsvd( projection, sval, lsvec, rsvec)
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
      call wannier_loc
      deallocate( projection, sval, lsvec, rsvec)
      call timesec( t1)

      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      write( wf_info, *)
      call flushifc( wf_info)

      return
    end subroutine wannier_gen_pro
    !EOC

!    subroutine wannier_projonwan
!      !!!!!
!      ! correct wf_nst and wf_nwf here
!      !!!!!
!      use m_wsweight
!      use mod_Gvector,              only: cfunig 
!      use mod_eigensystem,          only: oalo, ololo
!      integer :: ik, is, ia, ias, lmaxapw, l, l_, m, m_, lm, lm_, o, o_, lmo, lmo_, ir, ilo, ilo_, ngknr, ist, jst, ig, igk, ifg, nlmomax, nshell, nrpt
!      real(8) :: x, fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)
!      
!      integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:)
!      real(8), allocatable :: rptl(:,:)
!      complex(8), allocatable :: wanfmt(:,:,:,:), wanfir(:,:,:)
!      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:), match_combined(:,:), wfmt(:,:), wfir(:,:), cfun(:), zfft(:), wswgt(:)
!      complex(8), allocatable :: addmt(:,:), addir(:,:)
!
!      lmaxapw = input%groundstate%lmaxapw
!      nlmomax = (lmaxapw + 1)**2*apwordmax
!      nshell = 1
!      nrpt = (1 + 2*nshell)**3
!
!      allocate( rptl( 3, nrpt))
!      allocate( wswgt( nrpt))
!      ir = 0
!      do l = -nshell, nshell
!        do m = -nshell, nshell
!          do o = -nshell, nshell
!            ir = ir + 1
!            rptl( :, ir) = dble( (/o, m, l/))
!          end do
!        end do
!      end do
!
!      allocate( nlmo( nspecies))
!      allocate( lmo2l( nlmomax, nspecies), &
!                lmo2m( nlmomax, nspecies), &
!                lmo2o( nlmomax, nspecies))
!      allocate( wanfmt( nlmomax+nlotot, wf_fst:wf_lst, natmtot, nrpt))
!      allocate( wanfir( wf_Gset%ngrtot, wf_fst:wf_lst, nrpt))
!      
!      call readstate
!      call readfermi
!      call linengy
!      call genapwfr
!      call genlofr
!      call olprad
!      call genidxlo
!
!      call wannier_readfun( nshell, wanfmt, wanfir, nlmo, lmo2l, lmo2m, lmo2o)
!
!      allocate( evecfv( nmatmax_ptr, nstfv, nspnfv))
!      allocate( apwalm( ngkmax_ptr, apwordmax, lmmaxapw, natmtot))
!      allocate( match_combined( nlmomax, ngkmax_ptr))
!      allocate( wfmt( nlmomax+nlotot, wf_fst:wf_lst))
!      allocate( wfir( wf_Gset%ngrtot, wf_fst:wf_lst))
!      allocate( cfun( wf_Gset%ngrtot))
!      allocate( zfft( wf_Gset%ngrtot))
!      allocate( addmt( wf_fst:wf_lst, wf_nst))
!      allocate( addir( wf_fst:wf_lst, wf_nst))
!
!      write(*,*) nlmomax, nlotot
!      write(*,*) shape( wanfmt)
!      write(*,*) shape( wfmt)
!
!      wf_projection = zzero
!      wf_projused = 0
!      wf_projused( 1:wf_nst) = 1
!
!      cfun = zzero
!      do ig = 1, wf_Gset%ngrtot
!        cfun( igfft( ig)) = cfunig( ig)
!      end do
!      call zfftifc( 3, ngrid, 1, cfun)
!
!      do ik = 1, wf_kset%nkpt
!        ngknr = wf_Gkset%ngk( 1, ik)
!        write(*,*) ik
!
!        do ir = 1, nrpt
!          call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), wswgt( ir), kgrid=.true.)
!        end do
!
!        ! get matching coefficients
!        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)
!          
!        ! read eigenvector      
!        if( input%properties%wannier%input .eq. "gs") then
!          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
!        else if( input%properties%wannier%input .eq. "hybrid") then
!          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
!        else if( input%properties%wannier%input .eq. "gw") then
!          call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax_ptr, nstfv, nspnfv, evecfv)
!        else
!          stop
!        end if
!
!        do is = 1, nspecies
!          do ia = 1, natoms( is)
!            ias = idxas( ia, is)
!
!            match_combined = zzero
!            wfmt = zzero
!
!            do lmo = 1, nlmo( is)
!              l = lmo2l( lmo, is)
!              m = lmo2m( lmo, is)
!              o = lmo2o( lmo, is)
!              lm = idxlm( l, m)
!              match_combined( lmo, 1:ngknr) = apwalm( 1:ngknr, o, lm, ias)
!            end do
!
!            call zgemm( 'N', 'N', nlmo( is), wf_nst, ngknr, zone, &
!                 match_combined( 1:nlmo( is), 1:ngknr), nlmo( is), &
!                 evecfv( 1:ngknr, wf_fst:wf_lst, 1), ngknr, zzero, &
!                 wfmt( 1:nlmo( is), :), nlmo( is))
!
!            do ilo = 1, nlorb( is)
!              l = lorbl( ilo, is)
!              do m = -l, l
!                lm = idxlm( l, m)
!                wfmt( nlmo( is)+idxlo( lm, ilo, ias), :) = evecfv( ngknr+idxlo( lm, ilo, ias), wf_fst:wf_lst, 1)
!              end do
!            end do
!
!            do lmo = 1, nlmo( is)
!              l = lmo2l( lmo, is)
!              m = lmo2m( lmo, is)
!              o = lmo2o( lmo, is)
!              lm = idxlm( l, m)
!              do ir = 1, nrpt
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!!$OMP DO
!#endif    
!                do ist = wf_fst, wf_lst
!                  do jst = wf_fst, wf_lst
!                    wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
!                        conjg( wfmt( lmo, ist))*wanfmt( lmo, jst, ias, ir)*wswgt( ir)
!                  end do
!                end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif    
!              end do
!              do ilo_ = 1, nlorb( is)
!                l_ = lorbl( ilo_, is)
!                do m_ = -l_, l_
!                  lm_ = idxlm( l_, m_)
!                  if( lm .eq. lm_) then
!                    do ir = 1, nrpt
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!!$OMP DO
!#endif    
!                      do ist = wf_fst, wf_lst
!                        do jst = wf_fst, wf_lst
!                          wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
!                              conjg( wfmt( lmo, ist))*wanfmt( nlmo( is)+idxlo( lm_, ilo_, ias), jst, ias, ir)*cmplx( oalo( o, ilo_, ias), 0, 8)*wswgt( ir)
!                        end do
!                      end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif    
!                    end do
!                  end if
!                end do
!              end do
!            end do
!            do ilo = 1, nlorb( is)
!              l = lorbl( ilo, is)
!              do m = -l, l
!                lm = idxlm( l, m)
!                do lmo_ = 1, nlmo( is)
!                  l_ = lmo2l( lmo_, is)
!                  m_ = lmo2m( lmo_, is)
!                  o_ = lmo2o( lmo_, is)
!                  lm_ = idxlm( l_, m_)
!                  if( lm .eq. lm_) then
!                    do ir = 1, nrpt
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!!$OMP DO
!#endif    
!                      do ist = wf_fst, wf_lst
!                        do jst = wf_fst, wf_lst
!                          wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
!                              conjg( wfmt( nlmo( is)+idxlo( lm, ilo, ias), ist))*wanfmt( lmo_, jst, ias, ir)*cmplx( oalo( o_, ilo, ias), 0, 8)*wswgt( ir)
!                        end do
!                      end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif    
!                    end do
!                  end if
!                end do
!                do ilo_ = 1, nlorb( is)
!                  l_ = lorbl( ilo_, is)
!                  do m_ = -l_, l_
!                    lm_ = idxlm( l_, m_)
!                    if( lm .eq. lm_) then
!                      do ir = 1, nrpt
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst)
!!$OMP DO
!#endif    
!                        do ist = wf_fst, wf_lst
!                          do jst = wf_fst, wf_lst
!                            wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
!                                conjg( wfmt( nlmo( is)+idxlo( lm, ilo, ias), ist))*wanfmt( nlmo( is)+idxlo( lm_, ilo_, ias), jst, ias, ir)*cmplx( ololo( ilo, ilo_, ias), 0, 8)*wswgt( ir)
!                          end do
!                        end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif    
!                      end do
!                    end if
!                  end do
!                end do
!              end do
!            end do                      
!              
!          end do
!        end do
!
!        wfir = zzero
!        do igk = 1, ngknr
!          ifg = igfft( wf_Gkset%igkig( igk, 1, ik))
!          wfir( ifg, :) = evecfv( igk, wf_fst:wf_lst, 1)
!        end do
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist)
!!$OMP DO  
!#endif
!        do ist = wf_fst, wf_lst
!          call zfftifc( 3, ngrid, 1, wfir( :, ist))
!        end do
!#ifdef USEOMP                
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
!
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ig, x, ifg)
!!$OMP DO  
!#endif
!        do ig = 1, wf_Gset%ngrtot
!          x = wf_kset%vkl( 1, ik)*wf_Gset%ivg( 1, ig)/ngrid(1) + &
!              wf_kset%vkl( 2, ik)*wf_Gset%ivg( 2, ig)/ngrid(2) + &
!              wf_kset%vkl( 3, ik)*wf_Gset%ivg( 3, ig)/ngrid(3)
!          ifg = igfft( ig)
!          wfir( ifg, :) = wfir( ifg, :)*cmplx( cos( twopi*x), sin( twopi*x), 8)
!        end do      
!#ifdef USEOMP                
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
!
!        do ir = 1, nrpt
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, jst, ig, zfft)
!!$OMP DO  
!#endif
!          do ist = wf_fst, wf_lst
!            do jst = wf_fst, wf_lst
!              do ig = 1, wf_Gset%ngrtot
!                zfft( ig) = conjg( wfir( ig, ist))*wanfir( ig, jst, ir)*cfun( ig)
!              end do
!              call zfftifc( 3, ngrid, -1, zfft)
!              wf_projection( ist, jst-wf_fst+1, ik) = wf_projection( ist, jst-wf_fst+1, ik) + &
!                  zfft( igfft( wf_Gset%ivgig( 0, 0, 0)))*wswgt( ir)
!            end do
!          end do
!#ifdef USEOMP                
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
!        end do
!
!      end do   
!    
!      deallocate( wanfmt, wanfir, wfmt, wfir, evecfv, apwalm, match_combined, cfun, zfft, nlmo, lmo2l, lmo2m, lmo2o)
!
!      return
!    end subroutine wannier_projonwan
    
    subroutine wannier_gen_fromfile
      logical :: success

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      call wannier_readtransform( success)
      if( .not. success) then
        write(*, '(" Failed to read Wannier functions from file. Recalculate them.")')
        input%properties%wannier%do = "fromscratch"
        call wannier_init
        do wf_group = 1, wf_ngroups
          call wannier_gen_opf
          call wannier_maxloc
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
    
    !BOP
    ! !ROUTINE: wannier_geometry
    ! !INTERFACE:
    !
    subroutine wannier_geometry
      ! !USES:
      use m_linalg
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

      integer :: ix, iy, iz, i, j, nshell, nkmax, nposvec, iknr, d, nbzshell(3), minshell
      integer :: vi(3)
      real(8) :: vl(3), vc(3), dist
      real(8) :: nvlt( 3, wf_kset%nkpt, wf_kset%nkpt), nvct( 3, wf_kset%nkpt, wf_kset%nkpt)
      real(8) :: coeff(6,6), right(6), sval(6), coeffinv(6,6)
      real(8) :: mat1(3,3), mat2(3,3), mat3(3,3), umat(3,3), u(3), n(3), nvec(3,12),  c, s
      logical :: stopshell, cutshell

      real(8), allocatable :: work(:), dt(:)
      integer, allocatable :: iwork(:), mt(:), posvec(:,:)
      
      integer, allocatable :: tmp_n_n(:)
      real(8), allocatable :: tmp_n_dist(:), tmp_n_vl(:,:,:), tmp_n_vc(:,:,:)

      nbzshell = input%properties%wannier%nbzshell
      minshell = input%properties%wannier%minshell
      cutshell = input%properties%wannier%cutshell
      
      ! find possible distances (shells)
      nposvec = (nbzshell(3)*wf_kset%ngridk(3)-1)*(2*nbzshell(2)*wf_kset%ngridk(2)-1)*(2*nbzshell(1)*wf_kset%ngridk(1)-1) + &
                (nbzshell(2)*wf_kset%ngridk(2)-1)*(2*nbzshell(1)*wf_kset%ngridk(1)-1) + &
                 nbzshell(1)*wf_kset%ngridk(1)-1
      allocate( dt( nposvec))
      allocate( mt( nposvec))
      allocate( posvec( 3, nposvec))

      i = 0
      posvec(:,:) = 0
      do iz = nbzshell(3)*wf_kset%ngridk(3)-1, -nbzshell(3)*wf_kset%ngridk(3)+1, -1
        do iy = nbzshell(2)*wf_kset%ngridk(2)-1, -nbzshell(2)*wf_kset%ngridk(2)+1, -1
          do ix = nbzshell(1)*wf_kset%ngridk(1)-1, -nbzshell(1)*wf_kset%ngridk(1)+1, -1
            vi = (/ix, iy, iz/)
            if( (abs(ix)+abs(iy)+abs(iz)) .ne. 0) then
              stopshell = .true.
              d = maxloc( abs( vi), 1)
              do j = 1, i
                if( cutshell) then
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
      nshell = iz

      ! sort shells by distance
      allocate( tmp_n_dist( nshell))
      allocate( tmp_n_n( nshell))

      tmp_n_dist = 0.d0
      tmp_n_n = 0
      dist = minval( dt)
      do i = 1, nshell
        j = minloc( abs( dt( 1:iz) - dist), 1)
        tmp_n_dist( i) = dt( j)
        tmp_n_n( i) = mt( j)
        dt( j) = 1.d100 
        dist = minval( dt)
      end do
      nkmax = maxval( tmp_n_n)
      
      ! find all possible neighbors
      allocate( tmp_n_vl( 3, nkmax, nshell))
      allocate( tmp_n_vc( 3, nkmax, nshell))

      tmp_n_vl = 0.d0
      tmp_n_vc = 0.d0
      mt = 0
      do ix = 1, nposvec
        vl = dble( posvec(:,ix))/wf_kset%ngridk
        call r3mv( bvec, vl, vc)
        dist = norm2( vc)
        j = minloc( abs( tmp_n_dist(:) - dist), 1)
        if( abs( tmp_n_dist( j) - dist) .lt. input%structure%epslat) then
          mt( j) = mt( j) + 1
          tmp_n_vl( :, mt( j), j) = vl
          tmp_n_vc( :, mt( j), j) = vc
        end if
      end do

      ! find number of shells needed for gradient calculation and geometric weights
      if( .not. allocated( wf_n_wgts))          allocate( wf_n_wgts( nshell))
      if( .not. allocated( wf_n_usedshells))    allocate( wf_n_usedshells( nshell))

      j = 1
      stopshell = .false.
      coeff = 0.d0
      wf_n_wgts = 0.d0
      allocate( work( 1), iwork( 8*j))

      if( sum( wf_kset%ngridk) .eq. 1) then
        write(*,'(" Error (wannier_geometry): Wannier Functions for single k-point only not yet implemented.")')
        stop
      else if( minval( wf_kset%ngridk) .eq. 1) then
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
      wf_n_rot = 0.d0
      wf_n_rot(1,1) = 1.d0
      wf_n_rot(2,2) = 1.d0
      wf_n_rot(3,3) = 1.d0

      if( all( nbzshell*wf_kset%ngridk .gt. 1)) d = 6
      wf_n_usedshells = 0
      wf_n_ntot = 0
      nvec = 0.d0
      wf_n_nshells = 0

      if( d .gt. 1) then
        if( d .eq. 3) then
          ix = minloc( wf_kset%ngridk, 1)
          iy = mod( ix+1, 3)
          iz = mod( ix+2, 3)
          if( iy .eq. 0) iy = 3
          if( iz .eq. 0) iz = 3

          ! rotate into xy-plane
          mat3 = 0.d0
          mat3(1,1) = 1.d0
          mat3(2,2) = 1.d0
          mat3(3,3) = 1.d0
          call r3cross( bvec( :, iy), bvec( :, iz), n)
          n = n/norm2( n)
          if( norm2( n - (/0.d0, 0.d0, 1.d0/)) .gt. input%structure%epslat) then
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
            call r3mm( umat, umat, wf_n_rot)
            wf_n_rot = (1-c)*wf_n_rot + s*umat + mat3
          end if
          mat3(3,3) = 0.d0
          call r3mm( mat3, wf_n_rot, umat)
          call r3mm( transpose( wf_n_rot), umat, mat3)
          call r3mm( mat1, mat3, umat)
          call r3mm( umat, mat2, mat3)
          right(:) = (/mat3(1,1), mat3(2,2), mat3(1,2), 0.d0, 0.d0, 0.d0/)

        else
          right(:) = (/mat3(1,1), mat3(2,2), mat3(3,3), mat3(1,2), mat3(1,3), mat3(2,3)/)
        end if
        do while( (wf_n_nshells .lt. minshell) .or. (.not. stopshell))
          ! check if new shell is linear dependend
          wf_n_usedshells( j) = 1
          !write(*,'(i,3f13.6)') j, tmp_n_vl( :, 1, j) 
          !if( wf_n_ntot .gt. 0) then
          !  call rpinv( nvec( :, 1:wf_n_ntot), nveci( 1:wf_n_ntot, :))
          !  wf_n_usedshells( j) = 0
          !  do i = 1, tmp_n_n( j)
          !    if( norm2( tmp_n_vl( :, i, j) - matmul( matmul( nvec( :, 1:wf_n_ntot), nveci( 1:wf_n_ntot, :)), tmp_n_vl( :, i, j))) &
          !    .gt. input%structure%epslat) wf_n_usedshells( j) = 1
          !  end do
          !end if              

          if( wf_n_usedshells( j)) then
            wf_n_nshells = wf_n_nshells + 1
            if( d .eq. 3) then
              do i = 1, tmp_n_n( j)
                coeff( 1, wf_n_nshells) = coeff( 1, wf_n_nshells) + tmp_n_vl( iy, i, j)**2
                coeff( 2, wf_n_nshells) = coeff( 2, wf_n_nshells) + tmp_n_vl( iz, i, j)**2
                coeff( 3, wf_n_nshells) = coeff( 3, wf_n_nshells) + tmp_n_vl( iy, i, j)*tmp_n_vl( iz, i, j)
                wf_n_ntot = wf_n_ntot + 1
                nvec( :, wf_n_ntot) = tmp_n_vl( :, i, j)
              end do
            else
              do i = 1, tmp_n_n( j)
                coeff( 1, wf_n_nshells) = coeff( 1, wf_n_nshells) + tmp_n_vl( 1, i, j)**2
                coeff( 2, wf_n_nshells) = coeff( 2, wf_n_nshells) + tmp_n_vl( 2, i, j)**2
                coeff( 3, wf_n_nshells) = coeff( 3, wf_n_nshells) + tmp_n_vl( 3, i, j)**2
                coeff( 4, wf_n_nshells) = coeff( 4, wf_n_nshells) + tmp_n_vl( 1, i, j)*tmp_n_vl( 2, i, j)
                coeff( 5, wf_n_nshells) = coeff( 5, wf_n_nshells) + tmp_n_vl( 1, i, j)*tmp_n_vl( 3, i, j)
                coeff( 6, wf_n_nshells) = coeff( 6, wf_n_nshells) + tmp_n_vl( 2, i, j)*tmp_n_vl( 3, i, j)
                wf_n_ntot = wf_n_ntot + 1
                nvec( :, wf_n_ntot) = tmp_n_vl( :, i, j)
              end do
            end if
            call rpinv( coeff( 1:d, 1:wf_n_nshells), coeffinv( 1:wf_n_nshells, 1:d))
            sval( 1:d) = matmul( matmul( coeff( 1:d, 1:wf_n_nshells), coeffinv( 1:wf_n_nshells, 1:d)), right( 1:d))
            if( (norm2( sval( 1:d) - right( 1:d)) .lt. input%structure%epslat) .or. (j .eq. nshell)) then
              stopshell = .true.
              wf_n_wgts( 1:wf_n_nshells) = matmul( coeffinv( 1:wf_n_nshells, 1:d), right( 1:d))
            end if
          else
            write(*,'("shell ",i," is linear dependend")') j
          end if
          j = j+1
        end do
      else
        wf_n_wgts( 1) = 1.0d0/norm2( tmp_n_vc(:,1,1))**2
        wf_n_nshells = 1
      end if

      wf_n_wgts = 0.5d0*wf_n_wgts

      ! find k-point indices of neighbors and copy results to global arrays
      if( .not. allocated( wf_n_vl))            allocate( wf_n_vl( 3, wf_n_ntot))
      if( .not. allocated( wf_n_vc))            allocate( wf_n_vc( 3, wf_n_ntot))
      if( .not. allocated( wf_n_ik))            allocate( wf_n_ik( wf_n_ntot, wf_kset%nkptnr))
      if( .not. allocated( wf_n_ik2))           allocate( wf_n_ik2( wf_n_ntot, wf_kset%nkptnr))
      if( .not. allocated( wf_n_wgt))           allocate( wf_n_wgt( wf_n_ntot))
      if( .not. allocated( wf_n_dist))          allocate( wf_n_dist( wf_n_ntot))

      d = 0
      wf_n_nshells = 0
      do j = 1, nshell
        if( wf_n_usedshells( j)) then
          wf_n_nshells = wf_n_nshells + 1
          do i = 1, tmp_n_n( j)
            d = d + 1
            do iknr = 1, wf_kset%nkpt
              wf_n_vl( :, d) = tmp_n_vl( :, i, j)
              wf_n_vc( :, d) = tmp_n_vc( :, i, j)
              wf_n_wgt( d) = wf_n_wgts( wf_n_nshells)
              wf_n_dist( d) = tmp_n_dist( j)

              vl = wf_kset%vkl( :, iknr) + wf_n_vl( :, d)
              call r3frac( input%structure%epslat, vl, vi) 
              call findkptinset( vl, wf_kset, ix, iy)
              if( ix .gt. 0) then
                wf_n_ik( d, iknr) = iy
              else
                write( *, '(" Error (wannier_geometry): wrong neighboring vector.")')
                write( *, '(3F13.6)') wf_n_vl( :, d)
                stop
              end if

              vl = wf_kset%vkl( :, iknr) - wf_n_vl( :, d)
              call r3frac( input%structure%epslat, vl, vi) 
              call findkptinset( vl, wf_kset, ix, iy)
              if( ix .gt. 0) then
                wf_n_ik2( d, iknr) = iy
              else
                write( *, '(" Error (wannier_geometry): wrong neighboring vector.")')
                write( *, '(3F13.6)') wf_n_vl( :, d)
                stop
              end if
            end do
          end do
        end if
      end do
      
      deallocate( dt, mt, work, iwork, posvec)
      deallocate( tmp_n_dist, tmp_n_n, tmp_n_vl, tmp_n_vc)

      return
      !EOC
    end subroutine wannier_geometry
    !EOP
    
    subroutine wannier_readinput
      use mod_eigenvalue_occupancy, only: efermi
      integer :: igroup, ik, fst, lst, ist, un
      integer :: fst_, lst_, nst_, nwf_, nkpt_
      logical :: disentangle_, success
      real(8) :: e
      real(8), allocatable :: evalfv(:,:)
      character(256) :: fxt

      wf_fermizero = input%properties%wannier%fermizero
      wf_ngroups = size( input%properties%wannier%grouparray, 1)
      allocate( wf_groups( wf_ngroups))

      ! read input from input.xml
      if( input%properties%wannier%do .ne. "skip") then
        do igroup = 1, wf_ngroups
          wf_groups( igroup)%method = input%properties%wannier%grouparray( igroup)%group%method
          wf_groups( igroup)%win_o(1) = minval( input%properties%wannier%grouparray( igroup)%group%outerwindow)
          wf_groups( igroup)%win_o(2) = maxval( input%properties%wannier%grouparray( igroup)%group%outerwindow)

          if( norm2( wf_groups( igroup)%win_o) .gt. 1.d-20) then
            wf_groups( igroup)%win_i(1) = minval( input%properties%wannier%grouparray( igroup)%group%innerwindow)
            wf_groups( igroup)%win_i(2) = maxval( input%properties%wannier%grouparray( igroup)%group%innerwindow)

            if( norm2( wf_groups( igroup)%win_i) .gt. input%structure%epslat) then
              if( (wf_groups( igroup)%win_i(1) .lt. wf_groups( igroup)%win_o(1)) .or. &
                  (wf_groups( igroup)%win_i(2) .gt. wf_groups( igroup)%win_o(2))) then
                write( *, '(" Error (wannier_readinput): The inner window must be fully contained within the outer window for group ",I2,".")') igroup
                stop
              end if
            end if

            fxt = filext
            if( input%properties%wannier%input .eq. "gw") write( filext, '("_GW.OUT")')
            call readfermi
            filext = fxt
            call wannier_geteval( evalfv, fst, lst)
            allocate( wf_groups( igroup)%win_ii( lst-fst+1, wf_kset%nkpt), wf_groups( igroup)%win_io( lst-fst+1, wf_kset%nkpt))
            allocate( wf_groups( igroup)%win_ni( wf_kset%nkpt), wf_groups( igroup)%win_no( wf_kset%nkpt))
            wf_groups( igroup)%win_ii = 0
            wf_groups( igroup)%win_io = 0
            wf_groups( igroup)%win_ni = 0
            wf_groups( igroup)%win_no = 0

            wf_groups( igroup)%fst = 100000000
            wf_groups( igroup)%lst = 0
            do ik = 1, wf_kset%nkpt
              do ist = fst, lst
                e = evalfv( ist, ik)
                if( wf_fermizero) e = e - efermi
                if( (e .ge. wf_groups( igroup)%win_i(1)) .and. (e .le. wf_groups( igroup)%win_i(2))) then
                  wf_groups( igroup)%win_ni( ik) = wf_groups( igroup)%win_ni( ik) + 1
                  wf_groups( igroup)%win_ii( wf_groups( igroup)%win_ni( ik), ik) = ist
                else if( (e .ge. wf_groups( igroup)%win_o(1)) .and. (e .le. wf_groups( igroup)%win_o(2))) then
                  wf_groups( igroup)%win_no( ik) = wf_groups( igroup)%win_no( ik) + 1
                  wf_groups( igroup)%win_io( wf_groups( igroup)%win_no( ik), ik) = ist
                end if
              end do
              wf_groups( igroup)%fst = min( wf_groups( igroup)%fst, minval( wf_groups( igroup)%win_ii( 1:wf_groups( igroup)%win_ni( ik), ik)))
              wf_groups( igroup)%fst = min( wf_groups( igroup)%fst, minval( wf_groups( igroup)%win_io( 1:wf_groups( igroup)%win_no( ik), ik)))
              wf_groups( igroup)%lst = max( wf_groups( igroup)%lst, maxval( wf_groups( igroup)%win_ii( 1:wf_groups( igroup)%win_ni( ik), ik)))
              wf_groups( igroup)%lst = max( wf_groups( igroup)%lst, maxval( wf_groups( igroup)%win_io( 1:wf_groups( igroup)%win_no( ik), ik)))
            end do
            wf_groups( igroup)%nwf = input%properties%wannier%grouparray( igroup)%group%nwf
            if( wf_groups( igroup)%nwf .lt. 1) then
              wf_groups( igroup)%nwf = nint( 0.5d0*(maxval( wf_groups( igroup)%win_ni) + minval( wf_groups( igroup)%win_ni + wf_groups( igroup)%win_no)))
              write( *, '(" Warning (wannier_readinput): Number of Wannier functions (nwf) not set for group ",I2,". I chose nwf = ",I3)') igroup, wf_groups( igroup)%nwf
            end if

            do ik = 1, wf_kset%nkpt
              if( wf_groups( igroup)%win_no( ik) + wf_groups( igroup)%win_ni( ik) .lt. wf_groups( igroup)%nwf) then
                write( *, '(" Error (wannier_readinput): Outer window contains less than nwf (",I3,") bands for k-point ",3F13.6," in group ",I2,".")') &
                    wf_groups( igroup)%nwf, wf_kset%vkl( :, ik), igroup
                stop
              end if
              if( wf_groups( igroup)%win_ni( ik) .gt. wf_groups( igroup)%nwf) then
                write( *, '(" Error (wannier_readinput): Inner window contains more than nwf (",I3,") bands for k-point ",3F13.6," in group ",I2,".")') &
                    wf_groups( igroup)%nwf, wf_kset%vkl( :, ik), igroup
                stop
              end if
            end do

          else
            wf_groups( igroup)%fst = min( input%properties%wannier%grouparray( igroup)%group%fst, input%properties%wannier%grouparray( igroup)%group%lst)       
            wf_groups( igroup)%lst = max( input%properties%wannier%grouparray( igroup)%group%fst, input%properties%wannier%grouparray( igroup)%group%lst)       
            wf_groups( igroup)%nwf = wf_groups( igroup)%lst - wf_groups( igroup)%fst + 1
          end if  

          if( wf_groups( igroup)%fst .lt. 1) then
            write( *, '(" Error (wannier_readinput): The lowest band (fst) is smaller than 1 for group ",I2,".")') igroup
            stop
          end if

          if( wf_groups( igroup)%lst .gt. nstfv) then
            write( *, '(" Error (wannier_readinput): The highest band (lst = ",I4,") is greater than total number of states (nstfv = ",I4,") for group ",I2,".")') &
                wf_groups( igroup)%lst, nstfv, igroup
            stop
          end if

          wf_groups( igroup)%nst = wf_groups( igroup)%lst - wf_groups( igroup)%fst + 1

          if( input%properties%wannier%input .eq. "gw") then
            if( (wf_groups( igroup)%fst .lt. input%gw%ibgw) .or. (wf_groups( igroup)%fst .gt. input%gw%nbgw)) then
              write( *, '(" Error (wannier_readinput): lower band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")'), &
                  wf_groups( igroup)%fst, input%gw%ibgw, input%gw%nbgw-1, igroup
              stop
            end if
          else
            if( (wf_groups( igroup)%fst .lt. 1) .or. (wf_groups( igroup)%fst .gt. nstfv)) then
              write( *, '(" Error (wannier_readinput): lower band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")'), &
                  wf_groups( igroup)%fst, 1, nstfv-1, igroup
              stop
            end if
          end if

          if( input%properties%wannier%input .eq. "gw") then
            if( wf_groups( igroup)%lst .gt. input%gw%nbgw) then
              write( *, '(" Error (wannier_readinput): upper band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")'), &
                  wf_groups( igroup)%lst, wf_groups( igroup)%fst+1, input%gw%nbgw, igroup
              stop
            end if
          else
            if( wf_groups( igroup)%lst .gt. nstfv) then
              write( *, '(" Error (wannier_readinput): upper band-index (",I4,") out of range (",I4,":",I4,") for group ",I2,".")'), &
                  wf_groups( igroup)%lst, wf_groups( igroup)%fst+1, nstfv, igroup
              stop
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
          write(*,*) 'Error (wannier_readinput): File '//trim( wf_filename)//"_TRANSFORM"//trim( filext)//' does not exist.'
          return
        end if
        open( un, file=trim( wf_filename)//"_TRANSFORM"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
        read( un) fst_, lst_, nst_, nwf_, nkpt_, disentangle_
        !wf_disentangle = disentangle_
        if( (fst_ .ne. wf_fst) .or. (lst_ .ne. wf_lst)) then
          write( *, '(" Warning (wannier_readinput): different band-ranges in input (",I4,":",I4,") and file (",I4,":",I4,").")'), wf_fst, wf_lst, fst_, lst_
          write( *, '(" Use data from file.")')
          wf_fst = fst_
          wf_lst = lst_
          wf_nst = nst_
        end if
        if( nwf_ .ne. wf_nwf) then
          write( *, '(" Warning (wannier_readinput): different number of Wannier functions in input (",I4,") and file (",I4,").")') wf_nwf, nwf_
          write( *, '(" Use data from file.")')
          wf_nwf = nwf_
        end if
        close( un)
      end if

      return
    end subroutine wannier_readinput

    subroutine wannier_getevecphase
      integer :: ik, ist, idxp( wf_fst:wf_lst, wf_kset%nkpt)
      complex(8) :: evecfv( nmatmax_ptr, nstfv, nspinor)

      if( allocated( wf_evecphase)) deallocate( wf_evecphase)

      do ik = 1, wf_kset%nkpt
        call wannier_getevec( ik, evecfv)
        do ist = wf_fst, wf_lst
          idxp( ist, ik) = maxloc( abs( evecfv( 1:wf_Gkset%ngk( 1, ik), ist, 1)), 1)
        end do
      end do

      allocate( wf_evecphase( wf_fst:wf_lst, wf_kset%nkpt))
      wf_evecphase = idxp

      return
    end subroutine wannier_getevecphase

    subroutine wannier_updatephase
      integer :: ik, ist, idxp( wf_fst:wf_lst, wf_kset%nkpt)
      real(8) :: p1( wf_fst:wf_lst, wf_kset%nkpt), p2( wf_fst:wf_lst, wf_kset%nkpt)
      complex(8) :: evecfv( nmatmax_ptr, nstfv, nspinor)

      idxp = wf_evecphase
      if( allocated( wf_evecphase)) deallocate( wf_evecphase)

      do ik = 1, wf_kset%nkpt
        call wannier_getevec( ik, evecfv)
        do ist = wf_fst, wf_lst
          p1( ist, ik) = atan2( dble( aimag( evecfv( idxp( ist, ik), ist, 1))), dble( evecfv( idxp( ist, ik), ist, 1)))
        end do
      end do

      call wannier_getevecphase
      idxp = wf_evecphase
      if( allocated( wf_evecphase)) deallocate( wf_evecphase)

      do ik = 1, wf_kset%nkpt
        call wannier_getevec( ik, evecfv)
        do ist = wf_fst, wf_lst
          p2( ist, ik) = atan2( dble( aimag( evecfv( idxp( ist, ik), ist, 1))), dble( evecfv( idxp( ist, ik), ist, 1)))
        end do
      end do

      allocate( wf_evecphase( wf_fst:wf_lst, wf_kset%nkpt))
      wf_evecphase = idxp

      do ik = 1, wf_kset%nkpt
        do ist = wf_fst, wf_lst
          wf_transform( ist, :, ik) = wf_transform( ist, :, ik)*cmplx( cos( p2( ist, ik)-p1( ist, ik)), sin( p2( ist, ik)-p1( ist, ik)), 8)
        end do
      end do

      call wannier_writesetup
      call wannier_writetransform

      return
    end subroutine wannier_updatephase

end module mod_wannier
