Subroutine checkwanint( b)
  Use modinput
  Use modmain
  Use modmpi
  use mod_wannier
  use m_wannier_interpolate_eigsys

  implicit none

  logical, intent( in) :: b

  integer :: ik, iq, is, ia, ias, l, m, lm, ngkmaxint, ngkmaxwf, nmatmaxwf
  real(8), allocatable :: eval1(:,:), evalfv(:,:), evalint(:,:)
  complex(8), allocatable :: evecint(:,:,:,:), evec1(:,:,:), evec2(:,:,:), apwalm1(:,:,:,:), apwalm2(:,:,:,:), dmat(:,:,:)
  type( k_set) :: int_kset
  type( G_set) :: int_Gset
  type( Gk_set) :: int_Gkset
  character(256) :: fname

  Call init0
  
  allocate( eval1( nstfv, nkptnr), evalfv( nstfv, nspinor))
  ! read FV eigenvalues
  do ik = 1, wf_kset%nkpt
    call getevalfv( wf_kset%vkl( :, ik), evalfv)
    eval1( :, ik) = evalfv( :, 1)
  end do
  
  ! k-points for interpolation
  call generate_k_vectors( int_kset, bvec, (/8, 8, 8/), (/0.d0, 0.d0, 0.d0/), .false.)
  !call generate_G_vectors( int_Gset, bvec, intgv, input%groundstate%gmaxvr)
  call generate_Gk_vectors( int_Gkset, int_kset, wf_Gset, gkmax)
  write(*,*) int_kset%nkpt
  write(*,*) shape( int_Gkset%ngk)
  ngkmaxint = int_Gkset%ngkmax
  ngkmaxwf = wf_Gkset%ngkmax
  ngkmax = ngkmaxwf
  nmatmaxwf = ngkmaxwf + nlotot
  nmatmax = ngkmax + nlotot
  write(*, '("ngkmax wan: ",I5.5)') ngkmaxwf
  write(*, '("ngkmax int: ",I5.5)') ngkmaxint
  write(*, '("ngkmax sys: ",I5.5)') ngkmax

  if( b) then
    !allocate( evalint( wf_fst:wf_lst, int_kset%nkpt))
    !allocate( evecint( ngkmaxint+nlotot, wf_fst:wf_lst, nspinor, int_kset%nkpt))
    !call wannier_interpolate_eigsys( eval1( wf_fst:wf_lst, :), int_kset, int_Gkset, evalint( wf_fst:wf_lst, :), evecint( :, :, 1, :))
    !do iq = 1, int_kset%nkpt
    !  write( fname, '("evecint/evecint_",I3.3,"_",I3.3,"_",I3.3)') nint( int_kset%vkl( :, iq)*1000)
    !  call writematlab( evecint( 1:(int_Gkset%ngk( 1, iq)+nlotot), wf_fst:wf_lst, 1, iq), fname)
    !end do
    !deallocate( evecint)
    call readstate
    call linengy
    call genapwfr
    call genlofr
    allocate( apwalm1( ngkmax, apwordmax, lmmaxapw, natmtot))
    allocate( apwalm2( ngkmax, apwordmax, lmmaxapw, natmtot))
    allocate( evec1( nmatmaxwf, nstsv, nspinor))
    allocate( evec2( nmatmaxwf, nstsv, nspinor))
    allocate( dmat( (input%groundstate%lmaxapw+1)**2, wf_fst:wf_lst, wf_fst:wf_lst))
    do ik = 1, wf_kset%nkpt
      if( input%properties%wannier%input .eq. "groundstate") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evec1)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evec1)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmaxwf, nstfv, nspinor, evec1)
      else
        call terminate
      end if
      call match( wf_Gkset%ngk( 1, ik), &
                  wf_Gkset%gkc( :, 1, ik), &
                  wf_Gkset%tpgkc( :, :, 1, ik), &
                  wf_Gkset%sfacgk( :, :, 1, ik), &
                  apwalm1)
      do iq = 1, ik
        write(*,*) ik, iq
        if( input%properties%wannier%input .eq. "groundstate") then
          call getevecfv( wf_kset%vkl( :, iq), wf_Gkset%vgkl( :, :, :, iq), evec2)
        else if( input%properties%wannier%input .eq. "hybrid") then
          call getevecfv( wf_kset%vkl( :, iq), wf_Gkset%vgkl( :, :, :, iq), evec2)
        else if( input%properties%wannier%input .eq. "gw") then
          call getevecsvgw_new( "GW_EVECSV.OUT", iq, wf_kset%vkl( :, iq), nmatmaxwf, nstfv, nspinor, evec2)
        else
          call terminate
        end if
        call match( wf_Gkset%ngk( 1, iq), &
                    wf_Gkset%gkc( :, 1, iq), &
                    wf_Gkset%tpgkc( :, :, 1, iq), &
                    wf_Gkset%sfacgk( :, :, 1, iq), &
                    apwalm2)
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            call mygendmat( input%groundstate%lmaxapw, &
                            wf_fst, wf_lst, &
                            is, ia, &
                            wf_Gkset%ngk( 1, ik), &
                            wf_Gkset%ngk( 1, iq), &
                            apwalm1, apwalm2, &
                            evec1( :, :, 1), evec2( :, :, 1), &
                            dmat)
            do l = 0, input%groundstate%lmaxapw
              do m = -l, l
                lm = idxlm( l, m)
                write( fname, '("rmat/rmat",4("_",I3.3))') ik, iq, lm, ias
                call writematlab( dmat( lm, :, :), fname)
              end do
            end do
          end do
        end do
      end do
    end do 
  else
    allocate( evecint( nmatmaxwf, nstsv, nspinor, 1))
    do ik = 1, wf_kset%nkpt
      if( input%properties%wannier%input .eq. "groundstate") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecint( :, :, :, 1))
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecint( :, :, :, 1))
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmaxwf, nstfv, nspinor, evecint( :, :, :, 1))
      else
        call terminate
      end if
      write( fname, '("eveccal/eveccal_",I3.3,"_",I3.3,"_",I3.3)') nint( wf_kset%vkl( :, ik)*1000)
      call writematlab( evecint( 1:(wf_Gkset%ngk( 1, ik)+nlotot), wf_fst:wf_lst, 1, 1), fname)
    end do 
  end if

  return
End Subroutine checkwanint
!EOC
