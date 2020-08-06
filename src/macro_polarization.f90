subroutine macro_polarization( pol, mode)
  use modinput
  use modmpi
  use mod_pwmat
  use mod_kpointset
  use mod_lattice
  use mod_atoms
  use mod_constants
  use m_linalg, only: zdet
  use mod_charge_and_moment, only: chgval
  use mod_eigenvalue_occupancy, only: occmax, nstfv
  use mod_eigensystem, only: nmatmax_ptr 
  use mod_spin, only: nspinor
  use mod_Gvector, only: intgv
  use mod_Gkvector, only: gkmax

  real(8), intent( inout)  :: pol(3,3)
  character(*), intent( in):: mode

  integer :: ik1, ik2, d1, d2, d3, ix, iy, iz, nocc, nplane, is, ia, ist, vi(3)
  real(8) :: pel(3), pion(3), ptot(3), vr(3), phase, berrypel, berrypion, zion, vk(3), vkv(3)
  complex(8) :: det, pavg
  type( k_set) :: kset
  type( G_set) :: Gset
  type( Gk_set) :: Gkset
  
  complex(8), allocatable :: s(:,:), evec1(:,:,:), evec2(:,:,:), prod(:,:), prodv(:)

  if( trim( mode) == 'write') then
    call put( pol)
    return
  elseif( trim( mode) == 'read') then
    call get( pol)
    return
  end if

  call init0
  call init1

  nocc = nint( chgval/occmax)
  pel = 0.d0
  pion = 0.d0

  call generate_k_vectors( kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
  call generate_G_vectors( Gset, bvec, intgv, input%groundstate%gmaxvr)
  call generate_Gk_vectors( Gkset, kset, Gset, gkmax)

  allocate( evec1( nmatmax_ptr, nstfv, nspinor))
  allocate( evec2( nmatmax_ptr, nstfv, nspinor))
  allocate( s( nocc, nocc))

  call readstate
  call linengy
  call genapwfr
  call genlofr
  call olprad
  call pwmat_init( input%groundstate%lmaxapw, 8, kset, 1, nocc, 1, nocc)
  do ik1 = firstofset( mpiglobal%rank, kset%nkpt), lastofset( mpiglobal%rank, kset%nkpt)
    call getevecfv( kset%vkl( :, ik1), Gkset%vgkl( :, :, :, ik1), evec1)
    call pwmat_prepare( ik1, evec1( :, :, 1))
  end do
  call barrier

  do d1 = 1, 3
    if( d1 .eq. 1) then
      d2 = 2
      d3 = 3
    else if( d1 .eq. 2) then
      d2 = 1
      d3 = 3
    else
      d2 = 1
      d3 = 2
    end if
    nplane = kset%ngridk( d2)*kset%ngridk( d3)
    call r3cross( kset%bvec( :, d2), kset%bvec( :, d3), vr)
    vr = 0.d0
    vr(d1) = 1.d0/kset%ngridk( d1)
    call pwmat_init_qg( vr, (/0, 0, 0/), 1)
    if( allocated( prod)) deallocate( prod)
    allocate( prod( kset%ngridk( d2), kset%ngridk( d3)))
    if( allocated( prodv)) deallocate( prodv)
    allocate( prodv( kset%ngridk( d2)*kset%ngridk( d3)))

    prodv = zzero
    do ia = firstofset( mpiglobal%rank, kset%ngridk( d2)*kset%ngridk( d3)), lastofset( mpiglobal%rank, kset%ngridk( d2)*kset%ngridk( d3))
      iy = (ia-1)/kset%ngridk( d2)
      ix = ia - iy*kset%ngridk( d2) - 1
      prodv( ia) = zone
      do iz = 0, kset%ngridk( d1) - 1
        if( d1 .eq. 1) then
          vk = dble( (/iz, ix, iy/))/kset%ngridk
        else if( d1 .eq. 2) then
          vk = dble( (/ix, iz, iy/))/kset%ngridk
        else
          vk = dble( (/ix, iy, iz/))/kset%ngridk
        end if
        vkv = vk + vr
        call r3frac( input%structure%epslat, vkv, vi)
        call findkptinset( vk, kset, is, ik1)
        call getevecfv( kset%vkl( :, ik1), Gkset%vgkl( :, :, :, ik1), evec1)
        call findkptinset( vkv, kset, is, ik2)
        call getevecfv( kset%vkl( :, ik2), Gkset%vgkl( :, :, :, ik2), evec2)
        call pwmat_genpwmat( ik1, &
               evec1( :, 1:nocc, 1), &
               evec2( :, 1:nocc, 1), &
               s)
        call zdet( s, det)
        prodv( ia) = prodv( ia)*det
      end do
    end do
    call mpi_allgatherv_ifc( rlen=1, set=kset%ngridk( d2)*kset%ngridk( d3), zbuf=prodv)
    call barrier
    prod = reshape( prodv, (/kset%ngridk( d2), kset%ngridk( d3)/))
    pavg = zzero
    do ix = 1, kset%ngridk( d2)
      do iy = 1, kset%ngridk( d3)
        pavg = pavg + prod( ix, iy)
      end do
    end do
    pavg = pavg/dble( nplane)
    prod = prod/pavg
    berrypel = 0.d0
    do ix = 1, kset%ngridk( d2)
      do iy = 1, kset%ngridk( d3)
        phase = atan2( aimag( prod( ix, iy)), dble( prod( ix, iy)))
        !if( phase .lt. 0.d0) phase = phase + twopi
        berrypel = berrypel + phase
      end do
    end do
    berrypel = berrypel/dble( nplane)
    phase = atan2( aimag( pavg), dble( pavg))
    !if( phase .lt. 0.d0) phase = phase + twopi
    berrypel = berrypel + phase
    berrypel = 2.d0*berrypel

    berrypion = 0.d0
    do is = 1, nspecies
      zion = spze( is)
      do ist = 1, spnst( is)
        if( spcore( ist, is)) zion = zion - spocc( ist, is)
      end do
      do ia = 1, natoms( is)
        call r3mv( ainv, atposc( :, ia, is), vr)
        berrypion = berrypion + twopi*zion*vr( d1)
      end do
    end do
    det = cmplx( cos( berrypion), sin( berrypion), 8)
    berrypion = atan2( aimag( det), dble( det))
    if( berrypion .lt. -1.d-6) berrypion = berrypion + twopi
    pel = pel + berrypel/twopi*input%structure%crystal%basevect( :, d1)/omega
    pion = pion + berrypion/twopi*input%structure%crystal%basevect( :, d1)/omega
  end do
  ptot = pel + pion

  ! map result to shortest possible vectors
  call r3mv( ainv, ptot, vr)
  vr = vr*omega
  call r3ws( 1.d-16, input%structure%crystal%basevect, vr, vi)
  call r3mv( input%structure%crystal%basevect, vr, ptot)
  pol(:,1) = ptot/omega
  
  call r3mv( ainv, pel, vr)
  vr = vr*omega
  call r3ws( 1.d-16, input%structure%crystal%basevect, vr, vi)
  call r3mv( input%structure%crystal%basevect, vr, pel)
  pol(:,2) = pel/omega

  call r3mv( ainv, pion, vr)
  vr = vr*omega
  call r3ws( 1.d-16, input%structure%crystal%basevect, vr, vi)
  call r3mv( input%structure%crystal%basevect, vr, pion)
  pol(:,3) = pion/omega

  deallocate( evec1, evec2, s, prod, prodv)
  call delete_k_vectors( kset)
  call delete_G_vectors( Gset)
  call delete_Gk_vectors( Gkset)
  call pwmat_destroy
  return

  contains
    subroutine put( pol)
      use m_getunit
      use mod_misc, only: filext
      real(8), intent( in) :: pol(3,3)
      
      integer :: un

      if( mpiglobal%rank == 0) then
        call getunit( un)
        open( un, file='POLARIZATION'//trim( filext), action='write', form='formatted')
        write( un, '("# macroscopic polarization in cartesian directions")')
        write( un, '("# the results are the shortest possible vectors obtained by subtracting multiple polarization quanta")')
        write( un, '("# total polarization")')
        write( un, '(3g26.16e3)') pol(:,1)
        write( un, '("# electronic polarization")')
        write( un, '(3g26.16e3)') pol(:,2)
        write( un, '("# ionic polarization")')
        write( un, '(3g26.16e3)') pol(:,3)
        close( un)
      end if
      call barrier

      return
    end subroutine

    subroutine get( pol)
      use m_getunit
      use mod_misc, only: filext
      real(8), intent( out) :: pol(3,3)
      
      integer :: un, i, l
      character(256) :: fname, buf
      logical :: exist

      write( fname, '("POLARIZATION")')
      fname = trim( fname)//trim( filext)
      inquire( file=trim( fname), exist=exist)
      if( .not. exist) then
        if( mpiglobal%rank == 0) then
          write(*,*)
          write(*,'("Error (macro_polarization): File ",a," does not exist.")') trim( fname)
        end if
        call terminate
      end if
      call getunit( un)
      open( un, file=trim( fname), status='old', form='formatted')
      l = 0
      do
        read( un, '(a)', iostat=i) buf
        if( i /= 0) exit
        if( (buf( 1:1) == '#') .or. (len( trim( buf)) == 0)) cycle
        l = l + 1
        if( l <= 3) read( buf, *) pol(1,l), pol(2,l), pol(3,l)
      end do
      close( un)

      if( l < 3) then
        if( mpiglobal%rank == 0) then
          write(*,*)
          write(*,'("Error (macro_polarization): File ",a," contains not enough data.")') trim( fname)
        end if
        call terminate
      end if

      return
    end subroutine
end subroutine macro_polarization
