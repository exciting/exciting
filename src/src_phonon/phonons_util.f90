module phonons_util
  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: terminate_if_false

  implicit none

  contains

    !> Impose the acoustic sum rule on a set of dynamical matrices \({\bf D}({\bf p})\).
    !>
    !> On exit, the dynamical matrices are scaled such that the three smallest
    !> eigenvalues \({\bf p}=\Gamma\) are zero. \(\Gamma\) must be part of point set.
    subroutine ph_util_sumrule( pset, dynp )
      use mod_kpointset, only: k_set
      use math_utils, only: is_square
      use m_linalg, only: zhediag
      !> set of reciprocal space points \({\bf p}\)
      type(k_set), intent(in) :: pset
      !> dynamical matrices \({\bf D}({\bf p})\)
      complex(dp), intent(inout) :: dynp(:,:,:)

      integer :: ip, ip0, n, i, j, k
      
      real(dp), allocatable :: eval(:)
      complex(dp), allocatable :: dyn0(:,:), evec(:,:)

      call assert( is_square( dynp(:,:,1) ), 'Dynamical matrices are not square.' )
      n = size( dynp, dim=1 )

      allocate( eval(n), dyn0(n,n), evec(n,n) )

      call findkptinset( [0.d0, 0.d0, 0.d0], pset, i, ip0 )
      dyn0 = 0.5d0 * (dynp(:, :, ip0) + conjg( transpose( dynp(:, :, ip0) ) ))
      call zhediag( dyn0, eval, evec=evec )

      do ip = 1, pset%nkpt
        do i = 1, n
          do j = 1, n
            do k = 1, 3
              dynp(i, j, ip) = dynp(i, j, ip) - eval(k) * evec(i, k) * conjg( evec(j, k) )
            end do
          end do
        end do
      end do

      deallocate( eval, dyn0, evec )
    end subroutine ph_util_sumrule

    !> Symmetrize a dynamical matrix \({\bf D}({\bf p}\) with a given
    !> set of crystal symmetries.
    !> See also [[ph_util_symapp_dyn(subroutine)]].
    subroutine ph_util_symmetrize_dyn( vpl, dyn, nsym, isym )
      use constants, only: zzero
      use math_utils, only: is_square
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> dynamical matrix \({\bf D}({\bf p})\)
      complex(dp), intent(inout) :: dyn(:,:)
      !> number of crystal symmetries
      integer, intent(in) :: nsym
      !> indices of crystal symmetries in global symmetry arrays
      integer, intent(in) :: isym(:)

      integer :: n, jsym
      complex(dp), allocatable :: tmp(:,:)

      call assert( is_square(dyn), 'Dynamical matrix is not square.' )
      n = size( dyn, dim=1 )

      ! make Hermitian
      dyn = (dyn + conjg( transpose( dyn ) )) / 2

      ! symmetrize
      allocate( tmp(n, n), source=zzero )
      do jsym = 1, nsym
        call ph_util_symapp_dyn( isym(jsym), vpl, dyn, tmp, dir=1, matrix='T' )
      end do
      tmp = tmp / nsym

      dyn = zzero
      do jsym = 1, nsym
        call ph_util_symapp_dyn( isym(jsym), vpl, tmp, dyn, dir=-1, matrix='T' )
      end do
      dyn = dyn / nsym

      deallocate( tmp )
    end subroutine ph_util_symmetrize_dyn

    !> Apply a symmetry operation to a dynamical matrix \({\bf D}({\bf p})\).
    !>
    !> This is done by applying the sandwich product
    !> \[ {\bf S}({\bf p})\, {\bf D}({\bf p})\, {\bf S}^\dagger({\bf p}) \;, \]
    !> where \({\bf S}({\bf p})\) is either the symmetry matrix \({\bf \Gamma}({\bf p})\) or
    !> \({\bf T}({\bf p})\) according to *Maradudin, Vosko, Rev. Mod. Phys. **40**, 1 (1968)*
    !> if `dir >= 0` or the respective Hermitian conjugate if `dir < 0`.
    !> The result is added to the input argument `dyns`.
    !> See also [[ph_util_symmetry_G(subroutine)]] and [[ph_util_symmetry_T(subroutine)]].
    subroutine ph_util_symapp_dyn( isym, vpl, dyn, dyns, &
        dir, matrix )
      use constants, only: zzero, zone
      use mod_atoms, only: natmtot
      !> index of symmetry operation in global symmetry arrays
      integer, intent(in) :: isym
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> dynamical matrix \({\bf D}({\bf p})\)
      complex(dp), intent(in) :: dyn(:,:)
      !> symmetrized dynamical matrix
      complex(dp), intent(inout) :: dyns(:,:)
      !> direction of symmetry application (default: `1`)
      integer, optional, intent(in) :: dir
      !> type of symmetry matrix 
      !>(`G` for \({\bf \Gamma}\), `T` for \({\bf T}\); default: `G`)
      character, optional, intent(in) :: matrix

      integer :: direction
      character :: m, c1, c2
      complex(dp), allocatable :: aux(:,:), S(:,:)

      direction = 1
      if( present( dir ) ) direction = dir
      if( present( matrix ) ) m = matrix
      if( m /= 'G' .or. m /= 'T') m = 'G'

      allocate( aux(3*natmtot, 3*natmtot) )

      if( m == 'T') then
        S = ph_util_symmetry_T( isym, vpl )
      else
        S = ph_util_symmetry_G( isym, vpl )
      end if
      if( direction >= 0 ) then
        c1 = 'n'
        c2 = 'c'
      else
        c1 = 'c'
        c2 = 'n'
      end if
      call zgemm( c1, 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             S, 3*natmtot, &
             dyn, size( dyn, dim=1 ), zzero, &
             aux, 3*natmtot )
      call zgemm( 'n', c2, 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             aux, 3*natmtot, &
             S, 3*natmtot, zone, &
             dyns, size( dyns, dim=1 ) )

      deallocate( S, aux)
    end subroutine ph_util_symapp_dyn

    !> Get the symmetry matrix \({\bf \Gamma}({\bf p})\) for a given
    !> symmetry opertation
    !>
    !> according to *Maradudin, Vosko, Rev. Mod. Phys. **40**, 1 (1968)*.
    !> Note: In exciting, the translation is applied before the rotation!
    function ph_util_symmetry_G( isym, vpl ) result( G )
      use constants, only: zzero, zone, twopi
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use math_utils, only: all_zero
      use mod_symmetry, only: lsplsymc, isymlat, symlat, symlatc, ieqatom, vtlsymc
      use modinput
      !> index of symmetry operation in global symmetry arrays
      integer, intent(in) :: isym
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> symemtry matrix
      complex(dp), allocatable :: G(:,:)

      integer :: lspl, ilspl, is, ia, ias, ja, jas, i, j
      real(dp) :: sl(3, 3), sc(3, 3), v1(3), v2(3), phase
      complex(dp) :: aux(3*natmtot, 3*natmtot)
      character(256) :: errmsg

      if( allocated( G ) ) deallocate( G )
      allocate( G(3*natmtot, 3*natmtot), source=zzero )

      lspl = lsplsymc(isym)
      ilspl = isymlat(lspl)
      sl = dble( symlat(:, :, ilspl) )
      sc = symlatc(:, :, lspl)
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          i = (ias - 1) * 3 + 1
          ja = ieqatom(ia, is, isym)
          jas = idxas(ja, is)
          j = (jas - 1) * 3 + 1
          v1 = input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord + vtlsymc(:, isym)
          call r3mv( sl, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord, v2 )
          phase = twopi * dot_product( vpl, v2 - v1 )
          G(i:i+2, j:j+2) = cmplx( cos( phase ), sin( phase ), dp ) * sc
        end do
      end do
      ! check if Gamma is unitary
      call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             G, 3*natmtot, &
             G, 3*natmtot, zzero, &
             aux, 3*natmtot )
      do ias = 1, 3*natmtot
        aux(ias, ias) = aux(ias, ias) - zone
      end do
      write( errmsg, '("(ph_util_symmetry_G) Non-unitary Gamma for symmetry ",i2," and q = ",3f13.6)' ) isym, vpl
      call terminate_if_false( all_zero( aux, tol=input%structure%epslat ), trim (errmsg ) )
    end function ph_util_symmetry_G

    !> Get the symmetry matrix \({\bf T}({\bf p})\) for a given
    !> symmetry opertation
    !>
    !> according to *Maradudin, Vosko, Rev. Mod. Phys. **40**, 1 (1968)*.
    !> Note: In exciting, the translation is applied before the rotation!
    function ph_util_symmetry_T( isym, vpl ) result( T )
      use constants, only: zzero, zone, twopi
      use math_utils, only: all_zero
      use mod_symmetry, only: lsplsymc, symlat, vtlsymc
      use mod_atoms, only: natmtot
      use modinput
      !> index of symmetry operation in global symmetry arrays
      integer, intent(in) :: isym
      !> reciprocal space point \({\bf p}\)
      real(dp), intent(in) :: vpl(3)
      !> symemtry matrix
      complex(dp), allocatable :: T(:,:)

      integer :: lspl, ias
      real(dp) :: sl(3, 3), v1(3), phase
      complex(dp) :: aux(3*natmtot, 3*natmtot)
      character(256) :: errmsg

      T = ph_util_symmetry_G( isym, vpl )
      lspl = lsplsymc(isym)
      sl = dble( symlat(:, :, lspl) )
      call r3mv( sl, vtlsymc(:, isym), v1 )
      phase = twopi * dot_product( vpl, v1 )
      T = cmplx( cos( phase ), sin( phase ), dp ) * T
      ! check if T is unitary
      call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             T, 3*natmtot, &
             T, 3*natmtot, zzero, &
             aux, 3*natmtot )
      do ias = 1, 3*natmtot
        aux(ias, ias) = aux(ias, ias) - zone
      end do
      write( errmsg, '("(ph_util_symmetry_T): Non-unitary T for symmetry ",i2," and q = ",3f13.6)' ) isym, vpl
      call terminate_if_false( all_zero( aux, tol=input%structure%epslat ), trim( errmsg ))
    end function ph_util_symmetry_T

    !> Diagonalize dynamical matrix and obtain phonon frequencies
    !> and eigenvectors.
    !>
    !> Given a dynamical matrix \(D_{\alpha i, \beta j}\), this calculates
    !> the eigenvalues \(\omega_\nu^2\) and eigenvectors \(e_{\alpha i,\nu}\)
    !> such that 
    !> \[ {\bf M}\, {\bf D}\, {\bf M}
    !>    = {\bf e}\, \operatorname{diag}(\omega_\nu^2)\, {\bf e}^\dagger \;, \]
    !> where \({\bf M} = {\bf M}^0\) if `basis_transform` is not specified
    !> \({\bf M} = {\bf B}^\dagger \, {\bf M}^0\, {\bf B}\) otherwise. Here,
    !> \(M^0_{\alpha i,\beta j} = 1/\sqrt{M_\alpha}\, \delta_{\alpha \beta}\, \delta_{ij}\)
    !> and \(B_{\alpha i, n}\) describes a basis transform from the Cartesian
    !> \((\alpha,i)\) basis into another basis \(n\).
    !>
    !> If `normalized_evec = .false.`, \({\bf u} = {\bf M}\, {\bf e}\) is returned
    !> as the eigenvectors.
    subroutine ph_util_diag_dynmat( dyn, w, evec, &
        basis_transform, normalized_evec )
      use constants, only: zone, zzero
      use mod_atoms, only: natmtot, nspecies, natoms, spmass
      use m_linalg, only: zhediag, zhegauge
      !> dynamical matrix \({\bf D}\)
      complex(dp), intent(in) :: dyn(:,:)
      !> phonon freuqncies \(\omega_\nu\)
      real(dp), intent(out) :: w(:)
      !> phonon eigenvectors \({\bf e}_\nu\)
      complex(dp), intent(out) :: evec(:,:)
      !> unitary basis transformation (default: none)
      complex(dp), optional, intent(in) :: basis_transform(:,:)
      !> return normalized eigenvectors (default: `.true.`)
      logical, optional, intent(in) :: normalized_evec

      integer :: i, is, ia, ip
      complex(dp) :: M(3*natmtot, 3*natmtot), tmp(3*natmtot, 3*natmtot)
      logical :: norm

      norm = .true.
      if( present( normalized_evec ) ) norm = normalized_evec

      ! set up mass matrix
      if( .not. present( basis_transform ) ) then
        M = zzero
        i = 0
        do is = 1, nspecies
          if( spmass(is) <= 0.0_dp ) cycle
          do ia = 1, natoms(is)
            do ip = 1, 3
              i = i + 1
              M(i, i) = 1.0_dp / sqrt( spmass(is) )
            end do
          end do
        end do
      else
        tmp = zzero
        i = 0
        do is = 1, nspecies
          if( spmass(is) <= 0.0_dp ) cycle
          do ia = 1, natoms(is)
            do ip = 1, 3
              i = i + 1
              tmp(i, :) = basis_transform(i, :) / sqrt( spmass(is) )
            end do
          end do
        end do
        call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
               basis_transform, 3*natmtot, &
               tmp, 3*natmtot, zzero, &
               M, 3*natmtot )
      end if
      ! set up dynamical matrix (add mass term)
      tmp = 0.5_dp * (dyn + conjg( transpose( dyn ) ))
      call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             tmp, 3*natmtot, &
             M, 3*natmtot, zzero, &
             evec, 3*natmtot )
      call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
             M, 3*natmtot, &
             evec, 3*natmtot, zzero, &
             tmp, 3*natmtot )
      ! diagonalize
      call zhediag( tmp, w, evec )
      ! assign phonon frequencies
      do i = 1, 3*natmtot
        if( w(i) .ge. 0.0_dp ) then
          w(i) = sqrt( w(i))
        else
          w(i) = -sqrt( abs( w(i) ) )
        end if
      end do
      ! fix unique gauge of degenerate eigenvectors
      call zhegauge( w, evec, 1e-12_dp )

      if( .not. norm ) then
        tmp = evec
        call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
               M, 3*natmtot, &
               tmp, 3*natmtot, zzero, &
               evec, 3*natmtot )
      end if
    end subroutine ph_util_diag_dynmat

    !> write dynamical matrix row to file
    subroutine ph_util_write_dynmat( dyn, fxt, success, &
        directory )
      use m_getunit
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      !> dynamical matrix row
      complex(dp), intent(in) :: dyn(3, natmtot)
      !> file extension
      character(*), intent(in) :: fxt
      !> `.true.` if writing was successful
      logical, intent(out) :: success
      !> path to directory (default: current directory)
      character(*), optional, intent(in) :: directory

      integer :: stat, un, is, ia, ias, ip
      real(dp) :: a, b
      character(256) :: dirname

      dirname = '.'
      if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
      dirname = trim( dirname )//'/'

      call getunit( un )
      open( un, file=trim( dirname )//'DYN_'//trim( fxt )//'.OUT', action='write', form='formatted', iostat=stat )
      success = (stat == 0)
      if( .not. success ) return
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ip = 1, 3
            a = dble( dyn(ip, ias) )
            b = aimag( dyn(ip, ias) )
            if( abs(a) < 1e-12_dp ) a = 0.0_dp
            if( abs(b) < 1e-12_dp ) b = 0.0_dp
            write( un, '(2g18.10," : is = ",i4,", ia = ",i4,", ip = ",i4)', iostat=stat ) a, b, is, ia, ip
            success = success .and. (stat == 0)
          end do
        end do
      end do
      close( un )
    end subroutine ph_util_write_dynmat

    !> read dynamical matrix row from file
    subroutine ph_util_read_dynmat( dyn, fxt, success, &
        directory )
      use m_getunit
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      !> dynamical matrix row
      complex(dp), intent(out) :: dyn(3, natmtot)
      !> file extension
      character(*), intent(in) :: fxt
      !> `.true.` if reading was successful
      logical, intent(out) :: success
      !> path to directory (default: current directory)
      character(*), optional, intent(in) :: directory

      integer :: stat, un, is, ia, ias, ip, js, ja, jp
      real(dp) :: a, b
      character(256) :: dirname

      dirname = '.'
      if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
      dirname = trim( dirname )//'/'

      call getunit( un )
      open( un, file=trim( dirname )//'DYN_'//trim( fxt )//'.OUT', action='read', status='old', form='formatted', iostat=stat )
      success = (stat == 0)
      if( .not. success ) return
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ip = 1, 3
            js = 0; ja = 0; jp = 0;
            read( un, '(2g18.10,tr1,3(tr7,i4))', iostat=stat ) a, b, js, ja, jp
            success = success .and. (stat == 0)
            if( js /= is .or. ja /= ia .or. jp /= ip ) then
              write( *, * )
              write( *, '("Error (ph_io_read_dynmat): Incompatible file content in file")' )
              write( *, '(a)' ) trim(dirname)//'DYN_'//trim(fxt)//'.OUT'
              write( *, '("expected:       is = ",i4,", ia = ",i4,", ip = ",i4)', iostat=stat ) is, ia, ip
              write( *, '("read from file: is = ",i4,", ia = ",i4,", ip = ",i4)', iostat=stat ) js, ja, jp
              success = .false.
              return
            end if
            dyn(ip, ias) = cmplx( a, b, dp )
          end do
        end do
      end do
      close( un )
    end subroutine ph_util_read_dynmat

end module phonons_util
