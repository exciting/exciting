!--------------------------------------------!
!     Bare Coulomb potential related data    !
!--------------------------------------------!

module mod_coulomb_potential
#include "asserts.fpp"
    use constants, only: pi, twopi, fourpi, real_zero
    use gw_io, only: write_to_file, read_from_file, build_file_name
    use mod_lattice, only: avec
    use mod_product_basis, only: mbsiz, matsiz
    use modmpi, only: terminate_if_false
    use precision, only: dp, i32, max_length => str_32

    implicit none

    private

    character(len=*), parameter :: basename_barc = 'BARC_Q'
    
    ! The lattice summations matrix      
    complex(dp), allocatable, public :: sgm(:,:,:)
    
    ! The matrix representation of the bare coulomb potential in the mixed basis            
    complex(dp), allocatable, public :: barc(:,:)

    ! full set of the eigenvalues of barcoul matrix
    real(dp), allocatable, public :: barcev(:)
      
    ! full set of eigenvectors of barcoul matrix        
    complex(dp), allocatable, public :: vmat(:,:)
    
    ! use a truncation technique for the Coulomb potential
    logical, public :: vccut
    
    ! spherical integral over the Coulomb singularity
    real(dp), public :: rcut
    
    !> Singularity for 0D, 1D and 2D systems
    real(dp), public, protected :: low_dim_singularity

    public :: delete_coulomb_potential, &
              vcoul_q0_0d, &
              vcoul_q0_1d, &
              vcoul_q0_2d, &
              vcoul_q0_3d, &
              vcoul_0d, &
              vcoul_1d, &
              vcoul_2d, &
              vcoul_3d, &
              vcoul_3d_RIM, &
              calculate_singularities_coeff, &
              calculate_bare_coulomb, &
              calculate_sqrt_bare_coulomb, &
              write_barcev_vmat_to_file, &
              read_barcev_vmat_from_file, &
              index_of_first_element_above_threshold
    
contains
 
    subroutine delete_coulomb_potential()
      implicit none
      if (allocated(vmat)) deallocate(vmat)
      if (allocated(barcev)) deallocate(barcev)
    end subroutine
    
    subroutine vcoul_q0_0d(sing)
        implicit none
        real(dp), intent(out) :: sing
        rcut = 0.5d0*dsqrt(dot_product(avec(:,3),avec(:,3)))
        sing = 2.d0*pi*rcut**2
    end subroutine

    subroutine vcoul_q0_1d(nkpt, sing)
        implicit none
        integer(i32), intent(in)  :: nkpt
        real(dp),    intent(out) :: sing
        real(dp) :: v(3), omega_xy, omega_BZ
        real(dp) :: a, b, c, a2, b2, a2b2
        real(dp) :: t1, t2, t3, rws, beta
        real(dp), parameter :: gamma = -0.5772156649d0 + log(2.d0)
        real(dp), parameter :: small = 1.d-6
        !
        ! check if the unit cell orthorombic
        !
        t1 = dot_product(avec(:,1), avec(:,2))
        if (t1 > small) then
            write(*,*)
            write(*,*) 'Error(mod_coulomb_potential) Unitcell should be orthorombic and parallel to the cartesian vectors.'
            write(*,*) '    Fix the primitive cell geometry.'
            write(*,*)
            stop
        end if
        !
        ! \Omega_xy = ab-plane unit cell area
        !
        call r3cross(avec(:,1), avec(:,2), v)
        omega_xy = sqrt(dot_product(v, v))
        !
        ! 1D BZ volume / Nk
        !
        c = sqrt(dot_product(avec(:,3), avec(:,3)))
        beta =  2.d0*pi / c / dble(nkpt)
        !
        t1 = (gamma - log(0.5d0*beta) + 1.d0) * omega_xy
        !
        ! (1) Approximation: Integral over sphere with the same area as \Omega_xy
        ! rws = sqrt(omega_xy/pi)
        ! t2 = pi * ( rws**2 * log(rws) - 0.5d0 * rws**2 )
        !
        ! (2) Exact integral for rectangular cell: 
        !
        ! \int_{-a/2}^{a/2} \int_{-b/2}^{b/2} ln(\sqrt{x^2+y^2}) dx dy
        !
        a    = 0.5d0 * sqrt(dot_product(avec(:,1), avec(:,1)))
        b    = 0.5d0 * sqrt(dot_product(avec(:,2), avec(:,2)))
        a2   = a**2
        b2   = b**2
        a2b2 = a2+b2
        t2   = -pi*b2 + 2.d0 * ( &
                  2.d0 * b2 * datan(a/b) + &
                  a2b2 * datan(b/a) + &
                  a*b  * (-3.d0 + dlog(a2b2)) &
                  )
        ! Final value
        sing = 2.d0*(t1 - t2)
    end subroutine

    subroutine vcoul_q0_2d(nkpt, sing)
        use incgamma, only: incgam
        implicit none
        integer(i32), intent(in)  :: nkpt
        real(dp),    intent(out) :: sing
        real(dp) :: ab_plane, ab_norm(3), q0_vol
        real(dp), parameter :: eulergamma = 0.5772156649015329
        !--------------------------------------------------------
        ! Spherically averaged value of the integral around q->0
        !--------------------------------------------------------
        ! cutoff length
        rcut = 0.5d0*norm2(avec(:,3))
        ! ab-plane surface area
        call r3cross(avec(:,1), avec(:,2), ab_norm(:))
        ab_plane = norm2(ab_norm(:))
        q0_vol   = twopi / sqrt(pi*ab_plane*nkpt)
        sing     = incgam(0.d0, q0_vol*rcut) + eulergamma + log(q0_vol*rcut)
        sing     = 2.d0 * ab_plane * sing * dble(nkpt)
    end subroutine


    subroutine vcoul_q0_3d(nkpt, sing)
        use modmain, only: omega
        implicit none
        integer(i32), intent(in)  :: nkpt
        real(dp),    intent(out) :: sing
        real(dp) :: omega_BZ, V, beta
        !--------------------------------------------------------
        ! Spherically averaged value of the integral around q->0
        !--------------------------------------------------------
        omega_BZ = (2.d0*pi)**3 / omega
        V        = omega_BZ / dble(nkpt)
        beta     = ( 3.d0*V/(4.d0*pi) )**(1.d0/3.d0)
        sing     = 4.d0*pi/omega_BZ * beta
        ! 4pi/Nk prefactor is already accounted
    end subroutine


    subroutine vcoul_0d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        implicit none
        logical,      intent(in)  :: Gamma
        integer(i32),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(dp),      intent(out) :: vcoul(:)
        integer(i32) :: igk, igk0
        real(dp)    :: k
        if (Gamma) then
            igk0 = 2
            vcoul(1) = 0.d0
        else
            igk0 = 1
        end if
        do igk = igk0, Gkset%ngk(1,ik)
            k = Gkset%gkc(igk,1,ik)
            vcoul(igk) = 4.d0*pi/k**2 * (1.d0 - dcos( k*rcut ))
        end do
    end subroutine


    subroutine vcoul_1d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        use mod_quadrature
        implicit none
        logical,      intent(in)  :: Gamma
        integer(i32),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(dp),      intent(out) :: vcoul(:)
        ! local
        integer(i32) :: igk, igk0, n
        real(dp)    :: a, b, vgpk(3), intf, t1
        real(dp), parameter :: small = 1.d-6

        ! Romberg integration
        integer(i32), parameter :: dim_num = 2
        real(dp)    :: alim(dim_num), blim(dim_num)
        integer(i32) :: sub_num(dim_num)
        integer(i32) :: it_max, ind, eval_num
        real(dp)    :: tol

        ! check if the unit cell orthorombic
        t1 = dot_product(avec(:,1), avec(:,2))
        if (t1 > small) then
            write(*,*)
            write(*,*) 'Error(mod_coulomb_potential) Unitcell should be orthorombic and parallel to the cartesian vectors.'
            write(*,*) '    Fix the primitive cell geometry.'
            write(*,*)
            stop
        end if

        if (Gamma) then
            igk0 = 2
        else
            igk0 = 1
        end if

        a = 0.5d0*sqrt(dot_product(avec(:,1), avec(:,1)))
        b = 0.5d0*sqrt(dot_product(avec(:,2), avec(:,2)))
        
        alim(1) = 0 ; alim(2) = -b
        blim(1) = a ; blim(2) =  b
        
        sub_num(1) = nint(dble(64)*a/10.d0)
        sub_num(2) = nint(dble(64)*b/10.d0)

        it_max = 1000
        tol    = 0.1d0

        do igk = igk0, Gkset%ngk(1,ik)
            vgpk(:) = Gkset%vgkc(:,igk,1,ik)
            ! case q_z -> 0
            if (abs(vgpk(3)) < small) vgpk(3) = small
            call romberg_nd( func, alim, blim, dim_num, sub_num, it_max, tol, intf, ind, eval_num )
            if (ind < 0) print*, 'The error tolerance could not be achieved'
            vcoul(igk) = 2.d0 * 2.d0*intf ! extra factor 2 comes from the limits
        end do

       
    contains

        function func(dim_num, x)
            integer(i32) :: dim_num
            real(dp)    :: func
            real(dp)    :: x(dim_num)
            ! local
            real(dp) :: arg, t1, t2
            real(dp), external :: dbesk0
            t1 = sqrt( x(1)*x(1) + x(2)*x(2) )
            if (t1 > small) then
                arg  = t1 * abs(vgpk(3))
                t2   = dbesk0(arg)
                if (abs(t2) > small) then
                    func = t2 * dcos( vgpk(1)*x(1) + vgpk(2)*x(2) )
                else
                    func = 0.d0
                end if
            else
                func = 0.d0
            end if
            return
        end function

    end subroutine


    real(dp) function K0cosXY(vgpk, x, y)
            implicit none
            real(dp), intent(in) :: vgpk(3)
            real(dp), intent(in) :: x
            real(dp), intent(in) :: y
            ! local
            real(dp) :: arg, k0
            real(dp), external :: dbesk0
            arg = abs(vgpk(3)) * sqrt(x*x+y*y)
            K0cosXY = dbesk0(arg) * cos(vgpk(1)*x + vgpk(2)*y)
    end function


    subroutine vcoul_1d_Rozzi(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        implicit none
        logical,      intent(in)  :: Gamma
        integer(i32),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(dp),      intent(out) :: vcoul(:)
        ! local
        integer(i32) :: igk, igk0
        integer(i32) :: nr, ir
        real(dp)    :: k, kxy, kz, rkxy, rkz, r0
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: fr(:), gr(:), cf(:,:)
        real(dp), parameter :: small = 1.d-6
        real(dp), external :: dbesk0, dbesk1, dbesj0, dbesj1

        ! generate grid
        nr = 128
        allocate(r(nr))
        r0 = 1.d-4
        do ir = 1, nr
            r(ir) = r0 + (dble(ir-1)/dble(nr-1))**3*(rcut-r0)
        end do
        
        allocate(fr(nr), gr(nr), cf(3,nr))
        do igk = 1, Gkset%ngk(1,ik)
            k   = Gkset%gkc(igk,1,ik)
            kxy = sqrt( Gkset%vgkc(1,igk,1,ik)**2 +  Gkset%vgkc(2,igk,1,ik)**2 ) ! k_perpendicular
            kz  = abs(Gkset%vgkc(3,igk,1,ik))
            if ( kz > small ) then
                rkxy = rcut*kxy
                rkz  = rcut*kz
                vcoul(igk) = 4.d0*pi / k**2 * ( &
                             1.d0 + rkxy * dbesj1(rkxy) * dbesk0(rkz) - &
                             rkz * dbesj0(rkxy) * dbesk1(rkz) )
            else if ( (kz < small) .and. (abs(kxy) > small) ) then
                do ir = 1, nr
                    fr(ir) = r(ir)*log(r(ir))*dbesj0(kxy*r(ir))
                end do
                call fderiv(-1, nr, r, fr, gr, cf)
                vcoul(igk) = -4.d0*pi * gr(nr)
            else if ( (kz < small) .and. (abs(kxy) < small)) then
                vcoul(igk) = -pi * rcut**2 * (2.d0*log(rcut)-1.d0)
            end if
        end do
        deallocate(fr, gr, cf)
        deallocate(r)

    end subroutine


    subroutine vcoul_2d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        use modgw, only : Gset, kqset, Gqset, Gqbarc
        implicit none
        logical,      intent(in)  :: Gamma
        integer(i32),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(dp),      intent(out) :: vcoul(:)
        integer(i32) :: igk, igk0
        real(dp)    :: kxy, kz, g_plus_q2, g_plus_q(3)
 
        igk0 = 1
        if (Gamma) then
            igk0 = 2
            vcoul(1) = 0.d0            
        end if

        do igk = igk0, Gqbarc%ngk(1,ik)
            g_plus_q(1:3) = Gset%vgc(1:3, Gqbarc%igkig(igk,1,ik)) + kqset%vqc(1:3,ik)
            g_plus_q2 = dot_product(g_plus_q, g_plus_q)
            kxy = norm2(g_plus_q(1:2))
            kz = g_plus_q(3)
            vcoul(igk) = 4.d0*pi/g_plus_q2 * (1.d0 - exp(-kxy*rcut) * cos(kz*rcut))
        end do

    end subroutine

    
    subroutine vcoul_3d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset, only: Gk_set
        use modgw,         only : Gset, kqset, Gqset, Gqbarc
        !> Indicate if ik is Gamma
        logical,       intent(in)  :: Gamma
        !> k-point index
        integer(i32),  intent(in)  :: ik
        !> G + k vectors
        type(Gk_set),  intent(in)  :: Gkset
        !> 3D bare Coulomb potential
        real(dp),      intent(out) :: vcoul(:)
        !> G + q vector
        real(dp) :: g_plus_q(3)
        integer(i32)   :: igk, igk0

        igk0 = 1
        if (Gamma) then
            vcoul(1) = 0._dp
            igk0 = 2
        endif

        do igk = igk0, Gkset%ngk(1,ik)            
            g_plus_q(1:3) = Gset%vgc(1:3,Gqbarc%igkig(igk,1,ik)) + kqset%vqc(1:3,ik)
            vcoul(igk) =  fourpi / dot_product(g_plus_q, g_plus_q)
        end do

    end subroutine


    subroutine vcoul_3d_RIM(Gamma, ngridk, ik, Gkset, vc)
        use modmain, only: bvec, omega
        use mod_kpointset
        implicit none
        ! input/output
        logical,   intent(in)  :: Gamma
        integer(i32),   intent(in)  :: ngridk(3)
        integer(i32),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(dp),      intent(out) :: vc(Gkset%ngk(1,ik))
        ! local
        integer(i32) :: i, i1, i2, i3, nq, iq
        integer(i32) :: ngk, igk, igk0
        integer(i32) :: n(3), n0
        real(dp)    :: b(3), bmin, bmax, bvol
        real(dp)    :: vgpk(3), intf
        real(dp), parameter :: small = 1.d-6
        ! Romberg integration
        integer(i32), parameter :: dim_num = 3
        real(dp)    :: alim(dim_num), blim(dim_num)
        integer(i32) :: sub_num(dim_num)
        integer(i32) :: it_max, ind, eval_num
        real(dp)    :: tol

        ! Rectangular integration volume
        bvol = (2.d0*pi)**3 / omega / dble(product(ngridk))

        ! Determine the integration volume size
        b(1) = 0.5d0 * sqrt(dot_product(bvec(:,1),bvec(:,1))) / dble(ngridk(1))
        b(2) = 0.5d0 * sqrt(dot_product(bvec(:,2),bvec(:,2))) / dble(ngridk(2))
        b(3) = 0.5d0 * sqrt(dot_product(bvec(:,3),bvec(:,3))) / dble(ngridk(3))
        bmin = minval(b)
        bmax = maxval(b)

        n0 = 2
        if (Gamma) then
            n0 = 4*n0
        end if

        ! create uniform 3-d grid
        do i = 1, 3
            n(i) = nint(dble(n0)*b(i)/bmin)
        end do
    
        ! Integration over a small volume around k-point
        alim(1) = -b(1) ; alim(2) = -b(2) ; alim(3) = -b(3)
        blim(1) =  b(1) ; blim(2) =  b(2) ; blim(3) =  b(3)
        
        sub_num(:) = n(:)

        it_max = 1000
        tol    = 0.1d0

        ngk = Gkset%ngk(1,ik)
        do igk = 1, ngk
            vgpk(:) = Gkset%vgkc(:,igk,1,ik)
            call romberg_nd( func, alim, blim, dim_num, sub_num, it_max, tol, intf, ind, eval_num )
            if (ind < 0) print*, 'The error tolerance could not be achieved'
            vc(igk) = 4.d0*pi * intf / bvol
        end do

    contains

        function func(dim_num, x)
            integer(i32) :: dim_num
            real(dp)    :: func
            real(dp)    :: x(dim_num)
            ! local
            real(dp) :: v(3)
            v(1:3) = x(1:3) + vgpk(1:3)
            func   = 1.d0 / ( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )
            return
        end function

    end subroutine

    !> Compute the the coefficients needed to treat the singularities of the
    !> Coulomb potential and the self-energy
    subroutine calculate_singularities_coeff( cutoff_type, selfenergy_singularity_treatment, &
      & nkpt, coeff_s2_singularity )
      !> Type of Coulomb cutoff used
      character(len=*), intent(in) :: cutoff_type
      !> Treatment of the singularity for the computation of the self-energy
      character(len=*), intent(in) :: selfenergy_singularity_treatment
      !> Number of k/q points in the BZ 
      integer, intent(in) :: nkpt
      !> Coefficient for the integration of the self-energy singularity
      real(dp), intent(inout) :: coeff_s2_singularity
      
      select case ( trim(cutoff_type) )
        case('0d')
          call vcoul_q0_0d( low_dim_singularity )
    
        case('1d')
          call vcoul_q0_1d( nkpt, low_dim_singularity )
    
        case('2d')
          call vcoul_q0_2d( nkpt, low_dim_singularity )
        
        case('none')
          select case ( trim(selfenergy_singularity_treatment) )
            case('mpb')
              ! Auxiliary function method
              call setsingc
            case('crg')
              ! Auxiliary function method
              call calc_q0_singularities
            case('avg')
              ! Spherical average
              call vcoul_q0_3d( nkpt, coeff_s2_singularity )
            case('rim')
              ! Spherical average
            case default
              call calc_q0_singularities
          end select
      end select
    end subroutine

    !> Matrix with the bare Coulomb potential is calculated. 
    !> Then, it is diagonalized
    subroutine calculate_bare_coulomb( iq )
      integer(i32), intent(in) :: iq
      
      ! Get coulomb matrix im MB basis, its eigenvalues and eigenvectors
      call calcbarcmb( iq )
    end subroutine

    !> Take the square root of the matrix with the bare Coulomb potential
    !> Filter the results according to `eigenvalue_tol`
    subroutine calculate_sqrt_bare_coulomb( iq, eigenvalue_tol, remove_g_equal_zero )
      !> Index of the q-point
      integer(i32), intent(in) :: iq
      !> Eigenvalues smaller than `eigenvalue_tol` are discarded
      real(dp), intent(in) :: eigenvalue_tol
      !> If `.true.`, remove the eigenvectors closest to `G=0`
      logical, intent(in) :: remove_g_equal_zero
          
      ! Set v-diagonal MB and reduce its size
      if( (.not. vccut) .and. remove_g_equal_zero ) call setbarcev( real_zero, remove_g_equal_zero )
      call setbarcev( eigenvalue_tol, remove_g_equal_zero )
    end subroutine


    subroutine write_barcev_vmat_to_file( iq, file_format, threshold )
      integer(i32), intent(in) :: iq
      !> Format of the output file
      character(len=*), intent(in) :: file_format
      !> Only write to file the eigenvalues that are >= the threshold (and the corresponding eigenvectors)
      real(dp), intent(in) :: threshold

      character(len=max_length) :: file_name
      integer(i32) :: idx

      CALL_ASSERT( allocated(vmat), 'vmat must be allocated')
      CALL_ASSERT( allocated(barcev), 'barcev must be allocated')
      call build_file_name( basename_barc, iq, file_name )
      idx = index_of_first_element_above_threshold( barcev, threshold )
      call write_to_file( file_name, barcev(idx:), vmat(:, idx:), file_format )

    end subroutine


    subroutine read_barcev_vmat_from_file( iq, file_format )
      integer(i32), intent(in) :: iq
      character(len=*), intent(in) :: file_format

      character(len=max_length) :: file_name

      call build_file_name( basename_barc, iq, file_name )
      call read_from_file( file_name, barcev, vmat, file_format )
      matsiz = size( vmat, 1 )
      mbsiz = size( vmat, 2 )
      call terminate_if_false( size( barcev ) == mbsiz, 'Different number of eigenvalues and eigenvectors in ' // trim( file_name ) )

    end subroutine

    !(private) This can be easily replaced by findloc. But older versions of gfortran do not support it
    pure integer(i32) function index_of_first_element_above_threshold( array, threshold ) result(idx)
      !> Array to find the index of the first element >= threshold. It must be sorted in ascending order.
      real(dp), intent(in) :: array(:)
      real(dp), intent(in) :: threshold

      idx = count( array < threshold ) + 1
    end function


end module
