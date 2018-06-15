!--------------------------------------------!
!     Bare Coulomb potential related data    !
!--------------------------------------------!

module mod_coulomb_potential

    use modmain, only: pi, twopi, fourpi, avec
    implicit none
    
    ! The lattice summations matrix      
    complex(8), allocatable :: sgm(:,:,:)
    
    ! The matrix representation of the bare coulomb potential in the mixed basis            
    complex(8), allocatable :: barc(:,:)

    ! full set of the eigenvalues of barcoul matrix
    real(8), allocatable :: barcev(:)
      
    ! full set of eigenvectors of barcoul matrix        
    complex(8), allocatable :: vmat(:,:)
    
    ! use a truncation technique for the Coulomb potential
    logical :: vccut
    
    ! spherical integral over the Coulomb singularity
    real(8) :: rcut
    
contains
 
    subroutine delete_coulomb_potential()
      implicit none
      if (allocated(vmat)) deallocate(vmat)
      if (allocated(barcev)) deallocate(barcev)
    end subroutine
    
    subroutine vcoul_q0_0d(sing)
        implicit none
        real(8), intent(out) :: sing
        rcut = 0.5d0*dsqrt(dot_product(avec(:,3),avec(:,3)))
        sing = 2.d0*pi*rcut**2
    end subroutine

    subroutine vcoul_q0_1d(nkpt, sing)
        implicit none
        integer(4), intent(in)  :: nkpt
        real(8),    intent(out) :: sing
        real(8) :: v(3), omega_xy, omega_BZ
        real(8) :: a, b, c, a2, b2, a2b2
        real(8) :: t1, t2, t3, rws, beta
        real(8), parameter :: gamma = -0.5772156649d0 + log(2.d0)
        real(8), parameter :: small = 1.d-6
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
        ! 4pi/Nk prefactor is due to definition of the singular term in \Self_x
        sing = sing / (4.d0*pi*dble(nkpt))
        if (.true.) then
            print*, ''
            print*, '1D: LIMIT q->0'
            print*, 'omega_xy   =', omega_xy
            print*, 'c          =', c
            print*, 'beta       =', beta
            print*, 't1         =', t1
            print*, 't2         =', t2
            print*, 'sing       =', sing
            print*, ''
        end if
    end subroutine

    subroutine vcoul_q0_2d(nkpt, sing)
        implicit none
        integer(4), intent(in)  :: nkpt
        real(8),    intent(out) :: sing
        real(8) :: ab_plane, ab_norm(3), q0_vol
        !--------------------------------------------------------
        ! Spherically averaged value of the integral around q->0
        !--------------------------------------------------------
        ! cutoff length
        rcut = 0.5d0*dsqrt(dot_product(avec(:,3),avec(:,3)))
        ! ab-plane surface area
        call r3cross(avec(:,1), avec(:,2), ab_norm(:))
        ab_plane = sqrt(dot_product(ab_norm(:), ab_norm(:)))
        q0_vol   = 2.d0*pi / sqrt(pi*ab_plane*nkpt)
        sing     = q0_vol*rcut - ((q0_vol*rcut)**2.0d0)/4.0d0
        sing     = 2.d0 * ab_plane * sing * dble(nkpt)
        ! 4pi/Nk prefactor is due to definition of the singular term in \Self_x
        sing = sing / (4.d0*pi*dble(nkpt))
        if (.true.) then
            write(*,*)
            write(*,*) '2D: LIMIT q->0'
            write(*,*) ' nqpt     = ', nkpt
            write(*,*) ' rcut     = ', rcut
            write(*,*) ' ab_norm  = ', ab_norm
            write(*,*) ' ab_plane = ', ab_plane
            write(*,*) ' q0_vol   = ', q0_vol
            write(*,*) ' sing     = ', sing
            write(*,*)
        end if
    end subroutine


    subroutine vcoul_q0_3d(nkpt, sing)
        use modmain, only: omega
        implicit none
        integer(4), intent(in)  :: nkpt
        real(8),    intent(out) :: sing
        real(8) :: omega_BZ, V, beta
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
        integer(4),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(8),      intent(out) :: vcoul(:)
        integer(4) :: igk, igk0
        real(8)    :: k
        if (Gamma) then
            igk0 = 2
            vcoul(1) = 0.d0
        else
            igk0 = 1
        end if
        do igk = igk0, Gkset%ngk(1,ik)
            k = Gkset%gkc(igk,1,ik)
            vcoul(igk) = 4.d0*pi/k**2 * (1.d0 - cos( k*rcut ))
        end do
    end subroutine


    subroutine vcoul_1d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        use mod_quadrature
        implicit none
        logical,      intent(in)  :: Gamma
        integer(4),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(8),      intent(out) :: vcoul(:)
        ! local
        integer(4) :: igk, igk0
        integer(4) :: nx, ny, i, j
        real(8)    :: a, b, gpk(3), t1
        real(8), allocatable :: x(:), wx(:)
        real(8), allocatable :: y(:), wy(:)
        real(8), parameter :: small = 1.d-6
        ! integration
        integer(4) :: qopt, ntrial, npts, ierr
        real(8)    :: accuracy
        real(8) :: xx_, yy_, quad

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

        if (.true.) then

            !------------------------------------
            ! Double Gauss-Legendre quadrature
            !------------------------------------
            nx = 64
            ny = 64
            allocate(x(nx), wx(nx))
            allocate(y(ny), wy(ny))

            ! generate grid
            call gauleg(-a, a, x, wx, nx)
            call gauleg(-b, b, y, wy, ny)

            do igk = igk0, Gkset%ngk(1,ik)
                gpk(:) = Gkset%vgkc(:,igk,1,ik)
                ! case q_z -> 0
                if (abs(gpk(3)) < small) gpk(3) = small
                vcoul(igk) = 0.d0
                do i = 1, nx
                do j = 1, ny
                    vcoul(igk) = vcoul(igk) + &
                                 wx(i) * wy(j) * K0cosXY( x(i), y(j) )
                end do
                end do
                vcoul(igk) = 2.d0 * vcoul(igk)
            end do
            deallocate(x, wx)
            deallocate(y, wy)

        else
            
            ! ===================================================
            ! === Setup for the quadrature of matrix elements ===
            ! ===================================================
            qopt     = 6     ! Quadrature method, see quadrature routine.
            ntrial   = 30    ! Max number of attempts.
            accuracy = 0.001 ! Fractional accuracy required.
            npts     = 6     ! Initial number of point (only for Gauss-Legendre method).
 
            do igk = igk0, Gkset%ngk(1,ik)
                gpk(:) = Gkset%vgkc(:,igk,1,ik)
                if (gpk(3) < small) gpk(3) = small
                call quadrature(K0cos_dy, 0.d0, +a, qopt, quad, ierr, ntrial, accuracy, npts)
                if (ierr /= 0) then
                    write(*,*)
                    write(*,'(a,i3)') 'Error(mod_coulomb_potential::vcoul1d) Accuracy not reached'
                    write(*,*)
                    stop
                end if
                vcoul(igk) = 2.d0*(2.d0*quad)
            end do

        end if

    contains

        real(8) function K0cosXY(x, y)
            implicit none
            real(8), intent(in) :: x
            real(8), intent(in) :: y
            ! local
            real(8) :: arg, k0
            arg = abs(gpk(3)) * sqrt(x*x+y*y)
            call calck0(arg, k0, 1)
            K0cosXY = k0 * cos(gpk(1)*x + gpk(2)*y)
        end function

        real(8) function K0cos(yy)
            real(8), intent(in) :: yy
            real(8) :: k0, rho, arg
            ! K0cos(y) = K0(\rho*|qpg_z|)*COS(x.qpg_x+y*qpg_y)
            rho = sqrt(xx_**2+yy**2) 
            arg = abs(gpk(3)) * rho
            call CALCK0(arg, k0, 1)
            K0cos = k0 * cos(gpk(1)*xx_ + gpk(2)*yy)
        end function K0cos

        real(8) function K0cos_dy(xx)
            real(8), intent(in) :: xx
            integer(4) :: ierr
            real(8)    :: quad
            ! K0cos_dy(x)=\int_{-b/2}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
            xx_ = xx ! make it visible in K0cos
            call quadrature(K0cos, -b, +b, qopt, quad, ierr, ntrial, accuracy, npts)
            if (ierr /= 0) then
                write(*,*)
                write(*,'(a,i3)') 'Error(mod_coulomb_potential::K0cos_dy) Accuracy not reached'
                write(*,*)
                stop
            end if
            K0cos_dy = quad
        end function K0cos_dy

    end subroutine


    subroutine vcoul_2d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        implicit none
        logical,      intent(in)  :: Gamma
        integer(4),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(8),      intent(out) :: vcoul(:)
        integer(4) :: igk, igk0
        real(8)    :: vkc(3), k, kxy, kz
        if (Gamma) then
            igk0 = 2
            vcoul(1) = 0.d0
        else
            igk0 = 1
        end if
        do igk = igk0, Gkset%ngk(1,ik)
            k      = Gkset%gkc(igk,1,ik)
            vkc(:) = Gkset%vgkc(:,igk,1,ik)
            kxy    = sqrt(vkc(1)*vkc(1)+vkc(2)*vkc(2))
            kz     = abs(vkc(3))
            vcoul(igk) = 4.d0*pi/k**2 * (1.d0 - exp(-kxy*rcut) * cos(kz*rcut))
        end do
    end subroutine


    subroutine vcoul_3d(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        implicit none
        logical,      intent(in)  :: Gamma
        integer(4),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(8),      intent(out) :: vcoul(:)
        integer(4) :: igk, igk0
        real(8)    :: k
        if (Gamma) then
            igk0 = 2
            vcoul(1) = 0.d0
        else
            igk0 = 1
        end if
        do igk = igk0, Gkset%ngk(1,ik)
            k = Gkset%gkc(igk,1,ik)
            vcoul(igk) = 4.d0*pi / k**2
        end do
    end subroutine


    subroutine vcoul_3d_RIM(Gamma, ngridk, ik, Gkset, vc)
        use modmain, only: bvec, omega, pi
        use mod_kpointset
        implicit none
        ! input/output
        logical(4),   intent(in)  :: Gamma
        integer(4),   intent(in)  :: ngridk(3)
        integer(4),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(8),      intent(out) :: vc(Gkset%ngk(1,ik))
        ! local
        integer(4) :: i, i1, i2, i3, nq, iq
        integer(4) :: ngk, igk, igk0
        integer(4) :: n(3), n0
        real(8)    :: b(3), bmin, bmax, bvol
        real(8)    :: q, vgpq(3), gpq2, intf
        real(8), allocatable :: q1(:), q2(:), q3(:), vq(:,:)
        real(8), allocatable :: w1(:), w2(:), w3(:), wq(:)
        real(8), parameter :: small = 1.d-6

        ! Rectangular integration volume
        bvol = (2.d0*pi)**3 / omega / dble(product(ngridk))

        ! Determine the integration volume size
        b(1) = 0.5d0 * sqrt(dot_product(bvec(:,1),bvec(:,1))) / dble(ngridk(1))
        b(2) = 0.5d0 * sqrt(dot_product(bvec(:,2),bvec(:,2))) / dble(ngridk(2))
        b(3) = 0.5d0 * sqrt(dot_product(bvec(:,3),bvec(:,3))) / dble(ngridk(3))
        bmin = minval(b)
        bmax = maxval(b)

        n0 = 3
        if (Gamma) then
            n0 = 4*n0
        end if

        ! create uniform 3-d grid
        do i = 1, 3
            n(i) = dble(n0)*nint(b(i)/bmin)
        end do
        
        ! print*, 'grid=', n
        ! print*, 'bmin=', bmin
        ! print*, 'b=', b

        ! Generate the integration grid and corresponding weights
        allocate(q1(n(1)), w1(n(1)))
        allocate(q2(n(2)), w2(n(2)))
        allocate(q3(n(3)), w3(n(3)))
        call gauleg(-b(1), b(1), q1, w1, n(1))
        call gauleg(-b(2), b(2), q2, w2, n(2))
        call gauleg(-b(3), b(3), q3, w3, n(3))

        ! setup q-points
        nq = product(n)
        allocate(vq(3,nq), wq(nq))
        iq = 0
        do i1 = 1, n(1)
        do i2 = 1, n(2)
        do i3 = 1, n(3)
            iq = iq+1
            vq(1,iq) = q1(i1)
            vq(2,iq) = q2(i2)
            vq(3,iq) = q3(i3)
            wq(iq)   = w1(i1) * w2(i2) * w3(i3)
        end do
        end do
        end do
        deallocate(q1, q2, q3)
        deallocate(w1, w2, w3)

        ! test
        ! intf = 0.d0
        ! do iq = 1, nq
        !     intf = intf + wq(iq)
        ! end do
        ! print*, 'TEST integral=', intf, bvol

        ! Integration over a small volume around k-point
        ngk = Gkset%ngk(1,ik)
        do igk = 1, ngk
            intf = 0.d0
            do iq = 1, nq
                vgpq(:) = vq(:,iq) + Gkset%vgkc(:,igk,1,ik)
                gpq2    = dot_product(vgpq,vgpq)
                intf    = intf + wq(iq) / gpq2
            end do
            vc(igk) = 4.d0*pi * intf / bvol
        end do
        deallocate(vq, wq)

    end subroutine


end module
