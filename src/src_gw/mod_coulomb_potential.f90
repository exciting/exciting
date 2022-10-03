!--------------------------------------------!
!     Bare Coulomb potential related data    !
!--------------------------------------------!

module mod_coulomb_potential
    use precision, only: dp, i32
    use constants, only: pi, twopi, fourpi
    use modmain, only: avec
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
         ! 4pi/Nk prefactor is due to definition of the singular term in \Self_x
        sing = sing / (4.d0*pi)
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
            vcoul(igk) = 4.d0*pi/k**2 * (1.d0 - dcos( k*rcut ))
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
        integer(4) :: igk, igk0, n
        real(8)    :: a, b, vgpk(3), intf, t1
        real(8), parameter :: small = 1.d-6

        ! Romberg integration
        integer(4), parameter :: dim_num = 2
        real(8)    :: alim(dim_num), blim(dim_num)
        integer(4) :: sub_num(dim_num)
        integer(4) :: it_max, ind, eval_num
        real(8)    :: tol

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
            integer(4) :: dim_num
            real(8)    :: func
            real(8)    :: x(dim_num)
            ! local
            real(8) :: arg, t1, t2
            real(8), external :: dbesk0
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


    real(8) function K0cosXY(vgpk, x, y)
            implicit none
            real(8), intent(in) :: vgpk(3)
            real(8), intent(in) :: x
            real(8), intent(in) :: y
            ! local
            real(8) :: arg, k0
            real(8), external :: dbesk0
            arg = abs(vgpk(3)) * sqrt(x*x+y*y)
            K0cosXY = dbesk0(arg) * cos(vgpk(1)*x + vgpk(2)*y)
    end function


    subroutine vcoul_1d_Rozzi(Gamma, ik, Gkset, vcoul)
        use mod_kpointset
        implicit none
        logical,      intent(in)  :: Gamma
        integer(4),   intent(in)  :: ik
        type(Gk_set), intent(in)  :: Gkset
        real(8),      intent(out) :: vcoul(:)
        ! local
        integer(4) :: igk, igk0
        integer(4) :: nr, ir
        real(8)    :: k, kxy, kz, rkxy, rkz, r0
        real(8), allocatable :: r(:)
        real(8), allocatable :: fr(:), gr(:), cf(:,:)
        real(8), parameter :: small = 1.d-6
        real(8), external :: dbesk0, dbesk1, dbesj0, dbesj1

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
        !> G + q vector (I assume)
        real(dp) :: gpq(3)
        integer(i32)   :: igk

        do igk = 1, Gkset%ngk(1,ik)            
            gpq(1:3) = Gset%vgc(1:3,Gqbarc%igkig(igk,1,ik)) + kqset%vqc(1:3,ik)
            vcoul(igk) =  fourpi / dot_product(gpq, gpq)
        end do
        if (Gamma) vcoul(1) = 0._dp

    end subroutine


    subroutine vcoul_3d_RIM(Gamma, ngridk, ik, Gkset, vc)
        use modmain, only: bvec, omega
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
        real(8)    :: vgpk(3), intf
        real(8), parameter :: small = 1.d-6
        ! Romberg integration
        integer(4), parameter :: dim_num = 3
        real(8)    :: alim(dim_num), blim(dim_num)
        integer(4) :: sub_num(dim_num)
        integer(4) :: it_max, ind, eval_num
        real(8)    :: tol

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
        
        print*, 'grid=', n
        print*, 'bmin=', bmin
        print*, 'b=', b

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
            integer(4) :: dim_num
            real(8)    :: func
            real(8)    :: x(dim_num)
            ! local
            real(8) :: v(3)
            v(1:3) = x(1:3) + vgpk(1:3)
            func   = 1.d0 / ( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )
            return
        end function

    end subroutine

end module
