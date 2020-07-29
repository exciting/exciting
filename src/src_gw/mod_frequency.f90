MODULE mod_frequency

    implicit none

!-------------------------------------------------------------------------------
    type frequency
        character(40) :: fgrid           ! grid type
        character(40) :: fconv           ! real/imaginary frequency
        integer(4)    :: nomeg           ! grid size
        real(8)       :: freqmin         ! lower cutoff frequency
        real(8)       :: freqmax         ! upper cutoff frequency
        real(8), allocatable :: freqs(:) ! frequency grid
        real(8), allocatable :: womeg(:) ! integration weights
    end type frequency

    external gaulag
    external gauleg

CONTAINS

!-------------------------------------------------------------------------------
    subroutine delete_freqgrid(self)
        type(frequency), intent(INOUT) :: self
        if (allocated(self%freqs)) deallocate(self%freqs)
        if (allocated(self%womeg)) deallocate(self%womeg)
    end subroutine


    !-------------------------------------------------------------------------------
    function copy_freqgrid(self_in) result(self_out)
        implicit none
        type(frequency), intent(IN) :: self_in
        type(frequency)             :: self_out
        self_out%fgrid = self_in%fgrid
        self_out%fconv = self_in%fconv
        self_out%nomeg = self_in%nomeg
        self_out%freqmin = self_in%freqmin
        self_out%freqmax = self_in%freqmax
        if (allocated(self_out%freqs)) deallocate(self_out%freqs)
        allocate(self_out%freqs(self_out%nomeg))
        self_out%freqs =  self_in%freqs
        if (allocated(self_out%womeg)) deallocate(self_out%womeg)
        allocate(self_out%womeg(self_out%nomeg))
        self_out%womeg =  self_in%womeg
    end function


!-------------------------------------------------------------------------------
    subroutine generate_freqgrid(self,fgrid,fconv,nomeg,freqmin,freqmax)
        implicit none
        type(frequency), intent(OUT) :: self
        character*(*),   intent(IN)  :: fgrid
        character*(*),   intent(IN)  :: fconv
        integer,         intent(IN)  :: nomeg
        real(8),         intent(IN)  :: freqmin
        real(8),         intent(IN)  :: freqmax
        ! local variables
        integer :: i, ii, j, n, m
        real(8) :: t1, t2, step, sigma, delta, stretch, one
        real(8) :: t, ell, fmax
        real(8), allocatable :: u(:), x(:)
        real(8), allocatable :: wu(:), w1(:), w2(:)
        real(8), parameter   :: pi = 3.1415926535897932385d0

        ! grid type
        self%fgrid = trim(fgrid)

        ! frequency treatment
        self%fconv = trim(fconv)

        ! number of frequencies
        self%nomeg = nomeg

        ! cutoff frequencies
        self%freqmin = freqmin
        self%freqmax = freqmax
        if ( self%freqmin > self%freqmax ) then
            stop 'Error(mod_frequency::generate_freqgrid) freqmin > freqmax: Check your input!'
        end if

        if (allocated(self%freqs)) deallocate(self%freqs)
        if (allocated(self%womeg)) deallocate(self%womeg)

        ! generate frequency integration grid
        select case (self%fgrid)

        case('eqdist') ! Equaly spaced mesh (for tests purposes only)
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            do i = 1, self%nomeg
                self%freqs(i) = self%freqmin + dble(i-1)*(self%freqmax-self%freqmin)/dble(self%nomeg-1)
                self%womeg(i) = 1.d0/dble(self%nomeg)
            enddo

        case('gaulag') ! Grid for Gauss-Laguerre quadrature
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            allocate(wu(self%nomeg))
            call gaulag(self%freqs, wu, self%nomeg, -1.0d0)
            do i = 1, self%nomeg
                self%womeg(i) = wu(i)*dexp(self%freqs(i))*self%freqs(i)
            enddo
            deallocate(wu)

        case('gauleg') ! Grid for Gauss-Legendre quadrature from 0 to infinity
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            allocate( u(n), wu(n) )
            call gauleg(0.d0, 1.d0, u, wu, n)
            do i = 1, n
                self%freqs(i) = u(i) / ( 1.d0 - u(i) )
                self%womeg(i) = wu(i) * ( 1.d0 - u(i) )**(-2)
            end do
            deallocate(u, wu)

        case('gauleg2') ! Grid for Gauss-Legendre quadrature from 0 to freqmax
                                ! +
                                ! Gauss-Legendre quadrature from freqmax to infinity
            self%nomeg = nomeg - mod(nomeg,2)
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))

            n = self%nomeg/2
            allocate(u(n), wu(n))
            call gauleg(-1.0d0,1.0d0,u,wu,n)
            do i = 1, n
                self%freqs(i) = self%freqmax*(u(i)+1.0d0)/2.0d0
                self%freqs(2*n-i+1) = 2.0d0*self%freqmax/(u(i)+1.0d0)
                self%womeg(i) = self%freqmax*wu(i)/2.0d0
                self%womeg(2*n-i+1) = 2.0d0*self%freqmax*wu(i)/((u(i)+1.0d0)*(u(i)+1.0d0))
            enddo ! i
            deallocate(u, wu)

        case('exp')
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            self%freqmin = 1.d-4
            t1 = ( log(self%freqmax) - log(self%freqmin) ) / dble(n-1) ! step in log scale
            self%freqs(1) = self%freqmin
            do i = 2, n
                t2 = log(self%freqmin) + dble(i-1)*t1
                self%freqs(i) = exp(t2)
            end do
            self%womeg(:) = 1.d0 / dble(n)

        case('cubic')
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            do i = 1, n
                self%freqs(i) = self%freqmin + &
                    (dble(i-1)/dble(n-1))**3 * (self%freqmax-self%freqmin)
                self%womeg(i) = 1.d0 / dble(n)
            end do

        case('mix')
            n = self%nomeg     ! number of points in the first sub-interval
            m = 32             ! number of points in the second sub-interval
            self%nomeg = n + m ! total number of points
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            ! uniform grid from 0 to w_max
            step = self%freqmax / dble(n)
            do i = 1, n
                self%freqs(i) = dble(i-1)*step
                self%womeg(i) = 1.d0/dble(n)
            enddo
            ! Gauss-Legendre grid from w_max to infinity
            allocate(u(m), wu(m))
            call gauleg(-1.0d0, 1.0d0, u, wu, m)
            do i = 1, m
                self%freqs(self%nomeg-i+1) = 2.0d0*self%freqmax/(u(i)+1.0d0)
                self%womeg(self%nomeg-i+1) = 2.0d0*self%freqmax*wu(i)/((u(i)+1.0d0)*(u(i)+1.0d0))
            end do ! i
            deallocate(u, wu)

        case('mix2')
            n = self%nomeg     ! number of points in the first sub-interval
            m = 32             ! number of points in the second sub-interval
            self%nomeg = n + m ! total number of points
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            ! uniform grid from 0 to w_max
            step = self%freqmax / dble(n)
            do i = 1, n
                self%freqs(i) = dble(i-1)*step
                self%womeg(i) = 1.d0/dble(n)
            enddo
            ! Double Gauss-Legendre grid from w_max to infinity
            fmax = self%freqmax + 1.d0
            allocate(u(m/2), wu(m/2))
            call gauleg(-1.d0, 1.d0, u, wu, m/2)
            do i = 1, m/2
                self%freqs(n+i) = self%freqmax+fmax*(u(i)+1.d0)/2.d0
                self%freqs(self%nomeg-i+1) = self%freqmax+2.0d0*fmax/(u(i)+1.d0)
                self%womeg(n+i) = fmax*wu(i)/2.d0
                self%womeg(self%nomeg-i+1) = 2.d0*fmax*wu(i)/((u(i)+1.d0)*(u(i)+1.d0))
            enddo ! i
            deallocate(u, wu)


        case('tanh')
            ! Hyperbolic tangent stretching
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            delta = self%freqmin
            do i = 1, n
                sigma = dble(i)/dble(n) - 1.d0
                self%freqs(i) = self%freqmax * (1.d0 + tanh(delta*sigma) / tanh(delta))
                self%womeg(i) = 1.d0/dble(n)
            end do

        case('clencurt1')
            ! Clenshaw-Curtis grid [-1,1]
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            ell = self%freqmax
            self%womeg(:) = 0.d0
            do i = 1, n
                ii = n-i+1
                t = pi * dble(i) / dble(n+1)
                self%freqs(ii) = ell * cos(t)
                do j = 1, n
                    self%womeg(ii) = self%womeg(ii) + sin(dble(j)*t) * (1.d0 - cos(dble(j)*pi) ) / dble(j)
                end do
                self%womeg(ii) = self%womeg(ii) * ell * sin(t) * 2.0/dble(n + 1)
            end do

        case('clencurt2')
            ! Clenshaw-Curtis grid [0,infinity]
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            ell = self%freqmax
            self%womeg(:) = 0.d0
            do i = 1, n
                ii = n-i+1
                t = pi * dble(i) / dble(n+1)
                self%freqs(ii) = ell * ( 1.d0/tan(0.5d0*t) )**2
                do j = 1, n
                    self%womeg(ii) = self%womeg(ii) + sin(dble(j)*t) * (1.d0 - cos(dble(j)*pi) ) / dble(j)
                end do
                self%womeg(ii) = self%womeg(ii) * 2.d0 * ell * sin(t) / (1.d0 - cos(t))**2 * 2.0/dble(n + 1)
            end do

        case('clencurt3')
            ! Clenshaw-Curtis grid [-infinity,infinity]
            allocate(self%freqs(self%nomeg))
            allocate(self%womeg(self%nomeg))
            n = self%nomeg
            ell = self%freqmax
            self%womeg(:) = 0.d0
            do i = 1, n
                ii = n-i+1
                t = pi * dble(i) / dble(n+1)
                self%freqs(ii) = ell * 1.d0/tan(t)
                self%womeg(ii) = ell * pi / (sin(t)**2 * dble(n+1))
            end do

        case default
            stop 'Error(mod_frequency::generate_freqgrid) Unknown grid type!'

        end select

        return
    end subroutine

    !-------------------------------------------------------------------------------
    subroutine print_freqgrid(self,funit)
        implicit none
        type(frequency), intent(IN) :: self
        integer,         intent(IN) :: funit
        ! local variables
        integer :: i
        call boxmsg(funit,'-','frequency grid')
        write(funit,*) 'Type: < fgrid > ', self%fgrid
        write(funit,*) 'Frequency axis: < fconv > ', self%fconv
        write(funit,*) 'Number of frequencies: < nomeg > ', self%nomeg
        write(funit,*) 'Cutoff frequency: < freqmax > ', self%freqmax
        write(funit,*) 'frequency list: < #    freqs    weight > '
        do i = 1, self%nomeg
            write(funit,'(i4,1p,2g18.10)') i, self%freqs(i), self%womeg(i)
        enddo
    end subroutine

END MODULE
