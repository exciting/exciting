MODULE mod_frequency
    
    implicit none

!-------------------------------------------------------------------------------    
    type frequency
        character(40) :: fgrid          ! grid type
        character(40) :: fconv          ! type of the frequency dependence
        integer :: nomeg                ! grid size
        real(8) :: freqmin              ! lower cutoff frequency
        real(8) :: freqmax              ! upper cutoff frequency
        real(8), allocatable :: freqs(:)! frequency grid
        real(8), allocatable :: womeg(:)! integration weights
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
        real(8) :: t, ell
        real(8), allocatable :: u(:), x(:)
        real(8), allocatable :: wu(:), w1(:), w2(:)
        real(8), parameter   :: pi = 3.1415926535897932385d0

        ! grid type        
        self%fgrid = trim(fgrid)

        ! frequency dependence type
        self%fconv = trim(fconv)
        select case (self%fconv)
        case('nofreq')
            self%nomeg = 1
        case('refreq')
        case('imfreq')
        case default
            write(*,*) 'ERROR(mod_frequency::generate_freqgrid):&
            &  Unknown frequency convolution method!'
            write(*,*) '  Currently supported options are:'
            write(*,*) '  - nofreq : no frequecy dependence of the weights'
            write(*,*) '  - refreq : weights calculated for real frequencies'
            write(*,*) '  - imfreq : weights calculated for imaginary frequencies'
            stop
        end select
        
        ! number of frequencies (always even)
        self%nomeg = nomeg - mod(nomeg,2)

        ! cutoff frequencies
        self%freqmin = freqmin
        self%freqmax = freqmax

        ! generate frequency integration grid
        if (allocated(self%freqs)) deallocate(self%freqs)
        allocate(self%freqs(self%nomeg))
        if (allocated(self%womeg)) deallocate(self%womeg)
        allocate(self%womeg(self%nomeg))

        select case (self%fgrid) 

        case('eqdist','EQDIST') ! Equaly spaced mesh (for tests purposes only)
            do i = 1, self%nomeg
                self%freqs(i) = self%freqmin + dble(i-1)*(self%freqmax-self%freqmin)/dble(self%nomeg-1)
                self%womeg(i) = 1.d0/dble(self%nomeg)
            enddo  
      
        case('gaulag','GAULAG') ! Grid for Gauss-Laguerre quadrature
            allocate(wu(self%nomeg))
            call gaulag(self%freqs, wu, self%nomeg, -1.0d0)
            do i = 1, self%nomeg
                self%womeg(i) = wu(i)*dexp(self%freqs(i))*self%freqs(i)
            enddo  
            deallocate(wu)

        case('gaule2','GAULE2') ! Grid for Gauss-Legendre quadrature from 0 to freqmax 
                                ! +
                                ! Gauss-Legendre quadrature from freqmax to infinity
            self%nomeg = nomeg - mod(nomeg,2)
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
          
        case('gauleg','GAULEG') ! Grid for Gauss-Legendre quadrature from 0 to omegamax
            call gauleg(0.0d0,self%freqmax,self%freqs,self%womeg,self%nomeg)

        case('GL2') ! Grid for Gauss-Legendre quadrature from 0 to infinity
            n = self%nomeg
            allocate( u(n), wu(n) )
            call gauleg(0.d0, 1.d0, u, wu, n)
            do i = 1, n
                self%freqs(i) = u(i) / ( 1.d0 - u(i) )
                self%womeg(i) = wu(i) * ( 1.d0 - u(i) )**(-2)
            end do
            deallocate(u, wu)

        case('exp')
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
            n = self%nomeg
            do i = 1, n
                self%freqs(i) = self%freqmin + &
                    (dble(i-1)/dble(n-1))**3 * (self%freqmax-self%freqmin)
                self%womeg(i) = 1.d0 / dble(n)
            end do

        case('mix')
            n = self%nomeg/2
            step = (self%freqmax-self%freqmin) / dble(n-1)
            ! uniform grid from 0 to w_max
            do i = 1, n
                self%freqs(i) = self%freqmin + dble(i-1)*step
                self%womeg(i) = 1.d0/dble(self%nomeg)
            enddo
            ! non-uniform grid from w_max to infinity
            do i = 1, n
                self%freqs(n+i) = self%freqmax + dble(i)*step + &
                    (dble(i-1)/dble(n-1))**3 * (100.d0-self%freqmax-dble(n)*step)
                self%womeg(n+i) = 1.d0/dble(self%nomeg)
            end do

        case('tanh')
            ! Hyperbolic tangent stretching
            n = self%nomeg
            delta = 6.0
            do i = 1, n
                sigma = dble(i)/dble(n) - 1.d0
                self%freqs(i) = self%freqmax * (1.d0 + tanh(delta*sigma) / tanh(delta))
                self%womeg(i) = 1.d0/dble(n)
            end do

        case('kron2')
            ! Order of the quadrature (only odd)
            n = nomeg - mod(nomeg+1,2)
            m = (n-1)/2
            allocate(x(m+1), w1(m+1), w2(m+1))
            call kronrod(m, 1.d-8, x, w1, w2)
            do i = 1, m+1
                write(*,'(i4,2f12.4)') i, x(i), w1(i)
            end do
            write(*,*)
            
            ! Complete set in [-1,1] that contains 2M+1 points
            allocate(u(n), wu(n))
            do i = 1, m+1
                u(i) = -x(i);  u(i+m+1) = x(m-i+1)
                wu(i) = w1(i); wu(i+m+1) = w1(m-i+1)
            end do
            deallocate(x, w1, w2)

            self%nomeg = 2*n
            if (allocated(self%freqs)) deallocate(self%freqs)
            allocate(self%freqs(self%nomeg))
            if (allocated(self%womeg)) deallocate(self%womeg)
            allocate(self%womeg(self%nomeg))

            do i = 1, n
                self%freqs(i) = self%freqmax*(u(i)+1.0d0)/2.0d0
                self%freqs(2*n-i+1) = 2.0d0*self%freqmax/(u(i)+1.0d0)
                self%womeg(i) = self%freqmax*wu(i)/2.0d0
                self%womeg(2*n-i+1) = 2.0d0*self%freqmax*wu(i)/((u(i)+1.0d0)*(u(i)+1.0d0))
            enddo ! i
            deallocate(u, wu)

        case('clencurt2')
            ! order of the Clenshaw-Curtis rule + 1
            n = nomeg+1
            allocate(u(n), wu(n))
            call clenshaw_curtis_compute(n, u, wu)
            write(*,*)
            do i = 1, n
                write(*,*) i, u(i), wu(i)
            end do
            write(*,*)
            self%nomeg = 2*(n-2)
            if (allocated(self%freqs)) deallocate(self%freqs)
            allocate(self%freqs(self%nomeg))
            if (allocated(self%womeg)) deallocate(self%womeg)
            allocate(self%womeg(self%nomeg))
            do i = 2, n-1
                self%freqs(i-1) = self%freqmax*(u(i)+1.d0)/2.0d0
                self%freqs(2*n-i-2) = 2.0d0*self%freqmax/(u(i)+1.d0)
                self%womeg(i-1) = self%freqmax*wu(i)/2.0d0
                self%womeg(2*n-i-2) = 2.0d0*self%freqmax*wu(i)/((u(i)+1.d0)*(u(i)+1.d0))
            enddo ! i
            deallocate(u, wu)

        case('clencurt')
            ell = self%freqmax
            n = nomeg
            self%nomeg = n
            if (allocated(self%freqs)) deallocate(self%freqs)
            allocate(self%freqs(self%nomeg))
            if (allocated(self%womeg)) deallocate(self%womeg)
            allocate(self%womeg(self%nomeg))
            self%womeg(:) = 0.d0
            do i = 1, n
                ii = n-i+1
                t = pi * dble(i) / dble(n+1)
                self%freqs(ii) = ell * ( cotan(0.5d0*t) )**2
                do j = 1, n
                    self%womeg(ii) = self%womeg(ii) + sin(dble(j)*t) * (1.d0 - cos(dble(j)*pi) ) / dble(j)
                end do
                self%womeg(ii) = self%womeg(ii) * 2.d0 * ell * sin(t) / ( 1.d0 - cos(t) )**2 * 2.0 / dble(n + 1)
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
        write(funit,*) 'Type: < fgrid >'
        select case (self%fgrid)
            case('eqdist','EQDIST')
                write(funit,*) '  - Equaly spaced mesh (for tests purposes only)'
            case('gaulag','GAULAG')
                write(funit,*) '  - Grid for Gauss-Laguerre quadrature'
            case('gaule2','GAULE2')
                write(funit,*) '  - Grid for double Gauss-Legendre quadrature,'
                write(funit,*) '    from 0 to freqmax and from freqmax to infinity'
            case('gauleg','GAULEG')
                write(funit,*) '  - Grid for Gauss-Legendre quadrature, from 0 to freqmax'
            case('GL2')
                write(funit,*) '  - Gauss-Legendre from 0 to infinity'

        end select
        write(funit,*) 'Convolution method: < fconv >'
        select case (self%fconv)
            case('nofreq')
                write(funit,*) '  - nofreq : no frequecy dependence of the weights'
            case('refreq')
                write(funit,*) '  - refreq : weights calculated for real frequecies'
            case('imfreq')
                write(funit,*) '  - imfreq : weights calculated for imaginary frequecies'
            case default
                write(*,*) 'ERROR(mod_frequency::generate_freqgrid) &
               & Unknown frequency convolution method!'
                stop
        end select
        write(funit,*) 'Number of frequencies: < nomeg >', self%nomeg
        write(funit,*) 'Cutoff frequency: < freqmax >', self%freqmax
        write(funit,*) 'frequency list: < #    freqs    weight >'
        do i = 1, self%nomeg
            write(funit,'(i4,1p,2g18.10)') i, self%freqs(i), self%womeg(i)
        enddo    
        call linmsg(funit,'-','')   

    end subroutine

END MODULE
