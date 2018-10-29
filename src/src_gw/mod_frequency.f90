MODULE mod_frequency
    
    implicit none

!-------------------------------------------------------------------------------    
    type frequency
        character(8) :: fgrid           ! grid type
        character(8) :: fconv           ! type of the frequency dependence
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
        integer :: i, n
        real(8) :: t1, t2, step
        real(8), allocatable :: u(:)
        real(8), allocatable :: wu(:)

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
            write(*,*) '  - refreq : weights calculated for real frequecies'
            write(*,*) '  - imfreq : weights calculated for imaginary frequecies'
            stop
        end select
        
! number of frequencies is always even
        n = nomeg - mod(nomeg,2)
        self%nomeg = n

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
            ! non-uniform grid from w_max to infty
            do i = 1, n
                self%freqs(n+i) = self%freqmax + dble(i)*step + &
                    (dble(i-1)/dble(n-1))**3 * (100.d0-self%freqmax-dble(n)*step)
                self%womeg(n+i) = 1.d0/dble(self%nomeg)
            end do
        
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
