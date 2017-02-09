MODULE mod_frequency
    
    implicit none

!-------------------------------------------------------------------------------    
    type frequency
        character(8) :: fgrid           ! grid type
        character(8) :: fconv           ! type of the frequency dependence
        integer :: nomeg                ! grid size
        real(8) :: freqmax              ! cutoff frequency
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
    subroutine generate_freqgrid(self,fgrid,fconv,nomeg,freqmax)
        implicit none
        type(frequency), intent(OUT) :: self
        character*(*),   intent(IN)  :: fgrid
        character*(*),   intent(IN)  :: fconv
        integer,         intent(IN)  :: nomeg
        real(8),         intent(IN)  :: freqmax
! local variables
        integer :: i, n
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
        
! number of frequencies
        if ((self%fgrid=='gaule2').or.(self%fgrid=='GAULE2')) then
            ! always even for 'gaule2' grid type
            n = nomeg-mod(nomeg,2)
            self%nomeg = n
        else
            self%nomeg = nomeg
        end if

! cutoff frequency
        self%freqmax = freqmax

! generate frequency integration grid
        if (allocated(self%freqs)) deallocate(self%freqs)
        allocate(self%freqs(self%nomeg))
        if (allocated(self%womeg)) deallocate(self%womeg)
        allocate(self%womeg(self%nomeg))

        select case (self%fgrid) 

        case('eqdist','EQDIST') ! Equaly spaced mesh (for tests purposes only)
            do i = 1, self%nomeg
                self%freqs(i) = dble(i)*self%freqmax/dble(self%nomeg)
                self%womeg(i) = 1.d0/dble(self%nomeg)
            enddo  
      
        case('gaulag','GAULAG') ! Grid for Gauss-Laguerre quadrature
            allocate(wu(self%nomeg))
            call gaulag(self%freqs,wu,self%nomeg,-1.0d0)
            do i = 1, self%nomeg
                self%womeg(i) = wu(i)*dexp(self%freqs(i))*self%freqs(i)
            enddo  
            deallocate(wu)

        case('gaule2','GAULE2') ! Grid for Gauss-Legendre quadrature from 0 to freqmax 
                                ! +
                                ! Gauss-Legendre quadrature from freqmax to infinity
            n = self%nomeg/2
            allocate(u(n))
            allocate(wu(n))
            call gauleg(-1.0d0,1.0d0,u,wu,n)
            do i = 1, n
                self%freqs(i) = self%freqmax*(u(i)+1.0d0)/2.0d0
                self%freqs(2*n-i+1) = 2.0d0*self%freqmax/(u(i)+1.0d0)
                self%womeg(i) = self%freqmax*wu(i)/2.0d0
                self%womeg(2*n-i+1) = 2.0d0*self%freqmax*wu(i)/((u(i)+1.0d0)*(u(i)+1.0d0))
            enddo ! i  
            deallocate(u)
            deallocate(wu)
          
        case('gauleg','GAULEG') ! Grid for Gauss-Legendre quadrature from 0 to omegamax
            call gauleg(0.0d0,self%freqmax,self%freqs,self%womeg,self%nomeg)

        case default            
            write(*,*) 'ERROR(mod_frequency::generate_freqgrid) Unknown frequency grid type!'
            write(*,*) '  Available options:'
            write(*,*) '  - eqdist : equidistant frequencies from 0 to freqmax'
            write(*,*) '  - gaulag : grid for Gauss-Laguerre quadrature from 0 to infinity'
            write(*,*) '  - gauleg : grid for Gauss-Legendre quadrature from 0 to freqmax'
            write(*,*) '  - gaule2 : grid for Gauss-Legendre quadrature from 0 to freqmax and from freqmax to infinity'
            stop
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
            case default
                write(funit,*) 'ERROR(mod_frequency::generate_freqgrid) Unknown frequency grid type!'
                stop
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
