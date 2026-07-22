module m_makespectrum

  implicit none

  contains

    !> Read the data of EXCITON*.OUT files and compute the spectrum
    subroutine makespectrum(iqmt, nexc, nk, bevalre, oscsr, spectrum)
      use modinput
      use constants
      use modxs, only: unitout
      use modmpi
      use m_genwgrid

      implicit none

      ! I/O
      integer(4), intent(in) :: iqmt, nexc, nk
      real(8), intent(in) :: bevalre(:)
      complex(8), intent(in) :: oscsr(:,:,:)
      complex(8), allocatable, intent(out) :: spectrum(:,:,:)

      ! Local
      integer(4) :: nfreq
      real(8), allocatable :: freq(:)
      logical :: useoff, usechibar, fcoup
      integer(4) :: i, j, o1, o2, nopt, comp
      real(8) :: t1, t0, ts0, ts1, sign
      complex(8) :: zbrd
      complex(8), allocatable :: ns_spectr(:,:), buf(:)
      complex(8), allocatable :: enw(:,:), tmat(:,:)

      character(*), parameter :: thisname = "makespectrum"

      ! Allocate frequency array used in spectrum construction
      nfreq = input%xs%energywindow%points
      allocate(freq(nfreq))
      if(allocated(spectrum)) deallocate(spectrum)
      allocate(spectrum(3,3,nfreq))
      spectrum = zzero

      ! Generate an evenly spaced frequency grid 
      call genwgrid(nfreq, input%xs%energywindow%intv,&
        & input%xs%tddft%acont, 0.d0, w_real=freq)

      ! Offdiagonals of eps_m
      useoff = input%xs%dfoffdiag
      ! Broadending factor
      zbrd = zi*input%xs%broad
      
      fcoup = input%xs%bse%coupling
      usechibar = input%xs%bse%chibarq
      if(.not. usechibar) then 
        if(.not. fcoup .and. .not. usechibar) then
          write(unitout, '("Warning(",a,"):", a)') trim(thisname),&
            & " TDA and Chi not compatible!"
        end if
      ! Otherwise determine value with fcoup
      else
        usechibar = .not. fcoup
      end if

      write(unitout, '("Info(",a,"):&
       & Making spectrum using formula for coupling with time reversed ar basis.")')&
       & trim(thisname)
      if(iqmt == 1) then 
        write(unitout, '("Info(",a,"):&
         & Using formalism for vaninshing momentum transfer.")') trim(thisname)
        if(useoff) then
          write(unitout, '("Info(",a,"): Including off diagonals.")') trim(thisname)
        end if
      else
        write(unitout, '("Info(",a,"):&
         & Using formalism for finite momentum transfer. iqmt=",i4)')&
         & trim(thisname), iqmt
      end if

      ! Total time for spectrum construction
      call timesec(ts0)

      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Make energy denominator for each frequency !
      ! (resonant & anti-resonant)                 !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)

      ! enw_{w, \lambda} = 1/(w - E_\lambda + i\delta) + 1/(-w - E_\lambda - i\delta)
      allocate(enw(nfreq, nexc))

      ! duplication of the do loop to help better vectorization
      if (input%xs%BSE%aresbse .and. (fcoup .or. input%xs%BSE%chibar0)) then
        !$OMP PARALLEL DO &
        !$OMP& COLLAPSE(2),&
        !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
        do j = 1, nexc
          do i = 1, nfreq
            enw(i,j) = zone/(freq(i) - bevalre(j) + zbrd)&
                    &+ zone/(-freq(i) -bevalre(j) - zbrd)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO &
        !$OMP& COLLAPSE(2),&
        !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
        do j = 1, nexc
          do i = 1, nfreq
            enw(i,j) = zone/(freq(i) - bevalre(j) + zbrd)
          end do
        end do
        !$OMP END PARALLEL DO
      end if

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      !++++++++++++++++++++++++++++++++++++++++++++!
      
      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Helper matrix build form oscillator        !
      ! strengths                                  !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Making helper matrix tmat.")')
      call timesec(t0)

      sign = 1.0
      if (input%xs%bse%chibarq) sign = -1.0  ! accounts for sign in `writeoscillator`

      if(iqmt == 1) then 
        if(useoff) then
          nopt = 9
        else
          nopt = 3
        end if

        if (.not. input%xs%BSE%chibar0) sign = -sign

        ! tmat_{\lambda, j} = t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
        ! where j combines the Cartesian directions
        allocate(tmat(nexc,nopt))
        do j = 1, nopt
          if(useoff) then
            o2 = (j-1)/3 + 1
            o1 = j-(o2-1)*3
          else
            o2 = j
            o1 = j
          end if
          !$OMP PARALLEL DO &
          !$OMP& DEFAULT(SHARED), PRIVATE(i)
          do i = 1, nexc
            tmat(i,j) = sign*oscsr(i,o1,o2)
          end do
          !$OMP END PARALLEL DO
        end do

      else

        nopt = 1
        if (usechibar) sign = -sign

        ! tmat_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
        allocate(tmat(nexc,nopt))
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nexc
          tmat(i,1)= sign*oscsr(i,1,1)
        end do
        !$OMP END PARALLEL DO

      end if

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      !++++++++++++++++++++++++++++++++++++++++++++!
           
      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Make non-lattice-symmetrized response functions
      ! Chi_{Gmt,Gmt}(qmt,omega) for qmt/=0,
      ! \bar{P}^{ij}_{00}(0,omega) for qmt=0
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)

      allocate(ns_spectr(nfreq,nopt))

      ! qmt=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(w - E_\lambda + i\delta) + 1/(-w - E_\lambda - i\delta))
      !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! qmt/=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(w - E_\lambda + i\delta) + 1/(-w - E_\lambda - i\delta))
      !       t^*_{\lambda}(Gmt,qmt) t_{\lambda}(Gmt,qmt) 
      call zgemm('N','N', nfreq, nopt, nexc, zone, enw, nfreq,&
        & tmat, nexc, zzero, ns_spectr, nfreq)
      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Helper no longer needed
      deallocate(tmat, enw)

      ! add one to diagonal terms
      do j=1, nopt
        if (useoff) then
          o2=(j-1)/3+1
          o1=j-(o2-1)*3
          do i=1, nfreq
            if (j==1 .or. j==5 .or. j==9) then
              spectrum(o1,o2,i)=zone+ns_spectr(i,j)
            else
              spectrum(o1,o2,i)=ns_spectr(i,j)
            end if
          end do
        else
          o1=j
          o2=j
          do i=1,nfreq
            spectrum(o1,o2,i)=zone+ns_spectr(i,j)
          end do
        end if
      end do

      deallocate(ns_spectr)
      if (allocated(buf)) deallocate(buf)
      allocate(buf(nfreq))
      buf(:) = spectrum(1,1,:)

      ! no momentum transfer
      if (iqmt == 1 .and. .not. input%xs%BSE%chibar0) then
        spectrum = zzero
        comp = input%xs%BSE%chibar0comp
        spectrum(comp,comp,:) = 1.0d0 / buf(:)
      end if

      ! for finite momentum transfer
      if (iqmt /= 1 .and. .not. usechibar) then
        spectrum(1,1,:) = 1.0d0 / buf(:)
      end if

      deallocate(buf)

      ! Total time for spectrum construction
      call timesec(ts1)
      write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0
    end subroutine makespectrum

end module m_makespectrum
