module m_makespectrum

  implicit none

  contains

    ! Non-distributed versions

    subroutine makespectrum(iqmt, nexc, nk, bevalre, oscsr, spectrum)
      use modinput
      use mod_constants
      use modxs, only: unitout
      use modmpi
      use m_genwgrid

      implicit none

      ! I/O
      integer(4), intent(in) :: iqmt, nexc, nk
      real(8), intent(in) :: bevalre(:)
      complex(8), intent(in) :: oscsr(:,:)
      complex(8), allocatable, intent(inout) :: spectrum(:,:,:)

      ! Local
      integer(4) :: nfreq
      real(8), allocatable :: freq(:)
      logical :: useoff
      integer(4) :: i, j, o1, o2, nopt
      real(8) :: t1, t0, ts0, ts1
      complex(8) :: zbrd
      complex(8), allocatable :: ns_spectr(:,:)
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

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      !++++++++++++++++++++++++++++++++++++++++++++!
      
      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Helper matrix build form oscillator        !
      ! strengths                                  !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Making helper matrix tmat.")')
      call timesec(t0)

      if(iqmt == 1) then 

        if(useoff) then
          nopt = 9
        else
          nopt = 3
        end if

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
            tmat(i,j) = conjg(oscsr(i,o1))*oscsr(i,o2)
          end do
          !$OMP END PARALLEL DO
        end do

      else

        nopt = 1

        ! tmat_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
        allocate(tmat(nexc,nopt))
        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i)
        do i = 1, nexc
          tmat(i,1) = conjg(oscsr(i,1))*oscsr(i,1)
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

      ! Post process non-lattice-symmetrized spectrum
      call finalizespectrum(iqmt, nk, ns_spectr, spectrum)
      deallocate(ns_spectr)

      ! Total time for spectrum construction
      call timesec(ts1)
      write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0
    end subroutine makespectrum

    ! Distributed versions

    subroutine makespectrum_dist(iqmt, nexc, nk, bevalre, oscsr, symsp, binfo)
      use mod_constants, only: zone, zi, pi
      use modxs, only: unitout
      use modinput, only: input
      use modmpi
      use modscl
      use m_genwgrid
      use m_dzmatmult

      implicit none

      ! I/O
      integer(4), intent(in) :: iqmt, nexc, nk
      real(8), intent(in) :: bevalre(:)
      complex(8), intent(in) :: oscsr(:,:)
      complex(8), intent(inout) :: symsp(:,:,:)
      type(blacsinfo) :: binfo

      ! Local
      integer(4) :: nfreq
      real(8), allocatable :: freq(:)
      integer(4) :: i, j, o1, o2, ig, jg, nopt
      real(8) :: t1, t0, ts0, ts1
      complex(8) :: zbrd
      complex(8), allocatable :: ns_spectr(:,:)
      logical :: useoff, usechibar
      type(dzmat) :: denw, dtmat, dns_spectr

      character(*), parameter :: thisname = "makespectrum_dist"

      ! Allocate frequency array used in spectrum construction
      nfreq = input%xs%energywindow%points
      allocate(freq(nfreq))

      ! Generate an evenly spaced frequency grid 
      call genwgrid(nfreq, input%xs%energywindow%intv,&
        & input%xs%tddft%acont, 0.d0, w_real=freq)

      ! Offdiagonals of eps_m
      useoff = input%xs%dfoffdiag
      ! Broadending factor
      zbrd = zi*input%xs%broad

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
      call timesec(ts0)

      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Make energy denominator for each frequency !
      ! The signs of the complex shifts are such   !
      ! that the response is retarded.             !
      ! (resonant & anti-resonant)                 !
      ! resonant = poles for positive omega        !
      ! anti-resonant = poles for negative omega   !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)

      ! enw_{w, \lambda} = 1/(w - E_\lambda + i\delta) + 1/(-w - E_\lambda - i\delta)
      call new_dzmat(denw, nfreq, nexc, binfo)

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
#ifdef SCAL
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j,ig,jg)
#else
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
#endif
      do j = 1, denw%ncols_loc
        do i = 1, denw%nrows_loc
#ifdef SCAL
          ! Get corresponding global indices
          ig = denw%r2g(i)
          jg = denw%c2g(j)
          denw%za(i,j) = zone/(freq(ig)-bevalre(jg)+zbrd)&
                       &+ zone/(-freq(ig)-bevalre(jg)-zbrd)
#else
          denw%za(i,j) = zone/(freq(i)-bevalre(j)+zbrd)&
                       &+ zone/(-freq(i)-bevalre(j)-zbrd)
#endif
        end do
      end do
      !$OMP END PARALLEL DO

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      !++++++++++++++++++++++++++++++++++++++++++++!
      
      !++++++++++++++++++++++++++++++++++++++++++++!
      ! Helper matrix build form resonant          !
      ! oscillator strenghs                        !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Making helper matrix tmat.")')
      call timesec(t0)

      if(iqmt == 1) then 
        ! tmat_{\lambda, j} = t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
        ! where j combines the Cartesian directions
        if(useoff) then
          nopt = 9
        else
          nopt = 3
        end if
      else
        ! tmat_{\lambda}(G,qmt) = t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
        nopt = 1 
      end if

      call new_dzmat(dtmat, nexc, nopt, binfo,&
        & rblck=binfo%mblck, cblck=1)

      do j = 1, dtmat%ncols_loc 
#ifdef SCAL
        ! Get corresponding global indices
        jg = dtmat%c2g(j)
#else
        jg = j
#endif
        ! Get individual opical indices
        if(iqmt == 1) then 
          if(useoff) then
            o2 = (jg-1)/3 + 1
            o1 = jg-(o2-1)*3
          else
            o2 = jg
            o1 = jg
          end if
        ! Finite qmt is longitudinal only
        else
          o2=1
          o1=1
        end if

        !$OMP PARALLEL DO &
        !$OMP& DEFAULT(SHARED), PRIVATE(i,ig)
        do i = 1, dtmat%nrows_loc
#ifdef SCAL
          ! Get corresponding global indices
          ig = dtmat%r2g(i)
#else
          ig = i
#endif
          dtmat%za(i,j) = conjg(oscsr(ig,o1))*oscsr(ig,o2)
        end do
        !$OMP END PARALLEL DO
      end do

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

      call new_dzmat(dns_spectr, nfreq, nopt, binfo,&
        & rblck=binfo%mblck, cblck=1)

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
      call dzmatmult(denw, dtmat, dns_spectr)
      !++++++++++++++++++++++++++++++++++++++++++++!

      ! Helper no longer needed
      call del_dzmat(dtmat)
      call del_dzmat(denw)
      deallocate(freq)

      ! Send spectrum to a global matrix at rank 0 to 
      ! interface with non parallel post-processing routines.
      call dzmat_send2global_root(ns_spectr, dns_spectr, binfo)
      call del_dzmat(dns_spectr)

      ! Postprocess non-lattice-symmetrized spectrum
      if(binfo%isroot) then
        call finalizespectrum(iqmt, nk, ns_spectr, symsp)
        deallocate(ns_spectr)
      end if

      ! Total time for spectrum construction
      call timesec(ts1)
      write(unitout, '("  Spectrum made in", f12.3, "s")') ts1-ts0

    end subroutine makespectrum_dist

    subroutine finalizespectrum(iqmt, nk, nsp, sp)
      use modmpi
      use modxs, only: sptclg, ivgigq, ivgmt, unitout
      use mod_lattice, only: omega
      use modxs, only: symt2
      use mod_constants
      use modinput, only: input

      integer(4), intent(in) :: iqmt, nk
      complex(8), intent(inout) :: nsp(:,:)
      complex(8), intent(out) :: sp(:,:,:)

      complex(8), allocatable :: buf(:,:,:)
      integer(4) :: nfreq, nopt, i, j, o1, o2, igqmt
      real(8) :: pref, t0, t1
      logical :: foff, usechibar, fcoup
      character(*), parameter :: thisname = "finalizespectrum"

      ! For BSE with coupling at finite Q, always use the full Chi
      ! and not \bar{Chi}. For TDA warn if \bar{Chi} is not used.
      fcoup = input%xs%bse%coupling
      usechibar = input%xs%bse%chibarq
      ! If it is not the default of true
      if(.not. usechibar) then 
        if(.not. fcoup .and. .not. usechibar) then
          write(unitout, '("Waring(",a,"):", a)') trim(thisname),&
            & " TDA and Chi not compatible!"
        end if
      ! Otherwise determin value with fcoup
      else
        if(fcoup) usechibar = .false.
        if(.not. fcoup) usechibar = .true.
      end if

      ! Check input
      nfreq = size(nsp,1) 
      nopt = size(nsp,2)

      ! Zero momentum transfer:
      !
      !   Using \bar{P} formalism:
      !     epsilon^{ij}_M(w) = 
      !       \delta_{ij} - 4*pi 1/(V*nk) 
      !       * \sum_\lambda [1/(w-E_\lambda+i\delta) + 1/(-w-E_\lambda-i\delta)]
      !       * t^*_{\lambda,i} t_{\lambda,j}
      !  
      !   nsp(i,j) = \delta_{ij} + pref * nsp(i,j)
      if(iqmt==1) then 

        if(nopt == 9) then 
          foff = .true.
        else if(nopt == 3) then
          foff = .false.
        else
          write(*,'("Error(",a,"): nopt invalid", i8)') trim(thisname), nopt
          call terminate
        end if

      ! Finite momentum transfer
      !
      !   Using \chi formalism:
      !     epsilon_M(\vec{Q},w) = { 
      !       1 + 4*pi/|Q|^2 * 1/(V*nk) 
      !       * \sum_\lambda [1/(w - E_\lambda + i\delta) + 1/(-w - E_\lambda - i\delta)]
      !       * t^*_{\lambda}(Q) t_{\lambda}(Q)
      !       }^{-1}
      !
      !   OR use also \bar{P} formalism for finite q, needs zeroing of G+q
      !   component of coulomb interaction
      !  
      !   nsp(i,j) = 1 + pref' * nsp(i,j)
      else

        ! For finite qmt, only the component along q is calculated 
        if(nopt /= 1) then 
          write(*,'("Error(",a,"): nopt= ", i8, "while iqmt= ", i8)')&
            & trim(thisname), nopt, iqmt
          call terminate
        end if
        foff=.false.

      end if

      write(unitout, '("Info(",a,"): Finalizing spectrum.")') trim(thisname)
      if(usechibar) then 
        write(unitout, '("Info(",a,"): Using modified Chi to build spectrum.")') trim(thisname)
      end if
      call timesec(t0)

      ! Adjusting prfactor, the factor 2 accounts for spin degeneracy
      if(iqmt==1) then 
        ! -1 * 2 * 4 pi * 1/V * 1/nk
        pref = -2.d0*4.d0*pi/omega/nk
      else
        ! 2 * ( (4 pi/|Gmt+qmt|^2)^1/2 )^2 * 1/V * 1/nk
        igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
        pref = 2.0d0*sptclg(igqmt,iqmt)**2/omega/nk
        if(usechibar) then 
          pref = -pref
        end if
      end if

      ! Add 1 to diagonal elements
      do j = 1, nopt
        if(foff) then
          !$OMP PARALLEL DO &
          !$OMP& DEFAULT(SHARED), PRIVATE(i)
          do i = 1, nfreq
            if(j == 1 .or. j == 5 .or. j == 9) then 
              nsp(i,j) = zone + nsp(i,j)*pref
            else
              nsp(i,j) = nsp(i,j)*pref
            end if
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO &
          !$OMP& DEFAULT(SHARED), PRIVATE(i)
          do i = 1, nfreq
            nsp(i,j) = zone + nsp(i,j)*pref
          end do
          !$OMP END PARALLEL DO
        end if
      end do

      ! Zero momentum transfer --> tensor structure
      !   Symmetrize using the crystal symmetries
      if(iqmt == 1) then 

        ! Write to buffer for symmetry routine
        allocate(buf(3,3,nfreq))
        buf=zzero
        do j = 1, nopt
          if(foff) then
            o2 = (j-1)/3 + 1
            o1 = j-(o2-1)*3
          else
            o2 = j
            o1 = j
          end if
          !$OMP PARALLEL DO &
          !$OMP& DEFAULT(SHARED), PRIVATE(i)
          do i = 1, nfreq
            buf(o1,o2,i) = nsp(i,j)
          end do
          !$OMP END PARALLEL DO
        end do
      
        ! Symmetrize spectrum with respect to the crystal symmetry
        sp = zzero
        do o1=1,3
          do o2=1,3
            ! Symmetrize the macroscopic dielectric tensor
            call symt2app(o1, o2, nfreq, symt2, buf, sp(o1,o2,:))
          end do 
        end do

      ! Finite momentum transfer --> scalar function
      else

        sp = zzero

        if(usechibar) then 
          ! Use chibar also for finte q, needs zeroing of G+q in the coulomb
          ! interaction
          sp(1,1,:) = nsp(:,1)
        else
          ! As noted above nsp contains 1/epsm(Q,w), save epsm
          sp(1,1,:) = 1.0d0/nsp(:,1)
        end if

      end if

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end subroutine finalizespectrum

end module m_makespectrum
