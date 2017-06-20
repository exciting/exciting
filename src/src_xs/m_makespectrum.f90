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
       & Making spectrum using formula for coupling with time inverted ar basis.")')&
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

      ! enw_{w, \lambda} = 1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta)
      allocate(enw(nfreq, nexc))

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(i,j)
      do j = 1, nexc
        do i = 1, nfreq
          enw(i,j) = zone/(bevalre(j)-freq(i)-zbrd)&
                  &+ zone/(bevalre(j)+freq(i)+zbrd)
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
      ! Make non-lattice-symmetrized spectrum      !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)

      allocate(ns_spectr(nfreq,nopt))

      ! qmt=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! qmt/=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda}(G,qmt) t_{\lambda}(G,qmt) 
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
      use modinput
      use mod_constants, only: zone, zi, pi
      use modxs, only: unitout
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
      logical :: useoff
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
        & Making spectrum using formula for coupling with time inverted ar basis.")')&
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
      ! (resonant & anti-resonant)                 !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Making energy denominators ENW.")')
      call timesec(t0)

      ! enw_{w, \lambda} = 1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta)
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
          denw%za(i,j) = zone/(bevalre(jg)-freq(ig)-zbrd)&
                       &+ zone/(bevalre(jg)+freq(ig)+zbrd)
#else
          denw%za(i,j) = zone/(bevalre(j)-freq(i)-zbrd)&
                       &+ zone/(bevalre(j)+freq(i)+zbrd)
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
      ! Make non-lattice-symmetrized spectrum      !
      !++++++++++++++++++++++++++++++++++++++++++++!
      write(unitout, '("  Calculating spectrum.")')
      call timesec(t0)

      call new_dzmat(dns_spectr, nfreq, nopt, binfo,&
        & rblck=binfo%mblck, cblck=1)

      ! qmt=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
      !       t^*_{\lambda, o1_j} t_{\lambda, o2_j} 
      ! qmt/=0 case:
      ! nsspectr_{w,j} = \Sum_{\lambda} enw_{w,\lambda} tmat_{\lambda, j}
      !   i.e. nsspectr_{w,j} = 
      !     \Sum_{\lambda} (1/(E_\lambda - w - i\delta) + 1/(E_\lambda + w + i\delta))
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

      integer(4), intent(in) :: iqmt, nk
      complex(8), intent(inout) :: nsp(:,:)
      complex(8), intent(out) :: sp(:,:,:)

      complex(8), allocatable :: buf(:,:,:)
      integer(4) :: nfreq, nopt, i, j, o1, o2, igqmt
      real(8) :: pref, t0, t1
      logical :: foff
      character(*), parameter :: thisname = "finalizespectrum"

      ! Check input
      nfreq = size(nsp,1) 
      nopt = size(nsp,2)

      ! Zero momentum transfer
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
      else

        if(nopt /= 1) then 
          write(*,'("Error(",a,"): nopt= ", i8, "while iqmt= ", i8)')&
            & trim(thisname), nopt, iqmt
          call terminate
        end if
        foff=.false.

      end if

      write(unitout, '("Info(",a,"): Finalizing spectrum.")') trim(thisname)
      call timesec(t0)

      ! Adjusting prfactor 
      if(iqmt==1) then 
        pref = 2.d0*4.d0*pi/omega/nk
      else
        ! 2 * ( (4pi/|Gmt+qmt|^2)^1/2 )^2 * 1/V * 1/nk
        igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
        pref = 2.0d0*sptclg(igqmt,iqmt)**2/omega/nk
      end if

      ! Add 1 to diagonal elements
      do j = 1, nopt
        if(foff) then
          !$OMP PARALLEL DO &
          !$OMP& DEFAULT(SHARED), PRIVATE(i)
          do i = 1, nfreq
            if(j == 1 .or. j == 5 .or. j == 9) then 
              nsp(i,j) = nsp(i,j)*pref+zone
            else
              nsp(i,j) = nsp(i,j)*pref
            end if
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO &
          !$OMP& DEFAULT(SHARED), PRIVATE(i)
          do i = 1, nfreq
            nsp(i,j) = nsp(i,j)*pref+zone
          end do
          !$OMP END PARALLEL DO
        end if
      end do

      ! Zero momentum transfer --> tensor structure
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

      ! Finite momentum transfer --> function
      else

        sp = zzero
        sp(1,1,:) = nsp(:,1)

      end if

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end subroutine finalizespectrum

end module m_makespectrum
