
module mod_selfenergy

    use mod_frequency

    type(frequency) :: freq_selfc

    !--------------!
    ! self-energy  !
    !--------------!

    ! The exchange self-energy
    complex(8), allocatable :: selfex(:,:)

    ! Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
    complex(8), allocatable :: mwm(:,:,:)
    target mwm

    ! The correlation self-energy
    complex(8), allocatable :: selfeph(:,:,:)
    complex(8), allocatable :: selfeph0(:,:)
    real(8),    allocatable :: speceph(:,:,:)
    complex(8), allocatable :: selfec(:,:,:)

    ! Correction factors for (q^-1) and (q^-2) singularities
    real(8) :: singc1
    real(8) :: singc2

    !-------------!
    ! QP Energy   !
    !-------------!

    ! Original KS energies (evalfv will updated via self-consistent cycle)
    real(8) :: eferks
    real(8), allocatable :: evalks(:,:)

    ! QP energies
    real(8) :: eferqp
    real(8), allocatable :: evalqp(:,:)

    ! Chemical potential alignment
    real(8) :: deltaE

    ! Linearization (renormalization) factor
    real(8),    allocatable :: znorm(:,:)

    ! AC to the real axis of the correlation self-energy (selfec)
    complex(8), allocatable :: sigc(:,:)

    ! COHSEX approximation
    complex(8), allocatable :: sigsx(:,:) ! Screened exchange
    complex(8), allocatable :: sigch(:,:) ! Coulomb hole

contains

    !---------------------------------------------------------------------------
    subroutine init_selfenergy(ibgw,nbgw,nkpt)
        use modinput
        implicit none
        integer, intent(in) :: ibgw, nbgw
        integer, intent(in) :: nkpt
        ! local
        integer(4) :: nw

        ! KS eigenvalues
        if (allocated(evalks)) deallocate(evalks)
        allocate(evalks(ibgw:nbgw,nkpt))
        evalks(:,:) = 0.d0

        ! Quasi-Particle energy
        if (allocated(evalqp)) deallocate(evalqp)
        allocate(evalqp(ibgw:nbgw,nkpt))
        evalqp(:,:) = 0.d0

        ! Exchange self-energy
        if (allocated(selfex)) deallocate(selfex)
        allocate(selfex(ibgw:nbgw,nkpt))
        selfex(:,:) = 0.d0

        ! Correlation self-energy
        if (input%gw%selfenergy%method == 'cd') then
          if ( .not.associated(input%gw%selfenergy%wgrid) ) &
              input%gw%selfenergy%wgrid => getstructwgrid(emptynode)
          call generate_freqgrid(freq_selfc, &
                                 input%gw%selfenergy%wgrid%type, &
                                 'refreq', &
                                 input%gw%selfenergy%wgrid%size, &
                                 input%gw%selfenergy%wgrid%wmin, &
                                 input%gw%selfenergy%wgrid%wmax)
        else
          call generate_freqgrid(freq_selfc, &
                                 input%gw%freqgrid%fgrid, &
                                 input%gw%freqgrid%fconv, &
                                 input%gw%freqgrid%nomeg, &
                                 input%gw%freqgrid%freqmin, &
                                 input%gw%freqgrid%freqmax)
        end if
        ! call print_freqgrid(freq_selfc,6)

        nw = freq_selfc%nomeg

        if (input%gw%taskname.ne.'g0w0-x') then
          if (allocated(selfec)) deallocate(selfec)
          allocate(selfec(ibgw:nbgw,nw,nkpt))
          selfec(:,:,:) = 0.d0
          if (input%gw%taskname.ne.'cohsex') then
            ! Correlation self-energy at real frequencies after AC procedure
            if (allocated(sigc)) deallocate(sigc)
            allocate(sigc(ibgw:nbgw,nkpt))
            sigc(:,:) = 0.d0
            ! Renormalization (linearization) factors
            if (allocated(znorm)) deallocate(znorm)
            allocate(znorm(ibgw:nbgw,nkpt))
            znorm(:,:) = 0.d0
          else
            ! COHSEX approximation
            if (allocated(sigsx)) deallocate(sigsx)
            allocate(sigsx(ibgw:nbgw,nkpt))
            sigsx(:,:) = 0.d0
            if (allocated(sigch)) deallocate(sigch)
            allocate(sigch(ibgw:nbgw,nkpt))
            sigch(:,:) = 0.d0
          end if ! cohsex
        end if

    end subroutine

    !---------------------------------------------------------------------------
    subroutine delete_selfenergy
      if (allocated(evalks)) deallocate(evalks)
      if (allocated(evalqp)) deallocate(evalqp)
      if (allocated(selfex)) deallocate(selfex)
      if (allocated(selfec)) deallocate(selfec)
      if (allocated(znorm))  deallocate(znorm)
      if (allocated(sigc))   deallocate(sigc)
      if (allocated(sigsx))  deallocate(sigsx)
      if (allocated(sigch))  deallocate(sigch)
      ! frequency grid
      if (allocated(freq_selfc%freqs)) deallocate(freq_selfc%freqs)
      if (allocated(freq_selfc%womeg)) deallocate(freq_selfc%womeg)
    end subroutine

    !---------------------------------------------------------------------------
    subroutine write_selfenergy(ibgw,nbgw,nkpt,nw)
      use modinput
      implicit none
      integer, intent(in) :: ibgw, nbgw
      integer, intent(in) :: nkpt
      integer, intent(in) :: nw
      ! local variables
      integer :: fid, ie, ik, iom
      fid = 777
      ! exchange
      open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      write(fid) ibgw, nbgw, nkpt, selfex
      close(fid)
      ! correlation
      if (input%gw%taskname /= 'g0w0-x') then
        open(fid,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
        write(fid) ibgw, nbgw, nw, nkpt, selfec
        close(fid)
        if (input%gw%taskname == 'cohsex') then
          open(fid,file='COHSEX.OUT',form='UNFORMATTED',status='UNKNOWN')
          write(fid) ibgw, nbgw, nkpt, selfec, sigsx, sigch
          close(fid)
        end if
      end if
    end subroutine

    !> Write exchange self-energy in real text format
    ! TODO(Alex) Would be nicer to print the actual k-point, too
    ! NOTE. Not tested - how does it behave when running with MPI w.r.t. ik?
    subroutine write_exchange_selfenergy(ibgw, nbgw, nkpt)
      use precision, only: dp
      !> Band limits for which GW correction is applied
      integer, intent(in) :: ibgw, nbgw
      !>  Number of k-points 
      integer, intent(in) :: nkpt     
      !> Exchange self-energy   
      !complex(dp), intent(in) :: selfex(:, :)
      !> ile ID unit
      integer :: fid                      
      integer :: ik, ie

      open(newunit=fid, file='SELFX.DAT', form='FORMATTED', status='UNKNOWN')
      write(fid, *) '# first band, last band, N k-points'
      write(fid) ibgw, nbgw, nkpt
      write(fid,*) '# ie,    ik,    selfex'

      do ik = 1, nkpt
          do ie = ibgw, nbgw
            write(fid,'(2i6,2f18.6)') ie, ik, selfex(ie, ik)
          end do
      end do

      close(fid)
    end subroutine

    !> Write correlation self-energy in real text format
    ! TODO(Alex) Would be nicer to print the actual k-point, too
    ! NOTE. Not tested - how does it behave when running with MPI w.r.t. ik?
    subroutine write_correlation_selfenergy(ibgw, nbgw, nw, nkpt)
      use precision, only: dp
      !> Band limits for which GW correction is applied
      integer, intent(in) :: ibgw, nbgw
      !>  Number of frequency points 
      integer, intent(in) :: nw
      !>  Number of k-points 
      integer, intent(in) :: nkpt     
      !> Correlation self-energy   
      !complex(dp), intent(in) :: selfec(:, :, :)
      !> ile ID unit
      integer :: fid                      
      integer :: ik, ie

      open(newunit=fid, file='SELFC.DAT', form='FORMATTED', status='UNKNOWN')
      write(fid, *) '# first band, last band, N k-points, N frequencies'
      write(fid) ibgw, nbgw, nkpt, nw
      write(fid,*) '# ie,   iom,   ik,   selfec'

      do ik = 1, nkpt
        do iom = 1, nw
          do ie = ibgw, nbgw
            write(fid,'(2i6,2f18.6)') ie, iom, ik, selfec(ie, iom, ik)
          end do
        end do
      end do

      close(fid)
    end subroutine


    !---------------------------------------------------------------------------
    subroutine plot_selfc_iw()
      implicit none
      integer(4) :: ik, iw, nk, nb
      character(22) :: frmt
      nb = size(selfec,1)
      nk = size(selfec,3)
      !--------------------------------------
      ! Self-energy along the imaginary axis
      !--------------------------------------
      open(71, file='SelfC-Re-iW.dat', form='FORMATTED', status='UNKNOWN', action='WRITE')
      open(72, file='SelfC-Im-iW.dat', form='FORMATTED', status='UNKNOWN', action='WRITE')
      write(frmt, '("(",i8,"f14.6)")') 1+nb
      do ik = 1, nk
          write(71,*) '# ik = ', ik
          write(72,*) '# ik = ', ik
          do iw = -freq_selfc%nomeg, freq_selfc%nomeg
              if (iw < 0) then
                  write(71,trim(frmt)) -freq_selfc%freqs(abs(iw)), dble(conjg(selfec(:,abs(iw),ik)))
                  write(72,trim(frmt)) -freq_selfc%freqs(abs(iw)), aimag(conjg(selfec(:,abs(iw),ik)))
              else if (iw > 0) then
                  write(71,trim(frmt)) freq_selfc%freqs(iw), dble(selfec(:,iw,ik))
                  write(72,trim(frmt)) freq_selfc%freqs(iw), aimag(selfec(:,iw,ik))
              end if
          end do
          write(71,*); write(72,*)
          write(71,*); write(72,*)
      end do
      close(71)
      close(72)
    end subroutine

    !---------------------------------------------------------------------------
    subroutine plot_selfc()
      implicit none
      integer(4) :: ik, nk, nb, iw
      real(8) :: w
      character(22) :: frmt
      nb = size(selfec,1)
      nk = size(selfec,3)
      !--------------------------------------
      ! Self-energy along the real axis
      !--------------------------------------
      open(71, file='SelfC-Re.dat', form='FORMATTED', status='UNKNOWN', action='WRITE')
      open(72, file='SelfC-Im.dat', form='FORMATTED', status='UNKNOWN', action='WRITE')
      write(frmt, '("(",i8,"f14.6)")') 1+nb
      do ik = 1, nk
          write(71,*) '# ik = ', ik
          write(72,*) '# ik = ', ik
          do iw = 1, freq_selfc%nomeg
            w = freq_selfc%freqs(iw)
            write(71,trim(frmt)) w, dble(selfec(:,iw,ik))
            write(72,trim(frmt)) w, aimag(selfec(:,iw,ik))
        end do
        write(71,*); write(72,*)
        write(71,*); write(72,*)
      end do
      close(71)
      close(72)
    end subroutine

end module
