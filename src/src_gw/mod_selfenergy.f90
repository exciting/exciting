
module mod_selfenergy
    use constants, only: zzero
    use gw_io, only: build_file_name, read_from_file, write_to_file
    use mod_frequency
    use precision, only: i32, dp

    implicit none

    type(frequency) :: freq_selfc

    !--------------!
    ! self-energy  !
    !--------------!

    ! The exchange self-energy
    complex(dp), allocatable :: selfex(:,:)

    ! Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
    complex(dp), allocatable :: mwm(:,:,:)
    target mwm

    ! The correlation self-energy
    complex(dp), allocatable :: selfeph(:,:,:)
    complex(dp), allocatable :: selfeph0(:,:)
    real(dp),    allocatable :: speceph(:,:,:)
    complex(dp), allocatable :: selfec(:,:,:)

    ! Correction factors for (q^-1) and (q^-2) singularities
    real(dp) :: singc1
    real(dp) :: singc2

    !-------------!
    ! QP Energy   !
    !-------------!

    ! Original KS energies (evalfv will updated via self-consistent cycle)
    real(dp) :: eferks
    real(dp), allocatable :: evalks(:,:)

    ! QP energies
    real(dp) :: eferqp
    real(dp), allocatable :: evalqp(:,:)

    ! Chemical potential alignment
    real(dp) :: deltaE

    ! Linearization (renormalization) factor
    real(dp),    allocatable :: znorm(:,:)

    ! AC to the real axis of the correlation self-energy (selfec)
    complex(dp), allocatable :: sigc(:,:)

    ! COHSEX approximation
    complex(dp), allocatable :: sigsx(:,:) ! Screened exchange
    complex(dp), allocatable :: sigch(:,:) ! Coulomb hole

    !----------------------------------------------------------------------
    ! files to store the self-energy
    !----------------------------------------------------------------------
    character(len=*), parameter, private :: file_name_sigmax = 'SIGMAX_K'
    character(len=*), parameter, private :: file_name_sigmac = 'SIGMAC_K'

contains

    !---------------------------------------------------------------------------
    subroutine init_selfenergy(ibgw,nbgw,nkpt)
        use modinput, only: input
        implicit none
        integer, intent(in) :: ibgw, nbgw
        integer, intent(in) :: nkpt
        ! local
        integer(i32) :: nw

        ! KS eigenvalues
        if (allocated(evalks)) deallocate(evalks)
        allocate(evalks(ibgw:nbgw,nkpt), source=0.0_dp)

        ! Quasi-Particle energy
        if (allocated(evalqp)) deallocate(evalqp)
        allocate(evalqp(ibgw:nbgw,nkpt), source=0.0_dp)

        ! Exchange self-energy
        if (allocated(selfex)) deallocate(selfex)
        allocate(selfex(ibgw:nbgw,nkpt), source=zzero)

        ! Correlation self-energy
        call generate_frequency_grid_for_correlation_self_energy( input%gw )
        nw = freq_selfc%nomeg

        if (input%gw%taskname.ne.'g0w0-x') then
          if (allocated(selfec)) deallocate(selfec)
          allocate(selfec(ibgw:nbgw,nw,nkpt), source=zzero)

          if (input%gw%taskname.ne.'cohsex') then
            ! Correlation self-energy at real frequencies after AC procedure
            if (allocated(sigc)) deallocate(sigc)
            allocate(sigc(ibgw:nbgw,nkpt), source=zzero)
            
            ! Renormalization (linearization) factors
            if (allocated(znorm)) deallocate(znorm)
            allocate(znorm(ibgw:nbgw,nkpt), source=0.0_dp)
            
          else
            ! COHSEX approximation
            if (allocated(sigsx)) deallocate(sigsx)
            allocate(sigsx(ibgw:nbgw,nkpt), source=zzero)
            
            if (allocated(sigch)) deallocate(sigch)
            allocate(sigch(ibgw:nbgw,nkpt), source=zzero)
            
          end if ! cohsex
        end if

    end subroutine

    !> Generate the frequency grid needed by the self energy
    subroutine generate_frequency_grid_for_correlation_self_energy( gw_inp )
      use modinput, only: gw_type, emptynode, getstructwgrid
      !> GW input parameters
      type(gw_type), intent(inout) :: gw_inp

      if ( gw_inp%selfenergy%method == 'cd') then
        if ( .not.associated(gw_inp%selfenergy%wgrid) ) &
            gw_inp%selfenergy%wgrid => getstructwgrid(emptynode)
        call generate_freqgrid(freq_selfc, &
                               gw_inp%selfenergy%wgrid%type, &
                               'refreq', &
                               gw_inp%selfenergy%wgrid%size, &
                               gw_inp%selfenergy%wgrid%wmin, &
                               gw_inp%selfenergy%wgrid%wmax)
      else
        call generate_freqgrid(freq_selfc, &
                               gw_inp%freqgrid%fgrid, &
                               gw_inp%freqgrid%fconv, &
                               gw_inp%freqgrid%nomeg, &
                               gw_inp%freqgrid%freqmin, &
                               gw_inp%freqgrid%freqmax)
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
    subroutine write_selfenergy_binary(ibgw,nbgw,nkpt,nw)
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
      integer :: ik, ie, iom

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

    !> Write the correlation part of the self-energy for a given k-point
    subroutine write_selfec_single_kpoint( ik, file_format )
      !> Index of the current k-point
      integer(i32), intent(in) :: ik
      !> Format of the file where to print. It can be e.g. 'text' or 'binary'
      character(len=*), intent(in) :: file_format

      integer(i32), parameter :: maxlen = 30
      character(len=maxlen) :: file_name 
      
      call build_file_name( file_name_sigmac, ik, file_name )
      call write_to_file( file_name, selfec(:, :, ik), [ lbound( selfec, 1 ), lbound( selfec, 2 ) ], file_format )

    end subroutine

    !> Write the exchange part of the self-energy for a given k-point
    subroutine write_selfex_single_kpoint( ik, file_format )
      !> Index of the current k-point
      integer(i32), intent(in) :: ik
      !> Format of the file where to print. It can be e.g. 'text' or 'binary'
      character(len=*), intent(in) :: file_format

      integer(i32), parameter :: maxlen = 30
      character(len=maxlen) :: file_name 
      
      call build_file_name( file_name_sigmax, ik, file_name )
      call write_to_file( file_name, selfex(:, ik), lbound( selfex, 1 ), file_format )

    end subroutine


    !---------------------------------------------------------------------------
    subroutine plot_selfc_iw()
      implicit none
      integer(i32) :: ik, iw, nk, nb
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
      integer(i32) :: ik, nk, nb, iw
      real(dp) :: w
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
