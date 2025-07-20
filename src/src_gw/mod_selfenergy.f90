
module mod_selfenergy
    use asserts, only: assert
    use constants, only: zzero
    use gw_io, only: build_file_name, read_from_file, write_to_file, read_bounds_from_file
    use mod_frequency, only: frequency, generate_freqgrid
    use precision, only: i32, dp
    use to_char_conversion, only: to_char

    implicit none

    private 

    type(frequency), public :: freq_selfc

    !--------------!
    ! self-energy  !
    !--------------!

    ! The exchange self-energy
    complex(dp), allocatable, public :: selfex(:,:)

    ! Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
    complex(dp), allocatable, public, target :: mwm(:,:,:)

    ! The correlation self-energy
    complex(dp), allocatable, public :: selfeph(:,:,:)
    complex(dp), allocatable, public :: selfeph0(:,:)
    real(dp),    allocatable, public :: speceph(:,:,:)
    complex(dp), allocatable, public :: selfec(:,:,:)

    ! Correction factors for (q^-1) and (q^-2) singularities
    real(dp), public :: singc1
    real(dp), public :: singc2

    !-------------!
    ! QP Energy   !
    !-------------!

    ! Original KS energies (evalfv will updated via self-consistent cycle)
    real(dp), public :: eferks
    real(dp), allocatable, public :: evalks(:,:)

    ! QP energies
    real(dp), public :: eferqp
    real(dp), allocatable, public :: evalqp(:,:)

    ! Chemical potential alignment
    real(dp), public :: deltaE

    ! Linearization (renormalization) factor
    real(dp),    allocatable, public :: znorm(:,:)

    ! AC to the real axis of the correlation self-energy (selfec)
    complex(dp), allocatable, public :: sigc(:,:)

    ! COHSEX approximation
    complex(dp), allocatable, public :: sigsx(:,:) ! Screened exchange
    complex(dp), allocatable, public :: sigch(:,:) ! Coulomb hole

    !----------------------------------------------------------------------
    ! files to store the self-energy
    !----------------------------------------------------------------------
    character(len=*), parameter, private :: file_name_sigmax = 'SIGMAX_K'
    character(len=*), parameter, private :: file_name_sigmac = 'SIGMAC_K'

    public :: init_selfenergy, plot_selfc, plot_selfc_iw, &
              generate_frequency_grid_for_correlation_self_energy, &
              write_selfec_single_kpoint, write_selfex_single_kpoint, &
              write_selfenergy_binary, delete_selfenergy, &
              read_selfec_from_files, read_selfex_from_files

contains

    !---------------------------------------------------------------------------
    subroutine init_selfenergy(ibgw,nbgw,nkpt)
        use modinput, only: input
        integer(i32), intent(in) :: ibgw, nbgw
        integer(i32), intent(in) :: nkpt
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
      use modinput, only: input
      integer(i32), intent(in) :: ibgw, nbgw
      integer(i32), intent(in) :: nkpt
      integer(i32), intent(in) :: nw
      ! local variables
      integer(i32) :: fid, ie, ik, iom
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
      integer(i32), intent(in) :: ibgw, nbgw
      !>  Number of k-points 
      integer(i32), intent(in) :: nkpt
      !> ile ID unit
      integer(i32) :: fid                      
      integer(i32) :: ik, ie

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
      integer(i32), intent(in) :: ibgw, nbgw
      !>  Number of frequency points 
      integer(i32), intent(in) :: nw
      !>  Number of k-points 
      integer(i32), intent(in) :: nkpt
      
      integer(i32) :: fid                      
      integer(i32) :: ik, ie, iom

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


    !> Read the exchange part of the self-energy from files
    subroutine read_selfec_from_files( kpt_indexes, file_format )
      !> List of k-point indexes
      integer(i32), intent(in) :: kpt_indexes(:)
      !> Format of the file where to print. It can be e.g. 'text' or 'binary'
      character(len=*), intent(in) :: file_format

      integer(i32) :: i, lbounds(2), ubounds(2)
      integer(i32), parameter :: maxlen = 30
      character(len=maxlen) :: file_name
      
      call build_file_name( file_name_sigmac, kpt_indexes(1), file_name )
      call read_bounds_from_file( file_name, file_format, lbounds, ubounds )
      if( allocated(selfec) ) deallocate(selfec)
      allocate( selfec(lbounds(1):ubounds(1), lbounds(2):ubounds(2), 1:size(kpt_indexes)) )
      do i = 1, size( kpt_indexes )
        call build_file_name( file_name_sigmac, kpt_indexes(i), file_name )
        call read_from_file( file_name, selfec(:, :, i), lbound(selfec), file_format )
      end do

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


    !> Read the exchange part of the self-energy from files
    subroutine read_selfex_from_files( kpt_indexes, file_format )
      !> List of k-point indexes
      integer(i32), intent(in) :: kpt_indexes(:)
      !> Format of the file where to print. It can be e.g. 'text' or 'binary'
      character(len=*), intent(in) :: file_format

      integer(i32) :: i, l_bound(1), u_bound(1)
      integer(i32), parameter :: maxlen = 30
      character(len=maxlen) :: file_name 
      
      call build_file_name( file_name_sigmax, kpt_indexes(1), file_name )
      call read_bounds_from_file( file_name, file_format, l_bound, u_bound )
      if( allocated(selfex) ) deallocate(selfex)
      allocate( selfex(l_bound(1):u_bound(1), 1:size(kpt_indexes)) )
      do i = 1, size( kpt_indexes )
        call build_file_name( file_name_sigmax, kpt_indexes(i), file_name )
        call read_from_file( file_name, selfex(:, i), l_bound(1), file_format )
      end do

    end subroutine


    !---------------------------------------------------------------------------
    subroutine plot_selfc_iw()
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
    subroutine plot_selfc(frequencies, list_kpt_idx, sigmac, first_band)
      !> Array with the frequencies (in most cases, assumed to be along the real axis)
      real(dp), contiguous, intent(in) :: frequencies(:)
      !> Array containing the indexes of k-points
      integer(i32), contiguous, intent(in) :: list_kpt_idx(:)
      !> First band used to compute `sigmac`
      integer(i32), intent(in) :: first_band
      !> Correlation part of the self-energy \(\Sigma_C\)
      complex(dp), contiguous, intent(in) :: sigmac(first_band:, :, :)
      
      character(len=*), parameter :: file_name_selfc = 'SELFENERGY_C_K'
      character(len=*), parameter :: default_extension ='.OUT'
      integer(i32) :: i, ik, i_band, i_freq, i_unit

      call assert( size(frequencies) == size(sigmac, 2), "frequencies and sigmac have incompatible sizes")
      call assert( size(list_kpt_idx) == size(sigmac, 3), "list_kpt_idx and sigmac have incompatible sizes")
      do i = 1, size(sigmac, 3)
        ik = list_kpt_idx(i)
        open( newunit=i_unit, file=file_name_selfc//to_char(ik)//default_extension, action="write" )
        do i_band = first_band, size(sigmac, 1)
          do i_freq = 1, size(frequencies)
            write( i_unit, '(2I5, F10.6, 2F16.10)' ), i_band, i_freq, frequencies(i_freq), selfec(i_band, i_freq, i)
          end do
        end do
        close( i_unit )
      end do
    end subroutine

end module
