!> Postprocess and output optical properties files
!> Reads in EPSILON_11.OUT and EPSILON_12.OUT, and outputs KERR.OUT
subroutine moke
    use modmain
    use modinput
    use modmpi
    use unit_conversion, only: hartree_to_ev, radians_to_degrees
    use constants, only: pi
    implicit none

    ! Declarations

    integer :: iw, iostat
    !> Unit conversion factor
    real(8) :: conversion
    complex(8) :: exx, exy, zt1
    !> Grid size
    integer :: grid_size
    !> Energy string
    character(len=2) :: energy_str

    ! Main Routine

    real(8), allocatable :: w(:)
    real(8), allocatable :: eps11(:,:)
    real(8), allocatable :: eps12(:,:)
    complex(8), allocatable :: kerr(:)

    ! initialization for DIELMAT
    if (.not.associated(input%properties%dielmat)) then
        input%properties%dielmat => getstructdielmat(emptynode)
    end if
    
    input%properties%dielmat%intraband = input%properties%moke%intraband
    grid_size = input%properties%moke%wgrid
    input%properties%dielmat%wgrid = grid_size
    input%properties%dielmat%wmax = input%properties%moke%wmax
    input%properties%dielmat%scissor = input%properties%moke%scissor
    input%properties%dielmat%swidth = input%properties%moke%swidth
    input%properties%dielmat%drude = input%properties%moke%drude
    input%properties%dielmat%tevout = .false.

    ! for MOKE only two optical components required
    deallocate(input%properties%dielmat%epscomp)
    allocate(input%properties%dielmat%epscomp(2,2))
    !   XX
    input%properties%dielmat%epscomp(1,1) = 1
    input%properties%dielmat%epscomp(2,1) = 1
    !   XY
    input%properties%dielmat%epscomp(1,2) = 1
    input%properties%dielmat%epscomp(2,2) = 2
    
    call dielmat

    call barrier
    if (rank==0) then

        allocate(w(grid_size))
        allocate(eps11(grid_size,2), eps12(grid_size,2))
        allocate(kerr(grid_size))

        ! Read diagonal contribution to optical conductivity
        open(50, File='EPSILON_11.OUT', Action='READ', Status='OLD', Form='FORMATTED', IoStat=IoStat)
        if (iostat /= 0) then
            write (*,*)
            write (*, '("Error(moke): error opening EPSILON_11.OUT")')
            write (*,*)
            stop
        end if

        ! Skip header
        read(50, *)

        do iw = 1, grid_size
            read(50, *) w(iw), eps11(iw, 1), eps11(iw, 2)
        end do
        close(50)

        ! Read off-diagonal contribution to optical conductivity
        open (50, File='EPSILON_12.OUT', Action='READ', Status='OLD', Form='FORMATTED', IoStat=IoStat)
        if (iostat /= 0) then
            write (*,*)
            write (*, '("Error(moke): error opening EPSILON_12.OUT")')
            write (*,*)
            stop
        end if

        ! Skip header
        read(50, *)

        do iw = 1, grid_size
            read(50, *) w(iw), eps12(iw, 1), eps12(iw, 2)
        end do
        close(50)

        ! calculate the complex Kerr angle
        do iw = 1, grid_size
            exx = cmplx(eps11(iw,1),eps11(iw,2),8)
            exy = cmplx(eps12(iw,1),eps12(iw,2),8)
            zt1 = (exx-zone)*sqrt(exx)
            if (abs(zt1) > 1.d-8) then
                kerr(iw) = -exy / zt1
            else
                kerr(iw) = zzero
            end if
        end do

        ! Output results
        open(50, File='KERR.OUT', Action='WRITE', Form='FORMATTED')

        if (input%properties%moke%tevout) then
            conversion = hartree_to_ev
            energy_str = 'eV'
        else
            conversion = 1.0d0
            energy_str = 'Ha'
        end if

        write(50, '(A, A, A)') '# Energy (', energy_str, '), Re{Kerr} (degrees),  Im{Kerr} (degrees)'
        do iw = 1, grid_size
            write(50, '(3G18.10)') conversion * w(iw), kerr(iw) * radians_to_degrees
        end do

        close(50)

        write(*,*)
        write(*, '("Info(moke):")')
        write(*, '(" complex Kerr angle in degrees written to KERR.OUT")')
        write(*, '(A, A)') " Output energy is in ", energy_str
        write(*,*)

    end if ! rank

end subroutine
