!
! Description:
!
subroutine ldos()

    use modinput
    use modmain,        only: ngrid, ngrtot, ngkmax, apwordmax, lmmaxapw, nmatmax, nrmtmax, nspinor, &
                              nstfv, nstsv, vkl, vgkl, evalsv, ngk, gkc, tpgkc, sfacgk, omega, nkpt, &
                              ikmap, efermi, twopi, avec, task
    use mod_rgrid
    use mod_xsf_format
    use mod_cube_format
    use modmpi

    implicit none
    integer(4) :: ik, ib
    character(80) :: fname

    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: wfmt(:,:,:,:,:)
    complex(8), allocatable :: wfir(:,:,:)
    complex(8), allocatable :: zdata(:)

    integer(4) :: i1, i2, i3, ip, np, nptot, npstep, i_dim
    integer(4) :: nw , iw, nsk(3), iw_f, iw_vb, iw_cb
    real(8) :: v0, dw, vbz
    real(8) :: abc(3), z, delta
    real(8), allocatable :: rho(:,:,:)
    real(8), allocatable :: weight(:,:,:)
    real(8), allocatable :: w(:)
    real(8), allocatable :: g(:)
    real(8), allocatable :: e(:,:)

    type(rgrid) :: grid
    type(plot3d_type), pointer :: plot3d

    ! Initialize global exciting variables
    call init0
    call init1
    call readstate
    call readfermi
    call linengy
    call genapwfr
    call genlofr(.False.)

    ! lengths of the unit cell basis vectors
    do i1 = 1, 3
        abc(i1) = sqrt(avec(1,i1)**2+avec(2,i1)**2+avec(3,i1)**2)
    end do
    ! write(*,'(a,3f16.4)') "a, b, c =", abc

    ! allocate local arrays
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfv(nmatmax,nstfv))
    allocate(evecsv(nstsv,nstsv))
    allocate(wfmt(lmmaxapw,nrmtmax,natmtot,nspinor,nstsv))
    allocate(wfir(ngrtot,nspinor,nstsv))

    ! plot3d constructor
    allocate(plot3d)
    allocate(plot3d%box)
    allocate(plot3d%box%origin)
    allocate(plot3d%box%pointarray(3))
    allocate(plot3d%box%pointarray(1)%point)
    allocate(plot3d%box%pointarray(2)%point)
    allocate(plot3d%box%pointarray(3)%point)

    plot3d%box%grid(:) = input%properties%ldos%grid(:)
    plot3d%box%origin%coord = (/0.d0, 0.d0, 0.d0/)
    plot3d%box%pointarray(1)%point%coord = (/1.d0, 0.d0, 0.d0/)
    plot3d%box%pointarray(2)%point%coord = (/0.d0, 1.d0, 0.d0/)
    plot3d%box%pointarray(3)%point%coord = (/0.d0, 0.d0, 1.d0/)

    ! rgrid constructor
    grid = gen_3d_rgrid(plot3d, 1)
    ! call print_rgrid(grid)

    allocate(zdata(grid%npt))
    allocate(rho(grid%ngrid(1), grid%ngrid(2), grid%ngrid(3)))

    ! state weights / local charge (?)
    i_dim = 3 ! z-direction
    nptot = grid%ngrid(i_dim)
    delta = input%properties%ldos%delta
    if (delta > abc(3)) stop 'Error(ldos): Integration volume specified by delta is too big!'
    np = min(nint(abc(i_dim)/delta), nptot)
    allocate(weight(np, nstsv, nkpt))
    weight(:,:,:) = 0.d0
    npstep = max(nptot/np, 1)

    ! integration volume
    v0 = omega / dble(grid%npt)

    ! print*, 'delta=', delta
    ! print*, 'nptot=', nptot
    ! print*, 'np=', np
    ! print*, 'npstep=', npstep
    ! print*, 'v0=', v0

    !---------------------------------------
    ! begin parallel loop over k-points
    !---------------------------------------
    splittfile = .False.
#ifdef MPI
    call barrier()
    do ik = firstofset(rank,nkpt), lastofset(rank,nkpt)
#else
    do ik = 1, nkpt
#endif
        ! get the eigenvectors and values from file
        call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfv)
        call getevecsv(vkl(:,ik), evecsv)

        ! compute the matching coefficients
        call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), &
        &          sfacgk(:,:,1,ik), apwalm)

        do ib = 1, nstsv

            call genwfsv_new(ik, ib, ib, apwalm, evecfv, evecsv, wfmt, wfir)
            call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ib), wfir(:,1,ib), zdata)
            ! write(*,*) 'zdata=', ik, ib, sum(zdata)
            ! if (rank==0) then
            !     write(fname,'("wf3d-",i4.4,"-",i4.4,".xsf")') ik, ib
            !     call str_strip(fname)
            !     call write_structure_xsf(fname)
            !     call write_3d_xsf(fname, 'module squared', grid%boxl(1:4,:), grid%ngrid, grid%npt, abs(zdata)**2)
            ! end if

            ! real space partial density
            ip = 0
            do i3 = 1, grid%ngrid(3)
            do i2 = 1, grid%ngrid(2)
            do i1 = 1, grid%ngrid(1)
                ip = ip+1
                rho(i1,i2,i3) = abs(zdata(ip))**2
            end do
            end do
            end do
            ! write(*,*) 'rho=', ik, ib, sum(rho)

            ! integral over the volume slice
            ip = 1
            do i3 = 1, nptot
                weight(ip,ib,ik) = weight(ip,ib,ik) + v0*sum(rho(:,:,i3))
                if (mod(i3,npstep)==0) ip = ip+1
            end do

        end do ! ib

    end do ! ik
#ifdef MPI
    call mpi_allgatherv_ifc(nkpt, inplace=.False., rlen=np*nstsv, rbuf=weight)
    call barrier()
#endif

    ! memory cleaning
    deallocate(zdata)
    deallocate(rho)
    deallocate(apwalm)
    deallocate(evecfv)
    deallocate(evecsv)
    deallocate(wfmt, wfir)

    !------------------------------------------------
    ! finally compute the local density of states
    !------------------------------------------------
    if (rank == 0) then

        allocate(e(nstsv,nkpt))
        do ik = 1, nkpt
            call getevalsv(vkl(:,ik), evalsv(:, ik))
            do ib = 1, nstsv
                ! subtract the Fermi energy
                e(ib,ik) = evalsv(ib,ik)-efermi
                ! correction for scissors operator
                if (e(ib,ik) > 0.d0) &
                    e(ib,ik) = e(ib,ik) + input%properties%ldos%scissor
            end do
        end do

        ! generate energy grid
        nw = input%properties%ldos%nwdos
        allocate(w(nw))
        dw = (input%properties%ldos%winddos(2)-input%properties%ldos%winddos(1)) / dble(nw-1)
        iw_f = 0
        do iw = 1, nw
            w(iw) = dw*dble(iw-1) + input%properties%ldos%winddos(1)
            if (w(iw) <= 0.d0) iw_f = iw
        end do
        ! write(*,*) 'iw_f=', iw_f, w(iw_f)
        if (iw_f == 0) then
            write(*,*)
            write(*,*) "Warning(ldos) Fermi energy is outside of the chosen frequency window!"
            write(*,*) "              Impossible to determine locations of VBM and CBm!"
            write(*,*)
        end if

        ! number of subdivisions used for interpolation
        nsk(:) = max(input%properties%ldos%ngrdos/input%groundstate%ngridk(:), 1)

        open(30, file="ldos.out", action="write")
        open(31, file="band_edges.out", action="write")
        allocate(g(nw))
        do ip = 1, np
            call brzint_new(input%properties%ldos%nsmdos, &
                            input%groundstate%ngridk, nsk, ikmap, &
                            nw, input%properties%ldos%winddos, &
                            nstsv, nstsv, &
                            e, weight(ip,:,:), g)
            do iw = 1, nw
                write(30,'(2f16.6)') w(iw), g(iw)
            end do
            write(30,*); write(30,*)
            !------------------------
            ! search for VBM and CBm
            !------------------------
            call find_vbm_cbm_1(iw_f, w, g, iw_vb, iw_cb)
            ! call find_vbm_cbm_2(iw_f, w, g, iw_vb, iw_cb)
            ! output
            z = dble(ip-1)*abc(i_dim)/dble(np)
            write(31,'(3f16.6)') z, w(iw_vb), w(iw_cb)
        end do
        close(30)
        close(31)

        deallocate(g)
        deallocate(w)
        deallocate(e)
        call delete_rgrid(grid)

    end if ! rank==0

    deallocate(weight)

contains

    subroutine find_vbm_cbm_1(iw_f, w, g, iw_vb, iw_cb)
        implicit none
        integer(4), intent(in) :: iw_f
        real(8), intent(in) :: w(:)
        real(8), intent(in) :: g(:)
        integer(4), intent(out) :: iw_vb, iw_cb
        do iw_vb = iw_f, 1, -1
            if (g(iw_vb) > input%properties%ldos%tol) exit
        end do
        do iw_cb = iw_f, nw
            if (g(iw_cb) > input%properties%ldos%tol) exit
        end do
    end subroutine

    subroutine find_vbm_cbm_2(iw_f, w, g, iw_vb, iw_cb)
        implicit none
        integer(4), intent(in) :: iw_f
        real(8), intent(in) :: w(:)
        real(8), intent(in) :: g(:)
        integer(4), intent(out) :: iw_vb, iw_cb
        integer(4) :: i
        real(8) :: f, ftol, s
        ! VB search
        f = integrate(1, iw_f, w, g)
        ftol = input%properties%ldos%tol * f
        s = 0.d0
        do i = iw_f, 1, -1
            s = integrate(i, iw_f, w, g)
            if (s >= ftol) exit
        end do
        iw_vb = i
        ! CB search
        f = integrate(iw_f, size(w), w, g)
        ftol = input%properties%ldos%tol * f
        s = 0.d0
        do i = iw_f, size(w)
            s = integrate(iw_f, i, w, g)
            if (s >= ftol) exit
        end do
        iw_cb = i
    end subroutine

    function integrate(i1, i2, x, f) result(s)
        implicit none
        integer(4), intent(in) :: i1, i2
        real(8), intent(in) :: x(:)
        real(8), intent(in) :: f(:)
        real(8) :: s
        integer(4) :: i
        s = 0.d0
        do i = min(i1+1,i2+1), max(i1,i2)
            s = s + 0.5d0*(f(i)+f(i-1))*(w(i)-w(i-1))
        end do
        return
    end function

end subroutine