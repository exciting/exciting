subroutine task_band_specfunc()

    use modinput
    use modmain
    use modgw
    use mod_vxc, only: vxcnn, read_vxcnn
    use mod_frequency

    implicit none
    integer(4) :: ik, ib, ib0, ik_path
    character(80) :: fname, s
    logical :: exist
    integer(4) :: fid, recl
    integer(4) :: i, j, iw
    real(8) :: egap
    real(8), allocatable :: sf(:)
    real(8) :: w, sRe, sIm, div
    complex(8) :: sxc
    type(k_set) :: ksetnr

    !--------------------------------
    ! exciting global initialization
    !--------------------------------
    call init0()
    call init1()

    nvelgw = chgval-occmax*dble(ibgw-1)
    nbandsgw = nbgw-ibgw+1
    call init_kqpoint_set
    call generate_freqgrid(freq, &
                           input%gw%freqgrid%fgrid, &
                           input%gw%freqgrid%fconv, &
                           input%gw%freqgrid%nomeg, &
                           input%gw%freqgrid%freqmin, &
                           input%gw%freqgrid%freqmax)

    allocate(vxcnn(ibgw:nbgw,kset%nkpt))
    call init_selfenergy(ibgw,nbgw,kset%nkpt)

    ! real frequency grid
    if ( .not.associated(input%gw%selfenergy%wgrid) ) &
        input%gw%selfenergy%wgrid => getstructwgrid(emptynode)
    call delete_freqgrid(freq_selfc)
    call generate_freqgrid(freq_selfc, &
                        input%gw%selfenergy%wgrid%type, &
                        'refreq', &
                        input%gw%selfenergy%wgrid%size, &
                        input%gw%selfenergy%wgrid%wmin, &
                        input%gw%selfenergy%wgrid%wmax)
    deallocate(selfec)
    allocate(selfec(ibgw:nbgw,freq_selfc%nomeg,kset%nkpt))

    ! read data from files
    call readevalqp('EVALQP.OUT', kset, ibgw, nbgw, evalks, eferks, evalqp, eferqp)
    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(ibgw:nbgw,kset%nkpt))
    evalfv(:,:) = evalks(:,:)
    call read_vxcnn()
    call readselfx()
    call readselfc()

    ! KS states analysis
    call fermi_exciting(.false., &
    &                   nvelgw, &
    &                   nbandsgw, kset%nkpt, evalks(ibgw:nbgw,:), &
    &                   kset%ntet, kset%tnodes, kset%wtet,kset%tvol, &
    &                   efermi, egap, fermidos)
    call bandstructure_analysis('KS', &
                                ibgw, nbgw, kset%nkpt, &
                                evalks(ibgw:nbgw,:), efermi)

    ! non-reduced k-points
    call generate_k_vectors(ksetnr, &
                            bvec, &
                            input%gw%ngridq, &
                            input%gw%vqloff, &
                            .false.)

    !------------------------
    ! Read KS bandstructure
    !------------------------
    fname = 'bandstructure.dat'
    inquire(File=fname, Exist=exist)
    if (.not.exist) then
        write(*,*) 'ERROR(task_band): bandstructure.dat file is not found!'
        write(*,*) '    Run properties/bandstructure first to produce KS spectrum.'
        stop
    end if

    open(70, File='bandstructure.dat', Action='Read', Status='Old')
    read(70,*) s, ib0, nstfv, nkpt
    if (allocated(vkl)) deallocate(vkl)
    allocate(vkl(3,nkpt))
    if (allocated(dpp1d)) deallocate(dpp1d)
    allocate(dpp1d(nkpt))
    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,nkpt))
    do ib = ib0, nstfv
        do ik = 1, nkpt
            read(70,*) i, j, vkl(:,ik), dpp1d(ik), evalfv(ib,ik)
        end do
        read(70,*) ! skip line
    end do
    close(70)

    open(70,file='sf-band-heatMap.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    allocate(sf(ibgw:nbgw))
    do ik_path = 1, nkpt
        ! determine the nearest neighbor point in the selfenergy grid
        ik = find_nearest_grid_point(ik_path)
        do iw = 1, freq_selfc%nomeg
            w = freq_selfc%freqs(iw)
            ! compute spectral function
            do ib = ibgw, nbgw
                sxc = selfex(ib,ik) + selfec(ib,iw,ik) - vxcnn(ib,ik)
                sRe = dble(sxc)
                sIm = aimag(sxc) + input%gw%selfenergy%swidth
                div = (w-evalfv(ib,ik_path)-sRe)**2 + sIm**2
                sf(ib) = 1.d0/pi * abs(sIm) / div
            end do
            write(70,'(3f14.6)') dpp1d(ik_path), w, sum(sf(ibgw:nbgw))
        end do
        write(70,*)
    end do
    deallocate(sf)
    close(70)

    call delete_k_vectors(ksetnr)

contains

    function find_nearest_grid_point(ik_in) result(ik_out)
        implicit none
        integer(4), intent(in) :: ik_in
        integer(4) :: ik_out
        integer(4) :: ikr, ik_min
        real(8) :: vl_in(3), vc_in(3)
        real(8) :: vl(3), vc(3)
        integer(4) :: iv(3)
        real(8) :: dist, dist_min

        vl_in = vkl(:,ik_in)
        call r3frac( 1.d-6, vl_in, iv) ! translate to the first BZ
        call r3mv(bvec, vl_in, vc_in)  ! convert to cartesian

        dist_min = 100.d0
        ik_min = 0
        do ikr = 1, ksetnr%nkpt
            vc = ksetnr%vkc(:,ikr)
            dist = (vc(1)-vc_in(1))**2 + (vc(2)-vc_in(2))**2 + (vc(3)-vc_in(3))**2
            if (dist < dist_min) then
                ik_min = ikr
                dist_min = dist
            end if
        end do
        if (ikr == 0) stop 'Error(task_band_specfunc::find_nearest_grid_point) Search algorithm failed!'
        ik_out = kset%ik2ikp(ik_min)
        return
    end function

end subroutine