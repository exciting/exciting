
subroutine task_emac()

    use modinput
    use modmain
    use modgw
    use mod_coulomb_potential
    use mod_mpi_gw
    use m_getunit

    implicit none
    integer :: iq, iom, fid
    complex(8) :: e0, e1, e2

    if (rank==0) call boxmsg(fgw,'=','Calculate the macroscopic dielectric function')

    !==========================
    ! Perform initialization
    !==========================

    ! initialize local GW MPI environment
    call init_mpi_gw()

    ! prepare GW global data
    call init_gw()

    if (.not.input%gw%rpmat) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    else
      write(*,*)
      write(*,*) 'Info(task_emac): Momentum matrix elements read from files.'
      write(*,*)
    end if

    ! occupancy dependent BZ integration weights
    call kintw()

    !==========
    ! Set q=0
    !==========
    iq = 1
    Gamma = .true.

    !========================================
    ! Calculate interstitial basis functions
    !========================================
    matsiz = locmatsiz+Gqset%ngk(1,iq)
    call diagsgi(iq)
    call calcmpwipw(iq)

    !======================================
    ! Calculate the bare coulomb potential
    !======================================
    call calcbarcmb(iq)

    !========================================
    ! Set v-diagonal MB and reduce its size
    !========================================
    call setbarcev(input%gw%barecoul%barcevtol)
    call delete_coulomb_potential()

    !===================================
    ! Calculate the dielectric function
    !===================================

    call init_dielectric_function(mbsiz, 1, freq%nomeg, Gamma)

    select case (trim(input%gw%scrcoul%scrtype))

      case('ppm','PPM')
        call calcepsilon_ppm(iq, 1, freq%nomeg)

      case default
        call calcepsilon(iq, 1, freq%nomeg)
        call calcinveps(1, freq%nomeg)

    end select

    ! clean unused data
    if (allocated(fnm)) deallocate(fnm)
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)

    if (myrank==0) then
      call getunit(fid)
      open(fid, File='EPSMACRO.OUT', Form='Formatted', Action='Write', Status='Replace')
      write(fid,'(a)')'# frequency       eps_{00} (diag)            eps_{00}+LFE (diag)            <eps_{00}^{-1}>'
      do iom = 1, freq%nomeg
        e0 = eps00(1,1,iom)     ! isotropic average without LFE
        e1 = epsh(1,1,iom)+zone ! eps_00^-1
        e2 = 1.d0/e1            ! 1/eps_00^-1
        write(fid,10) freq%freqs(iom), e0, e1, e2
      end do
      close(fid)
      10 format(f12.6,4x,2G12.4,4x,2G12.4,4x,2G12.4)
    end if ! myrank

    call delete_dielectric_function(Gamma)

    return
end subroutine
