
subroutine task_analytic_continuation()

    use modinput
    use modmain
    use modgw
    use mod_frequency
    use mod_hdf5
    use mod_mpi_gw
    use m_getunit
    use mod_vxc, only: vxcnn, read_vxcnn
    implicit none
    ! local variables
    integer :: ik, ik_, ie, ie_, fid, recl
    real(8) :: egap
    character(20) :: s1, s2, v(3)

    call init0
    call init1

    nvelgw = chgval-occmax*dble(ibgw-1)
    nbandsgw = nbgw-ibgw+1
    call init_kqpoint_set
    call generate_freqgrid(freq, &
    &                      input%gw%freqgrid%fgrid, &
    &                      input%gw%freqgrid%fconv, &
    &                      input%gw%freqgrid%nomeg, &
    &                      input%gw%freqgrid%freqmin, &
    &                      input%gw%freqgrid%freqmax)

    if (myrank==0) then

      ! allocate the arrays

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
      call bandstructure_analysis('KS', ibgw, nbgw, kset%nkpt, &
                                  evalks(ibgw:nbgw,:), efermi)

      !======================================
      ! Calculate the quasiparticle energies
      !======================================
      call calcevalqp
      call plot_spectral_function()

      !------------------------------------------------------
      ! Write quasi-particle energies to file
      !------------------------------------------------------
      call write_qp_energies('EVALQP.DAT')
      call bandstructure_analysis('G0W0',ibgw,nbgw,kset%nkpt,&
      &                            evalqp(ibgw:nbgw,:),eferqp)

      !----------------------------------------
      ! Save QP energies into binary file
      !----------------------------------------
      call putevalqp('EVALQP.OUT', kset, ibgw, nbgw, evalks, eferks, evalqp, eferqp)

      ! clear memory
      deallocate(evalks, evalfv)
      deallocate(vxcnn)
      call delete_selfenergy

    end if ! myrank

    return
end subroutine
