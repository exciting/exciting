
subroutine task_evalqp()
    use mod_bands, only: evalfv, bandstructure_analysis
    use mod_frequency
    use mod_hdf5
    use mod_vxc, only: vxcnn, read_vxcnn, deallocate_vxcnn
    use modinput
    use modgw
    use modmain
    use modmpi, only: rank
    use quasiparticle_energies, only: write_qp_energies_text_format
    use precision, only: dp
    
    implicit none
    ! local variables
    integer :: ikp, ik, ie
    real(8) :: egap
    character(20) :: s1, s2, v(3)
    logical :: reducek

    input%groundstate%stypenumber = -1
    call init0
    reducek = input%groundstate%reducek
    input%groundstate%reducek = .false.
    call init1()
    input%groundstate%reducek = reducek

    nvelgw = chgval-occmax*dble(ibgw-1)
    nbandsgw = nbgw-ibgw+1
    call init_kqpoint_set
    call generate_freqgrid(freq, &
    &                      input%gw%freqgrid%fgrid, &
    &                      input%gw%freqgrid%fconv, &
    &                      input%gw%freqgrid%nomeg, &
    &                      input%gw%freqgrid%freqmin, &
    &                      input%gw%freqgrid%freqmax)

    if (rank==0) then

      ! allocate the arrays
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
      if (allocated(evalks)) deallocate(evalks)
      allocate(evalks(nstfv,kset%nkpt))
      filext = "_GW.OUT"
      do ikp = 1, kset%nkpt
        ik = kset%ikp2ik(ikp)
        call getevalfv(kqset%vkl(:,ik), evalks(:,ikp))
      end do

      if (allocated(evalfv)) deallocate(evalfv)
      allocate(evalfv(ibgw:nbgw,kset%nkpt))
      evalfv(ibgw:nbgw,:) = evalks(ibgw:nbgw,:)
      call read_vxcnn('binary')
      call readselfx()
      call readselfc()

      ! KS states analysis
      call fermi_exciting(.false., &
      &                   nvelgw, &
      &                   nbandsgw, kset%nkpt, evalks(ibgw:nbgw,:), &
      &                   kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
      &                   efermi, egap, fermidos)
      call bandstructure_analysis('KS', ibgw, evalks(ibgw:nbgw,:), efermi, .true.)

      !======================================
      ! Calculate the quasiparticle energies
      !======================================
      call calcevalqp
      if (input%gw%printSelfC)            call plot_selfc(freq_selfc%freqs, [(ik, ik=1,kset%nkpt)], selfec, first_band=1)
      if (input%gw%printSpectralFunction) call plot_spectral_function()

      !------------------------------------------------------
      ! Write quasi-particle energies to file
      !------------------------------------------------------
      call write_qp_energies_text_format( [(ik, ik=1,kset%nkpt)], kset%vkl, kset%wkpt, &
        ibgw, evalks, evalqp, real( vxcnn%diag_elements(ibgw:, :), dp ), selfex, sigc, znorm )
      call bandstructure_analysis('G0W0',ibgw,evalqp(ibgw:nbgw,:),eferqp, .true.)

      !----------------------------------------
      ! Save QP energies into binary file
      !----------------------------------------
      call putevalqp('EVALQP.OUT', kset, ibgw, nbgw, evalfv - efermi, 0.0, evalqp, eferqp)

      ! clear memory
      deallocate(evalks, evalfv)
      call delete_selfenergy
      call deallocate_vxcnn

    end if ! rank

    return
end subroutine
