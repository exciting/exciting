subroutine wannierlauncher
    use modinput
    use mod_wannier
    use mod_wfutil

    implicit none

    if( associated( input%properties%wannier)) then
      ! generate Wannier functions
      call wannier_init
      if( input%properties%wannier%do .eq. "fromscratch") then
        if( input%properties%wannier%printproj) call wannier_writeinfo_lo
        call wannier_writeinfo_overall
        call printbox( wf_info, '*', "Wannierization")
        do wf_group = 1, wf_ngroups
          call wannier_writeinfo_task
          if( wf_groups( wf_group)%method .eq. "pro") then
            call wannier_gen_pro
          !else if( input%properties%wannier%method .eq. "prowan") then
          !  call wannier_projonwan
          !  call wannier_gen_pro
          !  call wannier_writetransform
          else if( wf_groups( wf_group)%method .eq. "opf") then
            call wannier_gen_opf
          !else if( input%properties%wannier%method .eq. "opfwan") then
          !  call wannier_projonwan
          !  call wannier_gen_opf
          !  call wannier_writetransform
          else if( wf_groups( wf_group)%method .eq. "promax") then
            call wannier_gen_pro
            call wannier_maxloc
          !else if( input%properties%wannier%method .eq. "prowanmax") then
          !  call wannier_projonwan
          !  call wannier_gen_pro
          !  call wannier_maxloc
          !  call wannier_writetransform
          else if( wf_groups( wf_group)%method .eq. "opfmax") then
            call wannier_gen_opf
            call wannier_maxloc
          else if( wf_groups( wf_group)%method .eq. "scdm") then
            call wannier_scdm
          else if( wf_groups( wf_group)%method .eq. "scdmmax") then
            call wannier_scdm
            call wannier_maxloc
          else if( wf_groups( wf_group)%method .eq. "disentangle") then
            call wannier_subspace
            call wannier_gen_opf
            call wannier_maxloc
          else
            write(*,*) " Error (propertylauncher): invalid value for attribute method"
          end if
        end do
        call wannier_writetransform
      
      else if( input%properties%wannier%do .eq. "fromfile") then
        call wannier_gen_fromfile
        if( input%properties%wannier%printproj) call wannier_writeinfo_lo

      else
        call wannier_gen_fromfile
        if( input%properties%wannier%printproj) call wannier_writeinfo_lo
        call wannier_writeinfo_overall
        do wf_group = 1, wf_ngroups
          call wannier_writeinfo_task
          call wannier_maxloc
        end do
        call wannier_writetransform
      end if

      if( input%properties%wannier%do .ne. "fromfile") call wannier_writeinfo_finish

      ! further tasks if requested
      ! bandstructure
      if( associated( input%properties%bandstructure)) then
        if( input%properties%bandstructure%wannier) call wfutil_bandstructure
      end if
      ! density of states
      if( associated( input%properties%dos)) then
        if( input%properties%dos%wannier) call wfutil_dos
      end if
      ! band gap
      if( associated( input%properties%wanniergap)) then
        call wfutil_find_bandgap
      end if
    else
      write(*,*) " Error (wannierlauncher): Wannier element in input not found."
      stop
    end if
end subroutine wannierlauncher
