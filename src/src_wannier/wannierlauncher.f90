subroutine wannierlauncher
    use modinput
    use mod_wannier
    use mod_wannier_util
    use mod_wannier_maxloc
    use mod_wannier_interpolate
    use mod_manopt

    implicit none

    if( associated( input%properties%wannier)) then
      ! generate Wannier functions
      call wannier_init
      if( input%properties%wannier%do .eq. "fromscratch") then
        call wffile_writeinfo_overall
        if( mpiglobal%rank .eq. 0) call printbox( wf_info, '*', "Wannierization")
        do wf_group = 1, wf_ngroups
          call wffile_writeinfo_task
          if( wf_groups( wf_group)%method .eq. "pro") then
            if( mpiglobal%rank .eq. 0) call wannier_gen_pro
          else if( wf_groups( wf_group)%method .eq. "opf") then
            if( mpiglobal%rank .eq. 0) call wfopf_gen
          else if( wf_groups( wf_group)%method .eq. "promax") then
            if( mpiglobal%rank .eq. 0) call wannier_gen_pro
            if( mpiglobal%rank .eq. 0) call wfmax_gen
          else if( wf_groups( wf_group)%method .eq. "opfmax") then
            if( mpiglobal%rank .eq. 0) call wfopf_gen
            if( mpiglobal%rank .eq. 0) call wfmax_gen
          else if( wf_groups( wf_group)%method .eq. "disSMV" .or. wf_groups( wf_group)%method .eq. "disFull") then
            if( mpiglobal%rank .eq. 0) call wfopf_gen
            if( mpiglobal%rank .eq. 0) call wfdis_gen
            if( mpiglobal%rank .eq. 0) call wfopf_gen( subspace=.true.)
            if( mpiglobal%rank .eq. 0) call wfmax_gen
          else
            write(*,*) " Error (propertylauncher): invalid value for attribute method"
          end if
        end do
        call wffile_writetransform
      
      else if( input%properties%wannier%do .eq. "fromfile") then
        call wannier_gen_fromfile

      else
        call wannier_gen_fromfile
        call wffile_writeinfo_overall
        do wf_group = 1, wf_ngroups
          call wffile_writeinfo_task
          if( mpiglobal%rank .eq. 0) call wfmax_fromfile
        end do
        call wffile_writetransform
      end if

      if( input%properties%wannier%do .ne. "fromfile") call wffile_writeinfo_finish

      ! share results over all processes
#ifdef MPI
      call barrier
      if( input%properties%wannier%do .ne. "fromfile") then
        call mpi_bcast( wf_transform, wf_nst*wf_nwf*wf_kset%nkpt, mpi_double_complex, 0, mpiglobal%comm, ierr)
        do wf_group = 1, wf_ngroups
          call wfomega_gen
        end do
      end if
#endif

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
      ! visualization
      if( associated( input%properties%wannierplot)) then
        if( mpiglobal%rank .eq. 0) call wfutil_plot( input%properties%wannierplot%fst, input%properties%wannierplot%lst, input%properties%wannierplot%cell)
      end if
    else
      write(*,*) " Error (wannierlauncher): Wannier element in input not found."
      stop
    end if

end subroutine wannierlauncher
