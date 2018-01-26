subroutine wannierlauncher
    use modinput
    use mod_wannier
    use mod_wfutil
    implicit none

    if( associated( input%properties%wannier)) then
      ! generate Wannier functions
      call wannier_init
      if( input%properties%wannier%method .eq. "pro") then
        call wannier_gen_pro
      else if( input%properties%wannier%method .eq. "prowan") then
        call wannier_projonwan
        call wannier_gen_pro
      else if( input%properties%wannier%method .eq. "opf") then
        call wannier_gen_opf
      else if( input%properties%wannier%method .eq. "opfwan") then
        call wannier_projonwan
        call wannier_gen_opf
      else if( input%properties%wannier%method .eq. "promax") then
        call wannier_gen_pro
        call wannier_maxloc
      else if( input%properties%wannier%method .eq. "prowanmax") then
        call wannier_projonwan
        call wannier_gen_pro
        call wannier_maxloc
      else if( input%properties%wannier%method .eq. "opfmax") then
        call wannier_gen_opf
        call wannier_maxloc
      else if( input%properties%wannier%method .eq. "maxfromfile") then
        call wannier_gen_fromfile
        call wannier_maxloc
      else if( input%properties%wannier%method .eq. "fromfile") then
        call wannier_gen_fromfile
      else if( input%properties%wannier%method .eq. "disentangle") then

      else
        write(*,*) " Error (propertylauncher): invalid value for attribute method"
      end if
      call wannier_writeinfo_finish

      ! further tasks if requested
      ! bandstructure
      if( associated( input%properties%bandstructure)) then
        if( input%properties%bandstructure%wannier) call wfutil_bandstructure
      end if
      ! density of states
      if( associated( input%properties%dos)) then
        if( input%properties%dos%wannier) call wfutil_dos
      end if
    else
      write(*,*) " Error (wannierlauncher): Wannier element in input not found."
      call terminate
    end if
end subroutine wannierlauncher
