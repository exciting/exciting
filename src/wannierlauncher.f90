subroutine wannierlauncher
    use modinput
    use mod_wannier
    implicit none

    call wannier_init
    if( input%properties%wannier%method .eq. "pro") then
      call wannier_gen_pro
    else if( input%properties%wannier%method .eq. "opf") then
      call wannier_gen_opf
    else if( input%properties%wannier%method .eq. "promax") then
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
    else
      write(*,*) " Error (propertylauncher): invalid value for attribute method"
    end if
    call wannier_writeinfo_finish
end subroutine wannierlauncher
