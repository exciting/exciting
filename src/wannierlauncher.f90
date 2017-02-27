subroutine wannierlauncher
    use modinput
    use mod_wannier, only: wannier_gen_pro, wannier_gen_opf, wannier_gen_max, wannier_gen_fromfile
    implicit none

    if( input%properties%wannier%method .eq. "pro") then
      call wannier_gen_pro( input%properties%wannier%state, input%properties%wannier%nst, input%properties%wannier%projectors)
    else if( input%properties%wannier%method .eq. "opf") then
      call wannier_gen_opf( input%properties%wannier%state, input%properties%wannier%nst)
    else if( input%properties%wannier%method .eq. "promax") then
      call wannier_gen_max( input%properties%wannier%state, input%properties%wannier%nst, input%properties%wannier%projectors)
    else if( input%properties%wannier%method .eq. "opfmax") then
      call wannier_gen_max( input%properties%wannier%state, input%properties%wannier%nst)
    else if( input%properties%wannier%method .eq. "maxfromfile") then
      call wannier_gen_max( input%properties%wannier%state, input%properties%wannier%nst, fromfile=.true.)
    else if( input%properties%wannier%method .eq. "fromfile") then
      call wannier_gen_fromfile
    else
      write(*,*) " Error (propertylauncher): invalid value for attribute method"
    end if
end subroutine wannierlauncher
