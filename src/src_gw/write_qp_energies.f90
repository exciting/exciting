
subroutine write_qp_energies(fname)

  use modinput
  use modgw
  use mod_selfenergy, only : evalks, evalqp, selfex, sigc, znorm, sigsx, sigch
  use mod_vxc,        only : vxcnn
  use m_getunit
  implicit none
  character(len=*), intent(in) :: fname

  integer(4) :: ie
  integer(4) :: ikp
  integer(4) :: fid
  real(8)    :: de, dx
  real(8)    :: ehf, eks, egw
  real(8)    :: vxc, sx, scr, sci, z

  !---------------------------------------------------------------------------
  call getunit(fid)
  open(fid,file=trim(fname),action='WRITE',form='FORMATTED')

  do ikp = 1, kset%nkpt

    write(fid,1) ikp, kset%vkl(:,ikp), kset%wkpt(ikp)
    write(fid,2)

    do ie = ibgw, nbgw

      eks = evalks(ie,ikp)
      egw = evalqp(ie,ikp)
      de = egw-eks

      vxc = dble(vxcnn(ie,ikp))

      select case(input%gw%taskname)

        case('g0w0','gw0','evalqp')
          sx  = dble(selfex(ie,ikp))
          scr = dble(sigc(ie,ikp))
          sci = aimag(sigc(ie,ikp))
          z   = znorm(ie,ikp)

        case('g0w0-x')
          sx = dble(selfex(ie,ikp))
          scr = 0.d0
          sci = 0.d0
          z = 0.d0

        case('cohsex')
          sx  = dble(sigsx(ie,ikp))
          scr = dble(sigch(ie,ikp))
          sci = aimag(sigch(ie,ikp))
          z   = 0.d0

      end select

      dx = sx-vxc
      ehf = eks+dx

      write(fid,3) ie, eks, ehf, egw, &
      &            sx, scr, sci, vxc, dx, de, z

    enddo ! ie
    write(fid,*)

  enddo ! ikp
  close(fid)

  1 format('k-point #',i6,':',4f12.6)
  2 format(' state   E_KS       E_HF       E_GW       Sx         Re(Sc)     Im(Sc)     Vxc        DE_HF      DE_GW      Znk')
  3 format(i4,'  ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)

  return
end subroutine
