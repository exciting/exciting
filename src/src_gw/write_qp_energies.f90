
subroutine write_qp_energies(fname)

    use modinput
    use modgw
    use m_getunit
    implicit none 
    character(len=*), intent(in) :: fname
       
    integer(4) :: ie
    integer(4) :: ikp
    integer(4) :: fid
    real(8)    :: deltae, deltax
    real(8)    :: ehf, eks, egw
    real(8)    :: vxc, sx, sc, z

    !---------------------------------------------------------------------------
    call getunit(fid)
    open(fid,file=trim(fname),action='WRITE',form='FORMATTED')
      
    do ikp = 1, kset%nkpt

      write(fid,1) ikp, kset%vkl(:,ikp), kset%wkpt(ikp)
      write(fid,2)
      do ie = ibgw, nbgw

        eks = evalks(ie,ikp)
        egw = evalqp(ie,ikp)
        deltae = egw-eks
        
        vxc = dble(vxcnn(ie,ikp))
        
        select case(input%gw%taskname)
        
          case('g0w0','gw0','acon')
            sx = dble(selfex(ie,ikp))
            sc = dble(sigc(ie,ikp))
            z = znorm(ie,ikp)
          
          case('g0w0_x')
            sx = dble(selfex(ie,ikp))
            sc = 0.d0
            z = 0.d0
         
          case('cohsex')
            sx = dble(sigsx(ie,ikp))
            sc = dble(sigch(ie,ikp))
            z = 0.d0
        
        end select
          
        deltax = sx-vxc
        ehf = eks+deltax

        write(fid,3) ie, eks, ehf, egw, &
        &            sx, sc, vxc, deltax, deltae, z
        
      enddo ! ie
      write(fid,*)

    enddo ! ikp
    close(fid)

    1 format('k-point #',i6,':',4f12.6)
    2 format(' state    E_KS      E_HF       E_GW       Sx         Sc         Vxc         DE_HF        DE_GW       Znk')    
    3 format(i4,'  ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)

    return 
end subroutine
