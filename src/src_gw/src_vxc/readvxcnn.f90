
subroutine readvxcnn
    use modmain
    use modgw
    use m_getunit
    implicit none
    integer(4) :: ist, i  !(Counter) Runs over bands
    integer(4) :: ik      !(Counter) Runs over k-points
    real(8)    :: kvec(3)
    real(8)    :: vxc_r, vxc_i
    integer    :: fid

    call getunit(fid)
    open(fid,file='VXCNN.OUT',form='FORMATTED',status='OLD')
    do ik = 1, kset%nkpt
      read(fid,*)
      do ist = ibgw, nbgw
        read(fid,*) i, vxc_r, vxc_i
        vxcnn(ist,ik) = cmplx(vxc_r,vxc_i,8)
      end do
      read(fid,*)
    end do
    close(fid)
    
    return
end subroutine
