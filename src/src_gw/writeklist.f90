
subroutine writeklist(nkpt,vkl,vkc)
    
    implicit none
    
    integer(4), intent(in) :: nkpt          ! Number of k-points
    real(8),    intent(in) :: vkl(3,nkpt)   ! lattice coordinates
    real(8),    intent(in) :: vkc(3,nkpt)   ! cartesian coordinates
    integer(4) :: ik

    write(99,*)
    do ik = 1, nkpt
      write(99,100) ik, vkl(1:3,ik), vkc(1:3,ik)
    enddo
    write(99,*)
          
100 format(i10,3f8.4,4x,3f8.4)

end subroutine writeklist
