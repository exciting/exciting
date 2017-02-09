
subroutine vxcistl(ngp,igpig,v,h)
    use modgw,   only : Gset
    use mod_vxc, only : vxcig
    implicit none
    ! arguments
    integer, intent(in) :: ngp
    integer, intent(in) :: igpig(*)
    complex(8), Intent(in) :: v(*)
    complex(8), intent(inout) :: h(*)
    ! local variables
    integer :: i, j, ig, iv(3)
    complex(8) :: zt1

    do i = 1, ngp
      do j = i, ngp
        iv(:) = Gset%ivg(:,igpig(i))-Gset%ivg(:,igpig(j))
        ig = Gset%ivgig(iv(1),iv(2),iv(3))
        zt1 = vxcig(ig)
        h(i) = h(i)+zt1*v(j)
        if (i.ne.j) h(j) = h(j)+conjg(zt1)*v(i)
      end do
    end do
    
    return
end subroutine
