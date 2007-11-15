
module m_genqvkloff
  implicit none
contains

  subroutine genqvkloff(vq,voff)
    use modmain
    use modxs
    implicit none
    ! arguments
    real(8), intent(in) :: vq(3)
    real(8), intent(out) :: voff(3)
    ! local variables
    real(8) :: v2(3)

    if (any(vkloff/dble(ngridk)+vq.ge.1.d0)) then
       ! vector is outside Brillouine zone
       v2=vkloff/dble(ngridk)+vq
       call mapkto01(v2)
       voff=v2*dble(ngridk)
       if (any(v2*dble(ngridk).ge.1.d0)) then
          v2=v2*dble(ngridk)
          call mapkto01(v2)
          voff=v2
       end if
    else if (any(vkloff+vq*dble(ngridk).ge.1.d0)) then
       ! vector is inside Brillouine zone but outside k-point spacing
       v2=vkloff+vq*dble(ngridk)
       call mapkto01(v2)
       voff=v2
    else
       ! vector is inside k-point spacing
       voff=vkloff+vq*ngridk
    end if

!---voff=vkloff+vq*ngridk !SAG

  end subroutine genqvkloff

end module m_genqvkloff
