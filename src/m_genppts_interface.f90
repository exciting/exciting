module m_genppts_interface
  use precision, only: dp, i32
  implicit none

  interface
    subroutine genppts(reducep, tfbz, ngridp, boxl, nppt, ipmap, ivp, vpl, vpc, wppt)
      import :: dp, i32
      logical, intent(in) :: reducep, tfbz
      integer(i32), intent(in) :: ngridp(3)
      real(dp), intent(in) :: boxl(3,4)
      integer(i32), intent(out) :: nppt
      integer(i32), intent(out) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
      integer(i32), allocatable, intent(out) :: ivp(:,:)
      real(dp), intent(out) :: vpl(3,ngridp(1)*ngridp(2)*ngridp(3))
      real(dp), intent(out) :: vpc(3,ngridp(1)*ngridp(2)*ngridp(3))
      real(dp), intent(out) :: wppt(ngridp(1)*ngridp(2)*ngridp(3))
    end subroutine genppts
  end interface
end module m_genppts_interface
