
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine genptclg(cuttype,ngpmax,ngp,vgpc,gpc,sptcl)
  use modmain
  use modxs
  implicit none
  ! arguments
  character(*), intent(in) :: cuttype
  integer, intent(in) :: ngpmax,ngp
  real(8), intent(in) :: vgpc(3,ngpmax),gpc(ngpmax)
  real(8), intent(out) :: sptcl(ngpmax)
  ! local variables
  integer :: igp
  ! external functions
  real(8), external :: ptclg
  do igp=1,ngp
     sptcl(igp)=ptclg(cuttype,vgpc(:,igp),gpc(igp))
  end do
end subroutine genptclg


real(8) function ptclg(cuttype,vgpc,gpc)
  use modmain, only: fourpi
  implicit none
  ! arguments
  character(*), intent(in) :: cuttype
  real(8), intent(in) :: vgpc(3),gpc
  ! local variables
  real(8) :: t1
  select case(cuttype)
  case('nocutoff')
     ! set up the square root of the Coulomb potential from analytical
     ! expression (no cutoff)
     ptclg=sqrt(fourpi)/gpc
  case('0d')
     ! 0D spherical cutoff
     t1=vgpc(1)
     write(*,*)
     write(*,'("Error(genptclg): 0D cutoff to be implemented")')
     write(*,*)
     call terminate
  case('1d')
     ! 1D infinite cylinder
     write(*,*)
     write(*,'("Error(genptclg): 1D cutoff to be implemented")')
     write(*,*)
     call terminate
  case('2d')
     ! 2D infinite slab
     write(*,*)
     write(*,'("Error(genptclg): 2D cutoff to be implemented")')
     write(*,*)
     call terminate
  case default
     write(*,*)
     write(*,'("Error(genptclg): unknown type for Coulomb potential: ",a)') &
          cuttype
     write(*,*)
     call terminate
  end select

end function ptclg
