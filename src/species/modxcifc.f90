
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE:  modxcifc
! !DESCRIPTION:
!   Exchange-correlation functional interface block.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
module modxcifc
interface
subroutine xcifc(xctype,n,rho,rhoup,rhodn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn)
implicit none
! mandatory arguments
integer, intent(in) :: xctype
! optional arguments
integer, optional :: n
real(8), optional :: rho(*)
real(8), optional :: rhoup(*)
real(8), optional :: rhodn(*)
real(8), optional :: ex(*)
real(8), optional :: ec(*)
real(8), optional :: vx(*)
real(8), optional :: vc(*)
real(8), optional :: vxup(*)
real(8), optional :: vxdn(*)
real(8), optional :: vcup(*)
real(8), optional :: vcdn(*)
end subroutine
end interface
end module

!BOP
! !ROUTINE: xcifc
! !INTERFACE:
subroutine xcifc(xctype,n,rho,rhoup,rhodn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer)
!   n      : number of density points (in,integer,optional)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
!   ex     : exchange energy density (out,real(n),optional)
!   ec     : correlation energy density (out,real(n),optional)
!   vx     : spin-unpolarised exchange potential (out,real(n),optional)
!   vc     : spin-unpolarised correlation potential (out,real(n),optional)
!   vxup   : spin-up exchange potential (out,real(n),optional)
!   vxdn   : spin-down exchange potential (out,real(n),optional)
!   vcup   : spin-up correlation potential (out,real(n),optional)
!   vcdn   : spin-down correlation potential (out,real(n),optional)
! !DESCRIPTION:
!   Interface to the exchange-correlation routines. This makes it relatively
!   simple to add new functionals which do not necessarily depend only on $\rho$.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! mandatory arguments
integer, intent(in) :: xctype
! optional arguments
integer, optional :: n
real(8), optional :: rho(*)
real(8), optional :: rhoup(*)
real(8), optional :: rhodn(*)
real(8), optional :: ex(*)
real(8), optional :: ec(*)
real(8), optional :: vx(*)
real(8), optional :: vc(*)
real(8), optional :: vxup(*)
real(8), optional :: vxdn(*)
real(8), optional :: vcup(*)
real(8), optional :: vcdn(*)
! local variables
select case(abs(xctype))
case(1)
! No density-derived exchange-correlation energy or potential
  if (.not.(present(n))) goto 10
  if (n.le.0) goto 20
  if (present(ex)) ex(1:n)=0.d0
  if (present(ec)) ec(1:n)=0.d0
  if (present(vx)) vx(1:n)=0.d0
  if (present(vc)) vc(1:n)=0.d0
  if (present(vxup)) vxup(1:n)=0.d0
  if (present(vxdn)) vxdn(1:n)=0.d0
  if (present(vcup)) vcup(1:n)=0.d0
  if (present(vcdn)) vcdn(1:n)=0.d0
  return
case(2)
! Perdew-Zunger parameterisation of Ceperley-Alder electron gas
! J. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
  if (.not.(present(n).and.present(rho).and.present(ex).and.present(ec) &
   .and.present(vx).and.present(vc))) goto 10
  if (n.le.0) goto 20
  call xc_pzca(n,rho,ex,ec,vx,vc)
  return
case(3)
! Perdew-Wang parameterisation of the spin-polarised Ceperley-Alder electron gas
! J. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
  if (.not.(present(n).and.present(rhoup).and.present(rhodn).and.present(ex) &
   .and.present(ec).and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn))) goto 10
  if (n.le.0) goto 20
  call xc_pwca(n,rhoup,rhodn,ex,ec,vxup,vxdn,vcup,vcdn)
  return
case default
  write(*,*)
  write(*,'("Error(xcifc): xctype not defined : ",I8)') xctype
  write(*,*)
  stop
end select
10 continue
write(*,*)
write(*, &
 '("Error(xcifc): missing arguments for exchange-correlation type ",I5)') xctype
write(*,*)
stop
20 continue
write(*,*)
write(*,'("Error(xcifc): n <= 0 : ",I8)') n
write(*,*)
stop
end subroutine
!EOC

!BOP
! !ROUTINE: getxcdata
! !INTERFACE:
subroutine getxcdata(xctype,xcdescr,xcspin,xcgrad)
! !INPUT/OUTPUT PARAMETERS:
!   xctype  : type of exchange-correlation functional (in,integer)
!   xcdescr : description of functional (out,character(256))
!   xcspin  : spin treatment (out,integer)
!   xcgrad  : gradient treatment (out,integer)
! !DESCRIPTION:
!   Returns data on the exchange-correlation functional labelled by {\tt xctype}.
!   The character array {\tt xctype} contains a short description of the
!   functional including journal references. The variable {\tt xcspin} is set to
!   1 or 0 for spin-polarised or -unpolarised functionals, respectively. For
!   functionals which require the gradient of the density {\tt xcgrad} is set to
!   1, otherwise it is set to 0.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: xctype
character(256), intent(out) :: xcdescr
integer, intent(out) :: xcspin
integer, intent(out) :: xcgrad
! local variables
select case(abs(xctype))
case(1)
  xcdescr='No density-derived exchange-correlation energy or potential'
  xcspin=-1
  xcgrad=-1
  return
case(2)
  xcdescr='Perdew-Zunger/Ceperley-Alder, Phys. Rev. B 23, 5048 (1981)'
  xcspin=0
  xcgrad=0
  return
case(3)
  xcdescr='Perdew-Wang/Ceperley-Alder, Phys. Rev. B 45, 13244 (1992)'
  xcspin=1
  xcgrad=0
case default
  write(*,*)
  write(*,'("Error(getxcdata): xctype not defined : ",I8)') xctype
  write(*,*)
  stop
end select
return
end subroutine
!EOC
