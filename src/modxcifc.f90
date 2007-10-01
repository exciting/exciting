
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modxcifc
contains

!BOP
! !ROUTINE: xcifc
! !INTERFACE:
subroutine xcifc(xctype,n,rho,rhoup,rhodn,grho,gup,gdn,g2rho,g2up,g2dn,g3rho, &
 g3up,g3dn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer)
!   n      : number of density points (in,integer,optional)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
!   grho   : |grad rho| (in,real(n),optional)
!   gup    : |grad rhoup| (in,real(n),optional)
!   gdn    : |grad rhodn| (in,real(n),optional)
!   g2rho  : grad^2 rho (in,real(n),optional)
!   g2up   : grad^2 rhoup (in,real(n),optional)
!   g2dn   : grad^2 rhodn (in,real(n),optional)
!   g3rho  : (grad rho).(grad |grad rho|) (in,real(n),optional)
!   g3up   : (grad rhoup).(grad |grad rhoup|) (in,real(n),optional)
!   g3dn   : (grad rhodn).(grad |grad rhodn|) (in,real(n),optional)
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
!   simple to add new functionals which do not necessarily depend only on
!   $\rho$.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! mandatory arguments
integer, intent(in) :: xctype
! optional arguments
integer, optional, intent(in) :: n
real(8), optional, intent(in) :: rho(*)
real(8), optional, intent(in) :: rhoup(*)
real(8), optional, intent(in) :: rhodn(*)
real(8), optional, intent(in) :: grho(*)
real(8), optional, intent(in) :: gup(*)
real(8), optional, intent(in) :: gdn(*)
real(8), optional, intent(in) :: g2rho(*)
real(8), optional, intent(in) :: g2up(*)
real(8), optional, intent(in) :: g2dn(*)
real(8), optional, intent(in) :: g3rho(*)
real(8), optional, intent(in) :: g3up(*)
real(8), optional, intent(in) :: g3dn(*)
real(8), optional, intent(out) :: ex(*)
real(8), optional, intent(out) :: ec(*)
real(8), optional, intent(out) :: vx(*)
real(8), optional, intent(out) :: vc(*)
real(8), optional, intent(out) :: vxup(*)
real(8), optional, intent(out) :: vxdn(*)
real(8), optional, intent(out) :: vcup(*)
real(8), optional, intent(out) :: vcdn(*)
! local variables
real(8) kappa
! automatic arrays
real(8), allocatable :: ra(:,:)
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
case(2)
! Perdew-Zunger parameterisation of Ceperley-Alder electron gas
! J. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
  if (present(n).and.present(rho).and.present(ex).and.present(ec) &
   .and.present(vx).and.present(vc)) then
    if (n.le.0) goto 20
    call xc_pzca(n,rho,ex,ec,vx,vc)
  else
    goto 10
  end if
case(3)
! Perdew-Wang parameterisation of the spin-polarised Ceperley-Alder electron gas
! J. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
  if (present(n).and.present(rhoup).and.present(rhodn).and.present(ex) &
   .and.present(ec).and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn)) then
! spin-polarised density
    if (n.le.0) goto 20
    call xc_pwca(n,rhoup,rhodn,ex,ec,vxup,vxdn,vcup,vcdn)
  else if (present(n).and.present(rho).and.present(ex).and.present(ec) &
   .and.present(vx).and.present(vc)) then
! divide spin-unpolarised density into up and down
    if (n.le.0) goto 20
    allocate(ra(n,1))
    ra(1:n,1)=0.5d0*rho(1:n)
    call xc_pwca(n,ra(1,1),ra(1,1),ex,ec,vx,vx,vc,vc)
    deallocate(ra)
  else
    goto 10
  end if
case(4)
! X-alpha approximation
! J. C. Slater, Phys. Rev. 81, 385 (1951)
  if (present(n).and.present(rho).and.present(ex).and.present(ec) &
   .and.present(vx).and.present(vc)) then
    if (n.le.0) goto 20
    call xc_xalpha(n,rho,ex,vx)
! set correlation energy and potential to zero
    ec(1:n)=0.d0
    vc(1:n)=0.d0
  else
    goto 10
  end if
case(20,21)
  if (xctype.eq.20) then
! original
    kappa=0.804d0
  else
! Zhang-Yang
    kappa=1.245d0
  end if
! Perdew-Burke-Ernzerhof generalised gradient approximation
! Phys. Rev. Lett. 77, 3865 (1996); 78, 1396(E) (1997)
! Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998)
  if (present(n).and.present(rhoup).and.present(rhodn).and.present(grho) &
   .and.present(gup).and.present(gdn).and.present(g2up).and.present(g2dn) &
   .and.present(g3rho).and.present(g3up).and.present(g3dn).and.present(ex) &
   .and.present(ec).and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn)) then
    call xc_pbe(n,kappa,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn,ex, &
     ec,vxup,vxdn,vcup,vcdn)
  else if (present(n).and.present(rho).and.present(grho).and.present(g2rho) &
   .and.present(g3rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    allocate(ra(n,4))
    ra(1:n,1)=0.5d0*rho(1:n)
    ra(1:n,2)=0.5d0*grho(1:n)
    ra(1:n,3)=0.5d0*g2rho(1:n)
    ra(1:n,4)=0.25d0*g3rho(1:n)
    call xc_pbe(n,kappa,ra(1,1),ra(1,1),grho,ra(1,2),ra(1,2),ra(1,3),ra(1,3), &
     g3rho,ra(1,4),ra(1,4),ex,ec,vx,vx,vc,vc)
    deallocate(ra)
  else
    goto 10
  end if
case(26)
! Wu-Cohen exchange with PBE correlation generalised gradient functional
! Zhigang Wu and R. E. Cohen, Phys. Rev. B 73, 235116 (2006)
  if (present(n).and.present(rho).and.present(grho).and.present(g2rho) &
   .and.present(g3rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    call xc_wc06(n,rho,grho,g2rho,g3rho,ex,ec,vx,vc)
  else
    goto 10
  end if
case(30)
! Armiento-Mattsson generalised gradient functional
! R. Armiento and A. E. Mattsson, Phys. Rev. B 72, 085108 (2005)
  if (present(n).and.present(rho).and.present(grho).and.present(g2rho) &
   .and.present(g3rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    call xc_am05(n,rho,grho,g2rho,g3rho,ex,ec,vx,vc)
  else
    goto 10
  end if
case default
  write(*,*)
  write(*,'("Error(xcifc): xctype not defined : ",I8)') xctype
  write(*,*)
  stop
end select
! set exchange potential to zero for EXX
if (xctype.le.-2) then
  if (present(vx)) vx(1:n)=0.d0
  if (present(vxup)) vxup(1:n)=0.d0
  if (present(vxdn)) vxdn(1:n)=0.d0
end if
return
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
!   Returns data on the exchange-correlation functional labelled by
!   {\tt xctype}. The character array {\tt xctype} contains a short description
!   of the functional including journal references. The variable {\tt xcspin} is
!   set to 1 or 0 for spin-polarised or -unpolarised functionals, respectively.
!   For functionals which require the gradients of the density {\tt xcgrad} is
!   set to 1, otherwise it is set to 0.
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
select case(abs(xctype))
case(1)
  xcdescr='No density-derived exchange-correlation energy or potential'
! spin-polarisation or gradient status not required
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
case(4)
  xcdescr='X-alpha approximation, J. C. Slater, Phys. Rev. 81, 385 (1951)'
  xcspin=0
  xcgrad=0
case(20)
  xcdescr='Perdew-Burke-Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)'
  xcspin=1
  xcgrad=1
case(21)
  xcdescr='Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998)'
  xcspin=1
  xcgrad=1
case(26)
  xcdescr='Wu-Cohen exchange + PBE correlation, Phys. Rev. B 73, 235116 (2006)'
  xcspin=0
  xcgrad=1
case(30)
  xcdescr='Armiento-Mattsson functional, Phys. Rev. B 72, 85108 (2005)'
  xcspin=0
  xcgrad=1
case default
  write(*,*)
  write(*,'("Error(getxcdata): xctype not defined : ",I8)') xctype
  write(*,*)
  stop
end select
return
end subroutine
!EOC

end module

