
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modxcifc

use libxcifc

contains

!BOP
! !ROUTINE: xcifc
! !INTERFACE:
subroutine xcifc(xctype,n,rho,rhoup,rhodn,grho,gup,gdn,g2rho,g2up,g2dn,g3rho, &
 g3up,g3dn,grho2,gup2,gdn2,gupdn,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn,dxdg2,dxdgu2, &
 dxdgd2,dxdgud,dcdg2,dcdgu2,dcdgd2,dcdgud)
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer(3))
!   n      : number of density points (in,integer)
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
!   grho2  : |grad rho|^2 (in,real(n),optional)
!   gup2   : |grad rhoup|^2 (in,real(n),optional)
!   gdn2   : |grad rhodn|^2 (in,real(n),optional)
!   gupdn  : (grad rhoup).(grad rhodn) (in,real(n),optional)
!   ex     : exchange energy density (out,real(n),optional)
!   ec     : correlation energy density (out,real(n),optional)
!   vx     : spin-unpolarised exchange potential (out,real(n),optional)
!   vc     : spin-unpolarised correlation potential (out,real(n),optional)
!   vxup   : spin-up exchange potential (out,real(n),optional)
!   vxdn   : spin-down exchange potential (out,real(n),optional)
!   vcup   : spin-up correlation potential (out,real(n),optional)
!   vcdn   : spin-down correlation potential (out,real(n),optional)
!   dxdg2  : de_x/d(|grad rho|^2) (out,real(n),optional)
!   dxdgu2 : de_x/d(|grad rhoup|^2) (out,real(n),optional)
!   dxdgd2 : de_x/d(|grad rhodn|^2) (out,real(n),optional)
!   dxdgud : de_x/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
!   dcdg2  : de_c/d(|grad rho|^2) (out,real(n),optional)
!   dcdgu2 : de_c/d(|grad rhoup|^2) (out,real(n),optional)
!   dcdgd2 : de_c/d(|grad rhodn|^2) (out,real(n),optional)
!   dcdgud : de_c/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
! !DESCRIPTION:
!   Interface to the exchange-correlation routines.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Modiefied Januar 2013 (UW)
!EOP
!BOC
implicit none
! mandatory arguments
integer, intent(in) :: xctype(3)
integer, intent(in) :: n
! optional arguments
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
real(8), optional, intent(in) :: grho2(*)
real(8), optional, intent(in) :: gup2(*)
real(8), optional, intent(in) :: gdn2(*)
real(8), optional, intent(in) :: gupdn(*)
real(8), optional, intent(out) :: ex(*)
real(8), optional, intent(out) :: ec(*)
real(8), optional, intent(out) :: vx(*)
real(8), optional, intent(out) :: vc(*)
real(8), optional, intent(out) :: vxup(*)
real(8), optional, intent(out) :: vxdn(*)
real(8), optional, intent(out) :: vcup(*)
real(8), optional, intent(out) :: vcdn(*)
real(8), optional, intent(out) :: dxdg2(*)
real(8), optional, intent(out) :: dxdgu2(*)
real(8), optional, intent(out) :: dxdgd2(*)
real(8), optional, intent(out) :: dxdgud(*)
real(8), optional, intent(out) :: dcdg2(*)
real(8), optional, intent(out) :: dcdgu2(*)
real(8), optional, intent(out) :: dcdgd2(*)
real(8), optional, intent(out) :: dcdgud(*)
! local variables
real(8) kappa,mu,beta
! automatic arrays
real(8), allocatable :: ra(:,:)
if (n.le.0) then
  write(*,*)
  write(*,'("Error(xcifc): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
#ifdef DEBUGLIBXC
write(*,*)"(xcifc)xctype:", xctype
#endif
select case(abs(xctype(1)))
case(1)
! No density-derived exchange-correlation energy or potential
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
  if (present(rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    call xc_pzca(n,rho,ex,ec,vx,vc)
  else
    goto 10
  end if
case(3,407,408)
! Perdew-Wang parameterisation of the spin-polarised Ceperley-Alder electron gas
! J. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
  if (present(rhoup).and.present(rhodn).and.present(ex).and.present(ec) &
   .and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn)) then
! spin-polarised density
    call xc_pwca(n,rhoup,rhodn,ex,ec,vxup,vxdn,vcup,vcdn)
  else if (present(rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
! divide spin-unpolarised density into up and down
    allocate(ra(n,1))
    ra(1:n,1)=0.5d0*rho(1:n)
    call xc_pwca(n,ra(:,1),ra(:,1),ex,ec,vx,vx,vc,vc)
    deallocate(ra)
  else
    goto 10
  end if
case(4)
! X-alpha approximation
! J. C. Slater, Phys. Rev. 81, 385 (1951)
  if (present(rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    call xc_xalpha(n,rho,ex,vx)
! set correlation energy and potential to zero
    ec(1:n)=0.d0
    vc(1:n)=0.d0
  else
    goto 10
  end if
case(5)
! U. von Barth and L. Hedin parameterisation of LSDA
! J. Phys. C, 5, 1629 (1972)
  if (present(rhoup).and.present(rhodn).and.present(ex).and.present(ec) &
   .and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn)) then
! spin-polarised density
    call xc_vbh(n,rhoup,rhodn,ex,ec,vxup,vxdn,vcup,vcdn)
  else if (present(rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
! divide spin-unpolarised density into up and down
    allocate(ra(n,1))
    ra(1:n,1)=0.5d0*rho(1:n)
    call xc_vbh(n,ra(:,1),ra(:,1),ex,ec,vx,vx,vc,vc)
    deallocate(ra)
  else
    goto 10
  end if
case(20,21,22,300,406)
! original PBE kappa
  kappa=0.804d0
  if (xctype(1).eq.21) then
! Zhang-Yang kappa
    kappa=1.245d0
  end if
! original PBE mu and beta
  mu=0.2195149727645171d0
  beta=0.06672455060314922d0
  if (xctype(1).eq.22) then
! PBEsol parameters
    mu=10.d0/81.d0
    beta=0.046d0
  end if
  if (xctype(1).eq.300) then
! beta for acPBE
    mu=0.249d0 
  end if
! Perdew-Burke-Ernzerhof generalised gradient approximation
! Phys. Rev. Lett. 77, 3865 (1996); 78, 1396(E) (1997)
! Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998)
  if (present(rhoup).and.present(rhodn).and.present(grho).and.present(gup) &
   .and.present(gdn).and.present(g2up).and.present(g2dn).and.present(g3rho) &
   .and.present(g3up).and.present(g3dn).and.present(ex).and.present(ec) &
   .and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn)) then
    call xc_pbe(n,kappa,mu,beta,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up, &
     g3dn,ex,ec,vxup,vxdn,vcup,vcdn)
  else if (present(rho).and.present(grho).and.present(g2rho) &
   .and.present(g3rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    allocate(ra(n,6))
    ra(1:n,1)=0.5d0*rho(1:n)
    ra(1:n,2)=0.5d0*grho(1:n)
    ra(1:n,3)=0.5d0*g2rho(1:n)
    ra(1:n,4)=0.25d0*g3rho(1:n)
    call xc_pbe(n,kappa,mu,beta,ra(:,1),ra(:,1),grho,ra(:,2),ra(:,2),ra(:,3), &
     ra(:,3),g3rho,ra(:,4),ra(:,4),ex,ec,vx,ra(:,5),vc,ra(:,6))
    deallocate(ra)
  else
    goto 10
  end if
case(26)
! Wu-Cohen exchange with PBE correlation generalised gradient functional
! Zhigang Wu and R. E. Cohen, Phys. Rev. B 73, 235116 (2006)
  if (present(rho).and.present(grho).and.present(g2rho).and.present(g3rho) &
   .and.present(ex).and.present(ec).and.present(vx).and.present(vc)) then
    call xc_wc06(n,rho,grho,g2rho,g3rho,ex,ec,vx,vc)
  else
    goto 10
  end if
case(30)
! Armiento-Mattsson generalised gradient functional
! R. Armiento and A. E. Mattsson, Phys. Rev. B 72, 085108 (2005)
  if (present(rho).and.present(grho).and.present(g2rho).and.present(g3rho) &
   .and.present(ex).and.present(ec).and.present(vx).and.present(vc)) then
    call xc_am05(n,rho,grho,g2rho,g3rho,ex,ec,vx,vc)
  else
    goto 10
  end if
case(100)
! libxc library functionals
  if (present(rhoup).and.present(rhodn).and.present(gup2).and.present(gdn2) &
   .and.present(gupdn).and.present(ex).and.present(ec).and.present(vxup) &
   .and.present(vxdn).and.present(vcup).and.present(vcdn).and.present(dxdgu2) &
   .and.present(dxdgd2).and.present(dxdgud).and.present(dcdgu2) &
   .and.present(dcdgd2).and.present(dcdgd2).and.present(dcdgud)) then
    call xcifc_libxc(xctype,n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
     gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
     dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
     dcdgud=dcdgud)
  else if (present(rhoup).and.present(rhodn).and.present(ex).and.present(ec) &
   .and.present(vxup).and.present(vxdn).and.present(vcup) &
   .and.present(vcdn)) then
    call xcifc_libxc(xctype,n,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec,vxup=vxup, &
     vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (present(rho).and.present(grho2).and.present(ex).and.present(ec) &
   .and.present(vx).and.present(vc).and.present(dxdg2).and.present(dcdg2)) then
    call xcifc_libxc(xctype,n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
     dxdg2=dxdg2,dcdg2=dcdg2)
  else if (present(rho).and.present(ex).and.present(ec).and.present(vx) &
   .and.present(vc)) then
    call xcifc_libxc(xctype,n,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
  else
    goto 10
  end if
case default
  write(*,*)
  write(*,'("Error(xcifc): xctype not defined : ",I8)') xctype(1)
  write(*,*)
  stop
end select
! set exchange potential to zero for EXX
if (xctype(1).le.-2) then
  if (present(vx)) vx(1:n)=0.d0
  if (present(vxup)) vxup(1:n)=0.d0
  if (present(vxdn)) vxdn(1:n)=0.d0
end if
return
10 continue
write(*,*)
write(*,'("Error(xcifc): missing arguments for exchange-correlation type ",&
 &3I6)') xctype(:)
write(*,*)
stop
end subroutine
!EOC

!BOP
! !ROUTINE: getxcdata
! !INTERFACE:
subroutine getxcdata(xctype,xcdescr,xcspin,xcgrad,ex_coef)
! !INPUT/OUTPUT PARAMETERS:
!   xctype  : type of exchange-correlation functional (in,integer(3))
!   xcdescr : description of functional (out,character(256))
!   xcspin  : spin treatment (out,integer)
!   xcgrad  : gradient treatment (out,integer)
! !DESCRIPTION:
!   Returns data on the exchange-correlation functional labeled by {\tt xctype}.
!   The character array {\tt xcdescr} contains a short description of the
!   functional including journal references. The variable {\tt xcspin} is set to
!   1 or 0 for spin-polarised or -unpolarised functionals, respectively. For
!   functionals which require the gradients of the density {\tt xcgrad} is set
!   to 1, otherwise it is set to 0.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: xctype(3)
character(512), intent(out) :: xcdescr
integer, intent(out) :: xcspin
integer, intent(out) :: xcgrad
real(8), intent(out) :: ex_coef
select case(abs(xctype(1)))
case(1)
  xcdescr='No density-derived exchange-correlation energy or potential'
! spin-polarisation or gradient status not required
  xcspin=-1
  xcgrad=-1
  return
case(2)
  xcdescr='Perdew-Zunger, Phys. Rev. B 23, 5048 (1981)'
  xcspin=0
  xcgrad=0
  return
case(3)
  xcdescr='Perdew-Wang, Phys. Rev. B 45, 13244 (1992)'
  xcspin=1
  xcgrad=0
case(4)
  xcdescr='X-alpha approximation, J. C. Slater, Phys. Rev. 81, 385 (1951)'
  xcspin=0
  xcgrad=0
case(5)
  xcdescr='von Barth-Hedin, J. Phys. C 5, 1629 (1972)'
  xcspin=1
  xcgrad=0
case(20)
  xcdescr='Perdew-Burke-Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)'
  xcspin=1
  xcgrad=1
case(21)
  xcdescr='Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998)'
  xcspin=1
  xcgrad=1
case(22)
  xcdescr='PBEsol, Phys. Rev. Lett. 100, 136406 (2008)'
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
case(300)
  xcdescr='acPBE (asymptotically corrected PBE), arXiv:1409.4834 (2014)'
  xcspin=1
  xcgrad=1
case(406)
  xcdescr='PBE0, M. Ernzerhof, G. E. Scuseria, J. Chem. Phys. 110 , 5029 (1999)'
  xcspin=0
  xcgrad=1
case(407)
  xcdescr='LDA0 (test only)'
  xcspin=1
  xcgrad=0
case(100)
! libxc library functionals
  call xcdata_libxc(xctype,xcdescr,xcspin,xcgrad,ex_coef)
  !input%groundstate%Hybrid%excoeff=ex_coef
case default
  write(*,*)
  write(*,'("Error(getxcdata): xctype not defined : ",I8)') xctype(1)
  write(*,*)
  stop
end select
return
end subroutine
!EOC

end module

