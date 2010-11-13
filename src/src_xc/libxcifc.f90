! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
#define LIBXC
module libxcifc
#ifdef LIBXC
use xc_f90_lib_m
#endif
contains

!BOP
! !ROUTINE: xcifc_libxc
! !INTERFACE:
subroutine xcifc_libxc(xctype,n,rho,rhoup,rhodn,grho2,gup2,gdn2,gupdn,ex,ec, &
 vx,vc,vxup,vxdn,vcup,vcdn,dxdg2,dxdgu2,dxdgd2,dxdgud,dcdg2,dcdgu2,dcdgd2, &
 dcdgud)



! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer(3))
!   n      : number of density points (in,integer)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
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
!   Interface to the {\tt libxc} exchange-correlation functional library:
!   \newline{\tt http://www.tddft.org/programs/octopus/wiki/index.php/Libxc}.
!   The second and third integers in {\tt xctype} define the exchange and
!   correlation functionals in {\tt libxc}, respectively.
!
! !REVISION HISTORY:
!   Created April 2009 (Tyrel McQueen)
!   Modified September 2009 (JKD and TMQ)
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
#ifdef LIBXC
! local variables
integer nspin,xcf,id,i,k
real(8) r(2),v(2),sigma(3),vsigma(3)
type(xc_f90_pointer_t) p
type(xc_f90_pointer_t) info
if (present(rho)) then
  nspin=XC_UNPOLARIZED
else if (present(rhoup).and.present(rhodn)) then
  nspin=XC_POLARIZED
else
  write(*,*)
  write(*,'("Error(xcifc_libxc): missing arguments")')
  write(*,*)
  stop
end if
! loop over functional kinds (exchange or correlation)
do k=2,3
  id=xctype(k)
  if (id.gt.0) then
    xcf=xc_f90_family_from_id(id)
    select case(xcf)
    case(XC_FAMILY_LDA)
!-------------------------!
!     LDA functionals     !
!-------------------------!
      if (id.eq.XC_LDA_X) then
        call xc_f90_func_init(p,info,id,nspin)!,3,XC_NON_RELATIVISTIC)
      else if (id.eq.XC_LDA_C_XALPHA) then
        call xc_f90_func_init(p,info,id,nspin)!3,1.d0)
      else
        call xc_f90_func_init(p,info,id,nspin)
      end if
      if (k.eq.2) then
! exchange
        if (present(rho)) then
          do i=1,n
            call xc_f90_lda_exc_vxc(p,1,rho(i), ex(i),vx(i))
          end do
        else
          do i=1,n
            r(1)=rhoup(i); r(2)=rhodn(i)
            call xc_f90_lda_exc_vxc(p,1,r(1), ex(i),v(1))
            vxup(i)=v(1); vxdn(i)=v(2)
          end do
        end if
      else
! correlation
        if (present(rho)) then
          do i=1,n
            call xc_f90_lda_exc_vxc(p,1,rho(i), ec(i),vc(i))
          end do
        else
          do i=1,n
            r(1)=rhoup(i); r(2)=rhodn(i)
            call xc_f90_lda_exc_vxc(p,1,r(1), ec(i),v(1))
            vcup(i)=v(1); vcdn(i)=v(2)
          end do
        end if
      end if
! destroy functional
      call xc_f90_func_end(p)
    case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
!-------------------------!
!     GGA functionals     !
!-------------------------!
      call xc_f90_func_init(p,info,id,nspin)
      if (k.eq.2) then
! exchange
        if (present(rho)) then
            call xc_f90_gga_exc_vxc(p,n,rho(1),grho2(1),ex(1), vx(1),dxdg2(1))
        else
          do i=1,n
            r(1)=rhoup(i); r(2)=rhodn(i)
            sigma(1)=gup2(i); sigma(2)=gupdn(i); sigma(3)=gdn2(i)
            call xc_f90_gga_exc_vxc(p,1,r(1),sigma(1), ex(i),v(1),vsigma(1))
            vxup(i)=v(1); vxdn(i)=v(2)
            dxdgu2(i)=vsigma(1); dxdgud(i)=vsigma(2); dxdgd2(i)=vsigma(3)
          end do
        end if
      else
! correlation
        if (present(rho)) then
          do i=1,n
            call xc_f90_gga_exc_vxc(p,1,rho(i),grho2(i), ec(i),vc(i),dcdg2(i))
          end do
        else
          do i=1,n
            r(1)=rhoup(i); r(2)=rhodn(i)
            sigma(1)=gup2(i); sigma(2)=gupdn(i); sigma(3)=gdn2(i)
            call xc_f90_gga_exc_vxc(p,1,r(1),sigma(1), ec(i),v(1),vsigma(1))
            vcup(i)=v(1); vcdn(i)=v(2)
            dcdgu2(i)=vsigma(1); dcdgud(i)=vsigma(2); dcdgd2(i)=vsigma(3)
          end do
        end if
      end if
    case default
      write(*,*)
      write(*,'("Error(xcifc_libxc): unsupported libxc functional family : ",&
       &I8)') xcf
      write(*,*)
      stop
    end select
  else
! case when id=0
    if (k.eq.2) then
      ex(1:n)=0.d0
      if (present(rho)) then
        vx(1:n)=0.d0
      else
        vxup(1:n)=0.d0
        vxdn(1:n)=0.d0
      end if
    else
      ec(1:n)=0.d0
      if (present(rho)) then
        vc(1:n)=0.d0
      else
        vcup(1:n)=0.d0
        vcdn(1:n)=0.d0
      end if
    end if
  end if
end do
return

#endif

end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad)

implicit none
! arguments
integer, intent(in) :: xctype(3)
character(512), intent(out) :: xcdescr
integer, intent(out) :: xcspin
integer, intent(out) :: xcgrad
#ifdef LIBXC

! local variables
integer xcf,id,k
character(256) name
type(xc_f90_pointer_t) p
type(xc_f90_pointer_t) info
! unknown spin polarisation
xcspin=-1
! no gradients by default
xcgrad=0
do k=2,3
  id=xctype(k)
  if (id.gt.0) then
    xcf=xc_f90_family_from_id(id)
    select case(xcf)
    case(XC_FAMILY_LDA)
      if (id.eq.XC_LDA_X) then
        call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)!,3,XC_NON_RELATIVISTIC)
      else if (id.eq.XC_LDA_C_XALPHA) then
        call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)!,3,1.d0)
      else
        call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)
      end if
      call xc_f90_info_name(info,name)
      call xc_f90_func_end(p)
    case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
      call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)
      call xc_f90_info_name(info,name)
      call xc_f90_func_end(p)
! post-processed gradients required
      xcgrad=2
    case default
      write(*,*)
      write(*,'("Error(xcdata_libxc): unsupported libxc functional family : ",&
       &I8)') xcf
      write(*,*)
      stop
    end select
  else
    name='none'
  end if
  if (k.eq.2) then
    xcdescr='libxc; exchange: '//trim(name)
  else
    xcdescr=trim(xcdescr)//'; correlation: '//trim(name)
  end if
end do
xcdescr=trim(xcdescr)//' (see libxc for references)'

return
#endif

#ifndef LIBXC
 write(*,'("Error(xcdata_libxc): LIBXC not acivated : ",&
       &I8)')
      write(*,*)
      stop
#endif
end subroutine
!EOC


end module
