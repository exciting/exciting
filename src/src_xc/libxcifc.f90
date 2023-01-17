! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module libxcifc
#ifdef LIBXC
   use xc_f90_lib_m
#endif
   use precision, only: dp
   use asserts, only: assert

   implicit none
   private
   public :: libxc_lda_potential, &
             libxc_gga_potential, &
             xcdata_libxc

   interface libxc_lda_potential
      procedure libxc_lda_pot_spin_polarized
      procedure libxc_lda_pot_spin_unpolarized
   end interface

   interface libxc_gga_potential
      procedure libxc_gga_pot_spin_polarized
      procedure libxc_gga_pot_spin_unpolarized
   end interface 

contains

   !> Call the LDA libxc functional for the spin-unpolarised case which  
   !> returns the spin-unpolarised LDA potential.
   subroutine libxc_lda_pot_spin_unpolarized(id, n, rho, exc, vxc)
      !> Id of functional
      integer, intent(in) :: id
      !> Number of points
      integer, intent(in) :: n
      !> Spin-unpolarised charge density
      real(dp), intent(in) :: rho(:)
      !> Exchange correlation energy density
      real(dp), intent(out) :: exc(:)
      !> Spin-unpolarised exchange correlation potential
      real(dp), intent(out) :: vxc(:)

#ifdef LIBXC
      type(xc_f90_pointer_t) :: p
      type(xc_f90_pointer_t) :: info
      integer :: family_id 

      if (id == 0) then 
         exc(1:n) = 0.0_dp
         vxc(1:n) = 0.0_dp
         return
      end if 

      family_id = xc_f90_family_from_id(id)
      call assert(family_id == XC_FAMILY_LDA, "Error(libxcifc): Id of functional has & 
                  to belong to the family of LDA functionals.")
      call assert(id > 0, message="Id for xc functional has to be greater than zero.")

      call xc_f90_func_init(p, info, id, XC_UNPOLARIZED)
      call xc_f90_lda_exc_vxc(p, n, rho(1), exc(1), vxc(1))
#endif
   end subroutine

   !> Call the LDA libxc functional for the spin-polarised case which  
   !> returns the spin-polarised LDA potential.
   subroutine libxc_lda_pot_spin_polarized(id, n, rhoup, rhodn, exc, vxcup, vxcdn)
      !> Id of functional
      integer, intent(in) :: id
      !> Number of points
      integer, intent(in) :: n
      !> Spin-polarised charge density spin up
      real(dp), intent(in) :: rhoup(:)
      !> Spin-polarised charge density spin dn
      real(dp), intent(in) :: rhodn(:)
      !> Exchange correlation energy density
      real(dp), intent(out) :: exc(:)
      !> Spin-up exchange correlation potential
      real(dp), intent(out) :: vxcup(:)
      !> Spin-dn exchange correlation potential
      real(dp), intent(out) :: vxcdn(:)

      ! [rhoup, rhodn]
      real(dp) :: rho_total(2)
      ! [vxcup, vxcdn]
      real(dp) :: v_total(2)

#ifdef LIBXC
      type(xc_f90_pointer_t) :: p
      type(xc_f90_pointer_t) :: info
      integer :: i
      integer :: family_id 

      if (id == 0) then 
         exc(1:n) = 0.0_dp
         vxcup(1:n) = 0.0_dp
         vxcdn(1:n) = 0.0_dp
         return
      end if 

      family_id = xc_f90_family_from_id(id)
      call assert(family_id == XC_FAMILY_LDA, "Error(libxcifc): Id of functional has & 
                  to belong to the family of LDA functionals.")
      call assert(id > 0, message="Id for xc functional has to be greater than zero.")

      call xc_f90_func_init(p, info, id, XC_POLARIZED)
      do i = 1, n
         !TODO: pack density (see MR !330)
         rho_total(1) = rhoup(i)
         rho_total(2) = rhodn(i)
         call xc_f90_lda_exc_vxc(p, 1, rho_total(1), exc(i), v_total(1))
         !TODO: unpack potential
         vxcup(i) = v_total(1)
         vxcdn(i) = v_total(2)
      end do
#endif
   end subroutine

   !> Call the GGA libxc functional for the spin-unpolarised case which 
   !> returns the spin-unpolarised GGA potential.
   subroutine libxc_gga_pot_spin_unpolarized(id, n, rho, grho2, exc, vxc, dxcdg2)
      !> Id of functional
      integer, intent(in) :: id
      !> Number of points
      integer, intent(in) :: n
      !> Spin-unpolarised charge density
      real(dp), intent(in) :: rho(:)
      !> Gradient Density squared |grad rho|^2
      real(dp), intent(in) :: grho2(:)
      !> Exchange correlation energy density
      real(dp), intent(out) :: exc(:)
      !> Spin-unpolarised exchange correlation potential
      real(dp), intent(out) :: vxc(:)
      !> de_xc/d(|grad rho|^2)
      real(dp), intent(out) :: dxcdg2(:)

#ifdef LIBXC
      type(xc_f90_pointer_t) :: p
      type(xc_f90_pointer_t) :: info
      integer :: family_id 

      if (id .eq. 0) then  
         exc(1:n) = 0.0_dp
         vxc(1:n) = 0.0_dp
         dxcdg2(1:n) = 0.0_dp 
         return
      end if 

      family_id = xc_f90_family_from_id(id)
      call assert((family_id == XC_FAMILY_GGA) .or. (family_id == XC_FAMILY_HYB_GGA), & 
                  "Error(libxcifc): Id of functional has to belong to the family of GGA functionals.")
      call assert(id > 0, message="Id for xc functional has to be greater than zero.")

      call xc_f90_func_init(p, info, id, XC_UNPOLARIZED)
      call xc_f90_gga_exc_vxc(p, n, rho(1), grho2(1), exc(1), vxc(1), dxcdg2(1))
#endif
   end subroutine

   !> Call the GGA libxc functional for the spin-polarised case which 
   !> returns the spin-polarised GGA potential.
   subroutine libxc_gga_pot_spin_polarized(id, n, rhoup, rhodn, gup2, gdn2, &
                                       gupdn, exc, vxcup, vxcdn, dxcdgu2, dxcdgd2, dxcdgud)
      !> Id of functional
      integer, intent(in) :: id
      !> Number of points
      integer, intent(in) :: n
      !> Spin-polarised charge density spin up
      real(dp), intent(in) :: rhoup(:)
      !> Spin-polarised charge density spin down
      real(dp), intent(in) :: rhodn(:)
      !> Gradient Density spin up squared |grad rhoup|^2
      real(dp), intent(in) :: gup2(:)
      !> Gradient Density spin down squared |grad rhodn|^2
      real(dp), intent(in) :: gdn2(:)
      !> (grad rhoup).(grad rhodn)
      real(dp), intent(in) :: gupdn(:)
      !> Exchange correlation energy density
      real(dp), intent(out) :: exc(:)
      !> Spin-up exchange correlation potential
      real(dp), intent(out) :: vxcup(:)
      !> Spin-down exchange correlation potential
      real(dp), intent(out) :: vxcdn(:)
      !> de_xc/d(|grad rhoup|^2)
      real(dp), intent(out) :: dxcdgu2(:)
      !> de_xc/d(|grad rhodn|^2)
      real(dp), intent(out) :: dxcdgd2(:)
      !> de_xc/d((grad rhoup).(grad rhodn))
      real(dp), intent(out) :: dxcdgud(:)

      ! [rhoup, rhodn]
      real(dp) :: rho_total(2)
      ! [vxcup, vxcdn]
      real(dp) :: v_total(2)
      ! [rhoup*rhoup, rhodn*rhodn, rhoup*rhodn]
      real(dp) :: sigma(3)
      ! [de_xc/d(rhoup*rhoup), de_xc/d(rhodn*rhodn), de_xc/d(rhoup*rhodn)
      real(dp) :: vsigma(3)

#ifdef LIBXC
      type(xc_f90_pointer_t) :: p
      type(xc_f90_pointer_t) :: info
      integer :: i
      integer :: family_id

      if (id .eq. 0) then  
         exc(1:n) = 0.0_dp
         vxcup(1:n) = 0.0_dp
         vxcdn(1:n) = 0.0_dp
         dxcdgu2(1:n) = 0.0_dp 
         dxcdgd2(1:n) = 0.0_dp
         dxcdgud(1:n) = 0.0_dp
         return
      end if 

      family_id = xc_f90_family_from_id(id)
      family_id = 32
      call assert((family_id == XC_FAMILY_GGA) .or. (family_id == XC_FAMILY_HYB_GGA), & 
                  "Error(libxcifc): Id of functional has to belong to the family of GGA functionals.")
      call assert(id > 0, message="Id for xc functional has to be greater than zero.")

      call xc_f90_func_init(p, info, id, XC_POLARIZED)
      do i = 1, n
         !TODO: pack density
         rho_total(1) = rhoup(i)
         rho_total(2) = rhodn(i)
         sigma(1) = gup2(i)
         sigma(2) = gupdn(i)
         sigma(3) = gdn2(i)
         call xc_f90_gga_exc_vxc(p, 1, rho_total(1), sigma(1), exc(i), v_total(1), vsigma(1))
         !TODO: unpack potential
         vxcup(i) = v_total(1)
         vxcdn(i) = v_total(2)
         dxcdgu2(i) = vsigma(1)
         dxcdgud(i) = vsigma(2)
         dxcdgd2(i) = vsigma(3)
      end do
#endif
   end subroutine


subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad,ex_coef)

implicit none
! arguments
integer, intent(in) :: xctype(3)
character(512), intent(out) :: xcdescr
integer, intent(out) :: xcspin
integer, intent(out) :: xcgrad
real(dp), intent(out) :: ex_coef
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
    case(XC_FAMILY_GGA)
      call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)
      call xc_f90_info_name(info,name)
      call xc_f90_func_end(p)
! post-processed gradients required
      xcgrad=2
    case(XC_FAMILY_HYB_GGA)
      call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)
      call xc_f90_info_name(info,name)
! get mixing coefficient for exchange
      call xc_f90_hyb_exx_coef(p, ex_coef) 
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
