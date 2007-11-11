
! Copyright (C) 2002-2007 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!!$subroutine genpptsq(vql,reducep,ngridp,vploff,nppt,ipmap,ivp,vpl,vpc,wppt)
!!$  use modmain
!!$  use modtddft
!!$  implicit none
!!$  ! arguments
!!$  real(8), intent(in) :: vql(3)
!!$  logical, intent(in) :: reducep
!!$  integer, intent(in) :: ngridp(3)
!!$  real(8), intent(in) :: vploff(3)
!!$  integer, intent(out) :: nppt
!!$  integer, intent(out) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
!!$  integer, intent(out) :: ivp(3,ngridp(1)*ngridp(2)*ngridp(3))
!!$  real(8), intent(out) :: vpl(3,ngridp(1)*ngridp(2)*ngridp(3))
!!$  real(8), intent(out) :: vpc(3,ngridp(1)*ngridp(2)*ngridp(3))
!!$  real(8), intent(out) :: wppt(ngridp(1)*ngridp(2)*ngridp(3))
!!$  ! local variables
!!$  integer :: nsymcrys_,lsplsymc_(maxsymcrys),lsplsymct(maxsymcrys)
!!$  integer :: lspl,isym,jsym
!!$
!!$  ! use symmetries of small group of q
!!$  call findgroupq(vql,epslat,bvec,binv,symlat,nsymcrys,lsplsymc,&
!!$       nsymcrysq,scqmap,ivscwrapq)
!!$  
!!$  ! save global variables
!!$  nsymcrys_=nsymcrys; lsplsymc_(:)=lsplsymc(:)
!!$
!!$  ! set pointer to point group elements
!!$  lsplsymct(:)=0
!!$  jsym=0
!!$  do isym=1,nsymcrysq
!!$     jsym=jsym+1
!!$     lsplsymct(jsym)=lsplsymc(scqmap(isym))
!!$  end do
!!$
!!$  ! update global variables
!!$  nsymcrys=nsymcrysq; lsplsymc(:)=lsplsymct(:)
!!$
!!$  call genppts(reducep,ngridp,vploff,nppt,ipmap,ivp,vpl,vpc,wppt)
!!$
!!$  ! restore global varialbes
!!$  nsymcrys=nsymcrys_; lsplsymc(:)=lsplsymc_(:)
!!$
!!$end subroutine genpptsq
