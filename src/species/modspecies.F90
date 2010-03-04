
! Copyright (C) 2005-2010 S. Sagmeister, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module   modspecies

character(*), parameter :: version='0.9.224'
integer, parameter :: maxapword=4
character(256) :: bname,str,apwdescr,apwtype,suffix,inspecies
integer :: apword,apwdm(maxapword),nlorbsc
integer :: apwordx,apwdmx(maxapword)
logical :: locorb,locorbsc,searchlocorb,apwvex(maxapword),apwve(maxapword)
logical :: fullsearchlocorbsc

! order of predictor-corrector polynomial
integer, parameter :: np=4
! maximum number of states
integer, parameter :: maxspst=40
! maximum angular momentum allowed
integer, parameter :: lmax=50
! exchange-correlation type
integer, parameter :: xctype=3
integer, parameter :: xcgrad=0
integer nz,spnst,spnr,nrmt
integer nlx,nlorb,i,l,maxl,io
integer ist,jst,ir,iostat
real(8), parameter :: pi=3.1415926535897932385d0
! core-valence cut-off energy
real(8) :: ecvcut
! semi-core cut-off energy
real(8) :: esccut
! band offset energy
real(8), parameter :: boe=0.15d0
real(8) spmass,rmt,spzn,sprmin,sprmax,t1
character(256) spsymb,spname


! automatic arrays
logical spcore(maxspst)
integer spn(maxspst),spl(maxspst),spk(maxspst)
real(8) spocc(maxspst),eval(maxspst)
! allocatable arrays
real(8), allocatable :: r(:),rho(:),vr(:),rwf(:,:,:)
real(8), allocatable :: fr(:),gr(:),cf(:,:)

end module
