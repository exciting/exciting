


! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine xslinopt(iq)
  use modmain
use modinput
  use modxs
  use modtetra
  use modmpi
  use m_genwgrid
  use m_pade
  use m_genloss
  use m_gensigma
  use m_gensumrls
  use m_writeeps
  use m_writeloss
  use m_writesigma
  use m_writesumrls
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='xslinopt'
  character(256) :: filnam
  complex(8),allocatable :: mdf(:), mdf1(:),mdf2(:,:,:),w(:),wr(:),sigma(:)
  real(8),allocatable :: wplot(:),loss(:)
  real(8),allocatable :: eps1(:),eps2(:),cf(:,:)
  real(8) :: sumrls(3),brd
  integer :: n,m,recl,iw,wi,wf,nwdfp,nc,oct1,oct2,octl,octu,optcompt(3),i,j
  logical :: tq0
  logical, external :: tqgamma
  tq0=tqgamma(iq)
  ! number of components (3 for q=0)
  nc=1
  if (tq0) nc=3
  ! limits for w-points
  wi=wpari
  wf=wparf
  nwdfp=wparf-wpari+1
  ! matrix size for local field effects
  n=ngq(iq)
  allocate(mdf1(nwdf), mdf2(3, 3, nwdf), w(nwdf), wr(input%properties%dos%nwdos), wplot(input%properties%dos%nwdos), &
    &mdf(input%properties%dos%nwdos), &
       loss(input%properties%dos%nwdos), sigma(input%properties%dos%nwdos), cf(3, input%properties%dos%nwdos))
  allocate(eps1(input%properties%dos%nwdos), eps2(input%properties%dos%nwdos))
  mdf2(:,:,:)=zzero
  ! generate energy grids
  brd=0.d0
  if (input%xs%tddft%acont) brd=input%xs%broad
  call genwgrid(nwdf, wdos, input%xs%tddft%acont, 0.d0, w_cmplx=w)
  call genwgrid(input%properties%dos%nwdos, wdos, .false., brd, w_cmplx=wr)
  wplot=dble(wr)
  ! record length
  inquire(iolength=recl) mdf1(1)
  call getunit(unit1)
  ! neglect/include local field effects
  do m=1,n,max(n-1,1)
     ! loop over longitudinal components for optics
     do oct1=1,nc
	if (input%xs%dfoffdiag) then
           octl=1
           octu=nc
        else
           octl=oct1
           octu=oct1
        end if
        do oct2=octl,octu
           ! file name for inverse of dielectric function
           call genfilname(basename='IDF',asc=.false.,bzsampl=bzsampl,&
		acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),&
		&fxctype =input%xs%tddft%fxctypenumber, &
                tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=filnam)
           ! read macroscopic dielectric function (original frequencies)
           open(unit1,file=trim(filnam),form='unformatted', &
                action='read',status='old',access='direct',recl=recl)
           do iw=1,nwdf
              read(unit1,rec=iw) mdf1(iw)
           end do
           close(unit1)
           ! analytic continuation
	   if (input%xs%tddft%acont) then
	      call pade(input%properties%dos%nwdos, wr, nwdf, w, mdf1, mdf)
           else
              mdf(:)=mdf1(:)
           end if
           mdf2(oct1,oct2,:)=mdf(:)
        end do
     end do
     do oct1=1,nc
	if (input%xs%dfoffdiag) then
           octl=1
           octu=nc
        else
           octl=oct1
           octu=oct1
        end if
        do oct2=octl,octu
           optcompt(:)=(/oct1,oct2,0/)
           ! symmetrize the macroscopic dielectric function tensor
           call symt2app(oct1,oct2,nwdf,symt2,mdf2, mdf)
           ! file names for spectra
           call genfilname(basename='EPSILON',asc=.false.,bzsampl=bzsampl,&
		acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),&
        fxctype=input%xs%tddft%fxctypenumber, &
                tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fneps)
           call genfilname(basename='LOSS',asc=.false.,bzsampl=bzsampl,&
		acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),&
        fxctype =input%xs%tddft%fxctypenumber, &
                tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fnloss)
           call genfilname(basename='SIGMA',asc=.false.,bzsampl=bzsampl,&
		acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),&
       fxctype =input%xs%tddft%fxctypenumber, &
                tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fnsigma)
           call genfilname(basename='SUMRULES',asc=.false.,bzsampl=bzsampl,&
		acont = input%xs%tddft%acont, nar = .not.input%xs%tddft%aresdf, nlf = (m == 1),&
        fxctype=input%xs%tddft%fxctypenumber, &
                tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fnsumrules)
           ! generate optical functions
           call genloss(mdf,loss)
           call gensigma(dble(wr),mdf,optcompt,sigma)
           call gensumrls(dble(wr),mdf,sumrls)
           ! write optical functions to file
           call writeeps(iq,oct1,oct2,wplot,mdf,trim(fneps))
           call writeloss(iq,wplot,loss,trim(fnloss))
           call writesigma(iq,wplot,sigma,trim(fnsigma))
           call writesumrls(iq,sumrls,trim(fnsumrules))
           ! end loop over optical components
        end do
     end do
  end do ! m
  ! deallocate
  deallocate(mdf,mdf1,mdf2,w,wr,wplot,loss,sigma)
  deallocate(eps1,eps2,cf)
end subroutine xslinopt

!///////////////////////////////////////////////////////////////////////////////

subroutine symt2app(oct1,oct2,n,symt2,t,tsym)
  implicit none
  ! arguments
  integer, intent(in) :: oct1,oct2,n
  real(8), intent(in) :: symt2(3,3,3,3)
  complex(8), intent(in) :: t(3,3,n)
  complex(8), intent(out) :: tsym(n)
  ! local variables
  integer :: i,j
  ! symmetrize the macroscopic dielectric function tensor
  tsym(:)=(0.d0,0.d0)
  do i=1,3
     do j=1,3
        tsym(:)=tsym(:)+symt2(oct1,oct2,i,j)*t(i,j,:)
     end do
  end do
end subroutine
