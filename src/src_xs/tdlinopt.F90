
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdlinopt(iq)
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_findexciton
  use m_genwgrid
  use m_pade
  use m_genloss
  use m_gensigma
  use m_gensumrls
  use m_gensymdf
  use m_writeeps
  use m_writeloss
  use m_writesigma
  use m_writesumrls
  use m_writeexciton
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='tdlinopt'
  character(256) :: filnam,filnam2
  complex(8),allocatable :: mdf(:), mdf1(:),w(:),wr(:),sigma(:)
  real(8),allocatable :: wplot(:),loss(:)
  real(8),allocatable :: eps1(:),eps2(:),cf(:,:)
  real(8) :: sumrls(3),brd
  complex(8) :: zt1
  integer :: n,m,recl,iw,wi,wf,nwdfp,nc,oct,oct1,oct2,optcompt(3)
  logical :: tq0
  logical, external :: tqgamma
  integer, external :: octmap
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
  allocate(mdf1(nwdf),w(nwdf),wr(nwdos),wplot(nwdos),mdf(nwdos), &
       loss(nwdos),sigma(nwdos),cf(3,nwdos))
  allocate(eps1(nwdos),eps2(nwdos))
  ! generate energy grids
  brd=0.d0
  if (acont) brd=brdtd
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
  call genwgrid(nwdos,wdos,.false.,brd,w_cmplx=wr)
  wplot=dble(wr)
  ! record length
  inquire(iolength=recl) mdf1(1)
  call getunit(unit1)
  ! neglect/include local field effects
  do m=1,n,max(n-1,1)
     ! loop over longitudinal components for optics
     do oct1=1,nc
     do oct2=1,nc
        if (oct1.ne.oct2) goto 111
        oct=octmap(oct1,oct2)
        optcomp(1,1)=oct1
        optcomp(2,1)=oct2
        optcompt(:)=optcomp(:,1)
        ! symmetrization matrix for dielectric function
        call gensymdf(oct1,oct2)
        ! file name for inverse of dielectric function
        call genfilname(basename='IDF',asc=.false.,bzsampl=bzsampl,&
             acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
             tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=filnam2)
        ! read macroscopic dielectric function (original frequencies)
        open(unit1,file=trim(filnam2),form='unformatted', &
             action='read',status='old',access='direct',recl=recl)
        do iw=1,nwdf
           read(unit1,rec=iw) mdf1(iw)
        end do
        close(unit1)
        ! analytic continuation
        if (acont) then
           call pade(nwdos,wr,nwdf,w,mdf1,mdf)
        else
           mdf(:)=mdf1(:)
        end if
        ! file names for spectra
        call genfilname(basename='EPSILON',asc=.false.,bzsampl=bzsampl,&
             acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
             tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fneps)
        call genfilname(basename='LOSS',asc=.false.,bzsampl=bzsampl,&
             acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
             tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fnloss)
        call genfilname(basename='SIGMA',asc=.false.,bzsampl=bzsampl,&
             acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
             tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fnsigma)
        call genfilname(basename='SUMRULES',asc=.false.,bzsampl=bzsampl,&
             acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
             tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=fnsumrules)
        if (tetra.and.(fxctype/=0).and.(oct1.eq.oct2)) then
           call genfilname(basename='IDF',asc=.false.,bzsampl=bzsampl,&
                acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=0,&
                tq0=tq0,oc1=oct1,oc2=oct1,iqmt=iq,filnam=filnam)
           ! read macroscopic dielectric function (RPA)
           open(unit1,file=trim(filnam),form='unformatted', &
                action='read',status='old',access='direct',recl=recl)
           do iw=1,nwdf
              read(unit1,rec=iw) zt1
              mdfrpa(iw,oct1,1)=dble(zt1)
              mdfrpa(iw,oct1,2)=aimag(zt1)
           end do
           close(unit1)
           ! derivative of RPA dielectric function
           call fderiv(1,nwdos,wplot,mdfrpa(1,oct1,1),mdfrpad(1,oct1),cf)
           ! derivative of xc-kernel
           call fderiv(1,nwdos,wplot,fxc0(1,oct1),fxc0d(1,oct1),cf)
           ! file name for exciton file
           call genfilname(basename='EXCITON',asc=.false.,bzsampl=bzsampl,&
                acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
                tq0=tq0,oc1=oct1,oc2=oct1,iqmt=iq,filnam=fnexciton)
           ! determination of excitons in combination with tetrahedron method
           ! and neglect local field effects for kernel but consider for
           ! the RPA solution
           if (m.eq.max(n,1)) then
              call findexciton(oct1,nwdos,dble(w))
              call writeexciton(iq,oct1,wplot,mdf,trim(fnexciton))
           end if
        end if
        ! generate optical functions
        call genloss(mdf,loss)
        call gensigma(dble(wr),mdf,optcompt,sigma)
        call gensumrls(dble(wr),mdf,sumrls)
        ! write optical functions to file
        call writeeps(iq,wplot,mdf,trim(fneps))
        call writeloss(iq,wplot,loss,trim(fnloss))
        call writesigma(iq,wplot,sigma,trim(fnsigma))
        call writesumrls(iq,sumrls,trim(fnsumrules))
111     continue
     end do ! oct
     end do
  end do ! m
  ! deallocate
  if (allocated(mdfrpa)) deallocate(mdfrpa)
  deallocate(mdf,mdf1,w,wr,wplot,loss,sigma)
  deallocate(eps1,eps2,cf)
end subroutine tdlinopt
