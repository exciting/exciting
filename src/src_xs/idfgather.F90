
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine idfgather
  use modmain
  use modxs
  use modtetra
  use modmpi
  use m_filedel
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='idfgather'
  character(256) :: filnam,filnam2
  integer :: n,m,iq,iw,iproc,recl,nc,oct,oct1,oct2
  logical :: tq0
  complex(8), allocatable :: mdf1(:)
  logical, external :: tqgamma
  integer, external :: octmap
  allocate(mdf1(nwdf))
  inquire(iolength=recl) mdf1(1)
  call getunit(unit1)
  ! loop over q-points
  do iq=1,nqpt
     tq0=tqgamma(iq)
     ! number of components (3 for q=0)
     nc=1
     if (tq0) nc=3
     ! matrix size for local field effects
     n=ngq(iq)
     ! calculate k+q and G+k+q related variables
     call init1xs(qvkloff(1,iq))
     do m=1,n,max(n-1,1)
        ! loop over longitudinal components for optics
        do oct1=1,nc
        do oct2=1,nc
           oct=octmap(oct1,oct2)
           do iproc=0,procs-1
              wpari=firstofset(iproc,nwdf)
              wparf=lastofset(iproc,nwdf)
              ! filename for proc
              call genfilname(basename='IDF',bzsampl=bzsampl,acont=acont,&
                   nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,tq0=tq0,oc1=oct1,&
                   oc2=oct2,iqmt=iq,procs=procs,rank=iproc,filnam=filnam2)
              open(unit1,file=trim(filnam2),form='unformatted',&
                   action='read',status='old',access='direct',recl=recl)
              do iw=wpari,wparf
                 read(unit1,rec=iw-wpari+1) mdf1(iw)
              end do
              close(unit1)
           end do ! iproc
           ! write to file
           call genfilname(basename='IDF',bzsampl=bzsampl,&
                acont=acont,nar=.not.aresdf,nlf=(m==1),fxctype=fxctype,&
                tq0=tq0,oc1=oct1,oc2=oct2,iqmt=iq,filnam=filnam)
           open(unit1,file=trim(filnam),form='unformatted', &
                action='write',status='replace',access='direct',recl=recl)
           do iw=1,nwdf
              write(unit1,rec=iw) mdf1(iw)
           end do
           close(unit1)
           ! remove partial files
           do iproc=0,procs-1
              call genfilname(basename='IDF',bzsampl=bzsampl,acont=acont,&
                   nar=.not.aresdf,nlf=(m.eq.1),fxctype=fxctype,tq0=tq0, &
                   oc1=oct1,oc2=oct2,iqmt=iq,procs=procs,rank=iproc, &
                   filnam=filnam2)
              call filedel(trim(filnam2))
           end do
        end do ! oct
        end do
     end do ! m
     write(unitout,'(a,i8)') 'Info('//thisnam//'): inverse dielectric &
          &function gathered for q-point:',iq
  end do
  deallocate(mdf1)
end subroutine idfgather
