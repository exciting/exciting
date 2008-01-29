
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine idfq(iq)
  use modmain
  use modxs
  use modfxcifc
  use modtetra
  use modmpi
  use m_genwgrid
  use m_dyson
  use m_getx0
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(*), parameter :: thisnam='idfq'
  character(256) :: filnam,filnam2
  complex(8),allocatable :: chi0(:,:), fxc(:,:), idf(:,:), mdf1(:),w(:)
  complex(8),allocatable :: chi0hd(:),chi0wg(:,:,:),chi0h(:)
  integer :: n,m,recl,j,iw,wi,wf,nwdfp,nc,oct
  logical :: tq0
  integer, external :: l2int
  logical, external :: tqgamma
  ! sampling type for Brillouin zone sampling
  bzsampl=l2int(tetra)
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
  allocate(chi0(n,n),fxc(n,n),idf(n,n),w(nwdf),mdf1(nwdf),chi0hd(nwdf))
  allocate(chi0wg(n,2,3),chi0h(3))
  fxc=zzero
  ! generate energy grid
  call genwgrid(nwdf,wdos,acont,0.d0,w_cmplx=w)
  ! filename for response function file
  call genfilname(basename='X0',asc=.false.,bzsampl=bzsampl,&
       acont=acont,nar=.not.aresdf,iq=iq,filnam=filnam)
  ! record length
  inquire(iolength=recl) mdf1(1)
  call getunit(unit1)
  call getunit(unit2)
  ! neglect/include local field effects
  do m=1,n,max(n-1,1)
     ! The ALDA kernel does not depend on q in principle, but the G-mesh
     ! depends through its cutoff for G+q on q. It is independent of w.
     if (fxctype.eq.5) then
        call fxcifc(fxctype,iq=iq,ng=m,fxcg=fxc)
        ! add symmetrized Coulomb potential (is equal to unity matrix)
        forall(j=1:m) 
           fxc(j,j)=fxc(j,j)+1.d0
        end forall
     end if
     ! loop over longitudinal components for optics
     do oct=1,nc
        ! filename for output file
        call genfilname(basename='IDF',asc=.false.,bzsampl=bzsampl,&
             acont=acont,nar=.not.aresdf,nlf=(m.eq.1),fxctype=fxctype,&
             tq0=tq0,oc=oct,iq=iq,procs=procs,rank=rank,filnam=filnam2)
        open(unit1,file=trim(filnam2),form='unformatted', &
             action='write',access='direct',recl=recl)
        do iw=wi,wf
           ! read Kohn-Sham response function
           call getx0(tq0,iq,iw,trim(filnam),'',chi0,chi0wg,&
                chi0h)
           ! assign components to main matrix for q=0
           if (tq0) then
              ! head
              chi0(1,1)=chi0h(oct)
              ! wings
              if (m.gt.1) then
                 chi0(1,2:)=chi0wg(2:,1,oct)
                 chi0(2:,1)=chi0wg(2:,2,oct)
              end if
           end if
           ! generate xc-kernel
           if (fxctype.ne.5) then
              call fxcifc(fxctype,ng=m,w=w(iw),alrc=alphalrc,&
                   alrcd=alphalrcdyn,blrcd=betalrcdyn,fxcg=fxc)
              ! add symmetrized Coulomb potential (is equal to unity matrix)
              forall(j=1:m) 
                 fxc(j,j)=fxc(j,j)+1.d0
              end forall
              ! head of pure f_xc kernel
              if (m.eq.1) fxc0(iw,oct)=fxc(1,1)-1.d0
           end if
           ! solve Dyson's equation for the interacting response function
           call dyson(n,chi0,fxc,idf)
           ! symmetrized inverse dielectric function (add one)
           forall(j=1:m) 
              idf(j,j)=idf(j,j)+1.d0
           end forall
           ! Adler-Wiser treatment of macroscopic dielectric function
           mdf1(iw)=1.d0/idf(1,1)
           ! write macroscopic dielectric function to file
           write(unit1,rec=iw-wi+1) mdf1(iw)
        end do ! iw
        close(unit1)
     end do ! oct
  end do ! m
  ! deallocate
  deallocate(chi0,chi0wg,chi0h,fxc,idf,mdf1,w,chi0hd)
end subroutine idfq
