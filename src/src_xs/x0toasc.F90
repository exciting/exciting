
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine x0toasc
  use modmain
  use modxs
  use m_getx0
  use m_putx0
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='x0toasc'
  character(256) :: filnam, filnama
  integer :: n,iq,igq,igqp,iw,oct,noct,un
  complex(8) :: zt
  complex(8), allocatable :: chi0(:,:),chi0wg(:,:,:),chi0hd(:)
  logical :: tq0
  logical, external :: tqgamma
  if (calledxs.eq.1) call init0
  ! initialise universal variables
  call init1
  ! save Gamma-point variables
  call tdsave0
  ! initialize q-point set
  call init2xs
  ! loop over q-points
  do iq=1,nqpt
     tq0=tqgamma(iq)
     ! calculate k+q and G+k+q related variables
     call init1xs(qvkloff(1,iq))
     ! size of local field effects
     n=ngq(iq)
     ! allocate
     allocate(chi0(n,n),chi0wg(n,2,3),chi0hd(3))
     ! filenames
     call genfilname(asc=.true.,basename='X0',bzsampl=bzsampl,acont=acont,&
          nar=.not.aresdf,iq=iq,filnam=filnama)
     call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,&
          nar=.not.aresdf,iq=iq,filnam=filnam)
     ! open file to write ASCI
     call getunit(un)
     open(unit=un,file=trim(filnama),form='formatted', &
          action='write',status='replace')
     noct=1
     if (tq0) noct=3
     do oct=1,noct
        do iw=1,nwdf
           ! read from binary file
           call getx0(tq0,iq,iw,trim(filnam),'',chi0,chi0wg,chi0hd)
           ! write to ASCII file
           if (tq0) then
              ! head
              chi0(1,1)=chi0hd(oct)
              ! wings
              if (n.gt.1) then
                 chi0(1,2:)=chi0wg(2:,1,oct)
                 chi0(2:,1)=chi0wg(2:,2,oct)
              end if
           end if
           do igq=1,ngq(iq)
              do igqp=1,ngq(iq)
                 zt=chi0(igq,igqp)
                 write(un,'(5i6,3g18.10)') iq,oct,iw,igq,igqp,zt,abs(zt)
              end do
              write(un+1,'(100g12.4)') abs(chi0(igq,:))
              write(un+2,'(100g12.4)') dble(chi0(igq,:))
              write(un+3,'(100g12.4)') aimag(chi0(igq,:))
           end do
        end do
     end do
     ! close file
     close(un)
     deallocate(chi0,chi0wg,chi0hd)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sham response &
          &function converted to ASCII file for q-point:',iq
  end do
  call genfilname(setfilext=.true.)
end subroutine x0toasc
