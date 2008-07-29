
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getscreen(iqr,ngq,scrh,scrw,scrb)
  use modtetra
  use m_genfilname
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iqr,ngq
  complex(8), intent(out) :: scrb(ngq,ngq),scrw(ngq,2,3),scrh(9)
  ! local variables
  character(256) :: fname
  real(8) :: rm(2,9)
  integer :: igq1,igq2,j,it1,it2,it3,un,bzsampl
  ! sampling of Brillouin zone
  bzsampl=0
  if (tetra) bzsampl=1
  ! read in screening
  call genfilname(basename='SCREEN',iq=iqr,bzsampl=bzsampl,filnam=fname)
  call getunit(un)
  open(un,file=trim(fname),form='formatted',action='read',status='old')
  do igq1=1,ngq
     do igq2=1,ngq
        if (iqr.eq.1) then
           if ((igq1.eq.1).and.(igq2.eq.1)) then
              read(un,*) (it1,it2,it3,rm(1,j),rm(2,j),j=1,9)
              scrh(:)=cmplx(rm(1,:),rm(2,:),8)
           end if
           if ((igq1.eq.1).and.(igq2.ne.1)) then
              read(un,*) (it1,it2,it3,rm(1,j),rm(2,j),j=1,3)
              scrw(igq2,1,:)=cmplx(rm(1,:3),rm(2,:3),8)
           end if
           if ((igq1.ne.1).and.(igq2.eq.1)) then
              read(un,*) (it1,it2,it3,rm(1,j),rm(2,j),j=1,3)
              scrw(igq1,2,:)=cmplx(rm(1,:3),rm(2,:3),8)
           end if
           if ((igq1.ne.1).and.(igq2.ne.1)) read(un,*) it1,it2,it3,&
                rm(1,1),rm(2,1)
           scrb(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
        else
           read(un,*) it1,it2,it3,rm(1,1),rm(2,1)
           scrb(igq1,igq2)=cmplx(rm(1,1),rm(2,1),8)
        end if
     end do
  end do
  close(un)
end subroutine getscreen
