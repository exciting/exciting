



! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot3d
! !INTERFACE:


subroutine plot3d(fname, nf, lmax, ld, rfmt, rfir,plotdef)
! !USES:
use modinput
use modmain
use FoX_wxml
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngrtot,nf))
! !DESCRIPTION:
!   Produces a 3D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} in the parallelepiped defined by the corner vertices in the
!   global array {\tt vclp3d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Modified, October 2008 (F. Bultmark, F. Cricchio, L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
character(len=*), intent(in)  :: fname
integer, intent(in) :: nf
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld, nrmtmax, natmtot, nf)
real(8), intent(in) :: rfir(ngrtot, nf)
type(plot3d_type),intent(in)::plotdef
! local variables
integer::np, ip, ip1, ip2, ip3, i,fnum=50
real(8)::v1(3), v2(3), v3(3)
real(8)::t1, t2, t3
character(512)::buffer,buffer1
  type(xmlf_t), save::xf
! allocatable arrays
real(8), allocatable :: vpl(:, :)
real(8), allocatable :: fp(:, :)
buffer=fname//"3D.OUT"
open(fnum, file=trim(buffer), action='WRITE', form='FORMATTED')
call xml_OpenFile (fname//"3d.xml", xf, replace=.true.,pretty_print=.true.)
call xml_NewElement(xf,"plot3d")
call xml_NewElement(xf,"title")
call xml_AddCharacters(xf,trim(input%title))
call xml_endElement(xf,"title")
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*, *)
  write(*, '("Error(plot3d): invalid number of functions : ", I8)') nf
  write(*, *)
  stop
end if
! allocate local arrays
allocate(vpl(3, plotdef%box%grid(1)*plotdef%box%grid(2)*plotdef%box%grid(3)))
allocate(fp(plotdef%box%grid(1)*plotdef%box%grid(2)*plotdef%box%grid(3), nf))
! generate 3D grid
v1(:)=plotdef%box%pointarray(1)%point%coord-plotdef%box%origin%coord
v2(:)=plotdef%box%pointarray(2)%point%coord-plotdef%box%origin%coord
v3(:)=plotdef%box%pointarray(3)%point%coord-plotdef%box%origin%coord
ip=0
do ip3=0, plotdef%box%grid(3)-1
  t3=dble(ip3)/dble(plotdef%box%grid(3))
  do ip2=0, plotdef%box%grid(2)-1
    t2=dble(ip2)/dble(plotdef%box%grid(2))
    do ip1=0, plotdef%box%grid(1)-1
      t1=dble(ip1)/dble(plotdef%box%grid(1))
      ip=ip+1
      vpl(:, ip)=t1*v1(:)+t2*v2(:)+t3*v3(:)+plotdef%box%origin%coord
    end do
  end do
end do
np=ip
! evaluate the functions at the grid points
do i=1, nf
  call rfarray(lmax, ld, rfmt(:, :, :, i), rfir(:, i), np, vpl, fp(:, i))
end do
! write functions to file
  write(fnum, '(3I6, " : grid size")') plotdef%box%grid(:)
  write(buffer,'(3I6)') plotdef%box%grid(:)
   call xml_AddAttribute(xf, "grid", trim(adjustl(buffer)))
do ip=1, np
  call r3mv(input%structure%crystal%basevect, vpl(:, ip), v1)
  write(fnum, '(7G18.10)') v1(:), (fp(ip, i), i=1, nf)
    call xml_newElement(xf,"point")
    write(buffer,'(G18.10)') v1(1)
    call xml_AddAttribute(xf, "x", trim(adjustl(buffer)))
    write(buffer,'(G18.10)') v1(2)
    call xml_AddAttribute(xf, "y", trim(adjustl(buffer)))
    write(buffer,'(G18.10)') v1(3)
    call xml_AddAttribute(xf, "z", trim(adjustl(buffer)))
    do i=1,nf
       write(buffer,'(G18.10)')fp(ip, i)
       write(buffer1,*)i
       call xml_AddAttribute(xf, "function"//trim(adjustl(buffer1)), trim(adjustl(buffer)))
    end do
    call xml_endElement(xf,"point")

end do
deallocate(vpl, fp)
call xml_Close(xf)
close(fnum)
return
end subroutine
!EOC
