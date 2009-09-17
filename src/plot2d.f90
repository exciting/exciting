



! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot2d
! !INTERFACE:


subroutine plot2d(fname, nf, lmax, ld, rfmt, rfir,plotdef)
! !USES:
use modinput
use modmain
use FoX_wxml
! !INPUT/OUTPUT PARAMETERS:
!   fname : plot file name character(len=*)
!   nf   : number of functions (in,integer)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngrtot,nf))
!   plotdef:type(plot2d) defines plot region
! !DESCRIPTION:
!   Produces a 2D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} on the parallelogram defined by the corner vertices in the global
!   array {\tt vclp2d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
character(len=*) , intent(in):: fname
integer, intent(in) :: nf
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld, nrmtmax, natmtot, nf)
real(8), intent(in) :: rfir(ngrtot, nf)
type(plot2d_type)::plotdef
! local variables
integer::i, ip, ip1, ip2,fnum=50
real(8)::vl1(3), vl2(3), vc1(3), vc2(3)
real(8)::d1, d2, d12, t1, t2, t3, t4
character(128)::buffer,buffer1
  type(xmlf_t), save::xf
! allocatable arrays
real(8), allocatable :: vpl(:, :)
real(8), allocatable :: fp(:, :)
buffer=fname//"2D.OUT"
open(fnum, file=trim(buffer), action='WRITE', form='FORMATTED')
call xml_OpenFile (fname//"2d.xml", xf, replace=.true.,pretty_print=.true.)
call xml_NewElement(xf,"plot2d")

if ((nf.lt.1).or.(nf.gt.4)) then
  write(*, *)
  write(*, '("Error(plot2d): invalid number of functions : ", I8)') nf
  write(*, *)
  stop
end if
! allocate local arrays
allocate(vpl(3, plotdef%parallelogram%grid(1)*plotdef%parallelogram%grid(2)))
allocate(fp(plotdef%parallelogram%grid(1)*plotdef%parallelogram%grid(2), nf))
! generate 2D grid
vl1(:)=plotdef%parallelogram%pointarray(1)%point%coord-plotdef%parallelogram%origin%coord
vl2(:)=plotdef%parallelogram%pointarray(2)%point%coord-plotdef%parallelogram%origin%coord
vc1(:) = vl1(1) * input%structure%crystal%basevect(:, 1) + vl1(2) * input%structure%crystal%basevect(:, 2) + vl1(3) *&
    &input%structure%crystal%basevect(:, 3)
vc2(:) = vl2(1) * input%structure%crystal%basevect(:, 1) + vl2(2) * input%structure%crystal%basevect(:, 2) + vl2(3) *&
    &input%structure%crystal%basevect(:, 3)
d1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
d2=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
if ((d1.lt.input%structure%epslat).or.(d2.lt.input%structure%epslat)) then
  write(*, *)
  write(*, '("Error(plot2d): zero length plotting vectors")')
  write(*, *)
  stop
end if
d12=(vc1(1)*vc2(1)+vc1(2)*vc2(2)+vc1(3)*vc2(3))/(d1*d2)
ip=0
do ip2=0, plotdef%parallelogram%grid(2)-1
  do ip1=0, plotdef%parallelogram%grid(1)-1
    ip=ip+1
    t1=dble(ip1)/dble(plotdef%parallelogram%grid(1))
    t2=dble(ip2)/dble(plotdef%parallelogram%grid(2))
    vpl(:, ip)=t1*vl1(:)+t2*vl2(:)+vclp2d(:, 1)
  end do
end do
! evaluate the functions at the grid points
do i=1, nf
  call rfarray(lmax, ld, rfmt(:, :, :, i), rfir(:, i), ip, vpl, fp(:, i))
end do
! write the functions to file
write(fnum, '(2I6, " : grid size")') plotdef%parallelogram%grid(:)
  write(buffer,'(2I6)') plotdef%parallelogram%grid(:)
   call xml_AddAttribute(xf, "grid", trim(adjustl(buffer)))
call xml_NewElement(xf,"title")
call xml_AddCharacters(xf,trim(input%title))
call xml_endElement(xf,"title")
ip=0
do ip2=0, plotdef%parallelogram%grid(2)-1
  do ip1=0, plotdef%parallelogram%grid(1)-1
    ip=ip+1
    t1=dble(ip1)/dble(plotdef%parallelogram%grid(1))
    t2=dble(ip2)/dble(plotdef%parallelogram%grid(2))
    t3=t1*d1+t2*d2*d12
    t4=t2*d2*sqrt(abs(1.d0-d12**2))
    write(fnum, '(6G18.10)') t3, t4, (fp(ip, i), i=1, nf)
    call xml_newElement(xf,"point")
    write(buffer,'(6G18.10)') t3
    call xml_AddAttribute(xf, "x", trim(adjustl(buffer)))
    write(buffer,'(6G18.10)') t4
    call xml_AddAttribute(xf, "y", trim(adjustl(buffer)))
    do i=1,nf
       write(buffer,'(5G18.10)')fp(ip, i)
       write(buffer1,*)i
       call xml_AddAttribute(xf, "function"//trim(adjustl(buffer1)), trim(adjustl(buffer)))
    end do
    call xml_endElement(xf,"point")
  end do
end do
close(fnum)
call xml_Close(xf)
deallocate(vpl, fp)
return
end subroutine
!EOC
