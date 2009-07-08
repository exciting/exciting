



! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot1d
! !INTERFACE:


subroutine plot1d(fname,  nf, lmax, ld, rfmt, rfir,plotdef)
! !USES:
use modinput
use modmain
use FoX_wxml
! !INPUT/OUTPUT PARAMETERS:
!   fnum1 : plot file name (character*)
!   nf    : number of functions (in,integer)
!   lmax  : maximum angular momentum (in,integer)
!   ld    : leading dimension (in,integer)
!   rfmt  : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir  : real intersitial function (in,real(ngrtot,nf))
!   plotdef: type(plot1d) defining plotregion
! !DESCRIPTION:
!   Produces a 1D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} along the lines connecting the vertices in the global array
!   {\tt vvlp1d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
character(len=*), intent(in) :: fname
integer, intent(in) :: nf
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld, nrmtmax, natmtot, nf)
real(8), intent(in) :: rfir(ngrtot, nf)
type(plot1d_type)::plotdef
! local variables
integer::i, ip, iv,fnum1=50,fnum2=51
real(8)::fmin, fmax, t1
character(128)::buffer,buffer1
  type(xmlf_t), save::xf
! allocatable arrays
real(8), allocatable :: fp(:, :)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*, *)
  write(*, '("Error(plot1d): invalid number of functions : ", I8)') nf
  write(*, *)
  stop
end if
buffer=fname//"1D.OUT"
open(fnum1, file=trim(buffer), action='WRITE', form='FORMATTED')
buffer=fname//"LINES.OUT"
open(fnum2, file=trim(buffer), action='WRITE', form='FORMATTED')
call xml_OpenFile (fname//"1d.xml", xf, replace=.true.,pretty_print=.true.)
call xml_NewElement(xf,"plot1d")
call xml_NewElement(xf,"title")
call xml_AddCharacters(xf,trim(input%title))
call xml_endElement(xf,"title")
! connect the plotting vertices
  nvp1d=size(plotdef%path%pointarray)
  npp1d=plotdef%path%steps
  allocate(fp(npp1d, nf))
  if (allocated(dvp1d)) deallocate(dvp1d)
  allocate(dvp1d(nvp1d))
  if (allocated(vplp1d)) deallocate(vplp1d)
  allocate(vplp1d(3, npp1d))
  if (allocated(dpp1d)) deallocate(dpp1d)
  allocate(dpp1d(npp1d))
call connect(input%structure%crystal%basevect, plotdef, &
size(plotdef%path%pointarray),plotdef%path%steps,vplp1d, dvp1d, dpp1d)
do i=1, nf
! evaluate function at each point
  call rfarray(lmax, ld, rfmt(:, :, :, i), rfir(:, i), npp1d, vplp1d, fp(:, i))
end do
fmin=fp(1, 1)
fmax=fp(1, 1)
do ip=1, npp1d
  do i=1, nf
    fmin=min(fmin, fp(ip, i))
    fmax=max(fmax, fp(ip, i))
  end do
! write the point distances and function to file
  write(fnum1, '(5G18.10)') dpp1d(ip), (fp(ip, i), i=1, nf)
  call xml_NewElement(xf,"point")
  write(buffer,'(5G18.10)')dpp1d(ip)
   call xml_AddAttribute(xf, "distance", trim(adjustl(buffer)))
   do i=1,nf
    write(buffer,'(5G18.10)')fp(ip, i)
     write(buffer1,*)i
  call xml_AddAttribute(xf, "function"//trim(adjustl(buffer1)), trim(adjustl(buffer)))
   end do
   call xml_endElement(xf,"point")
end do
! write the vertex location lines
t1=0.5d0*(fmax-fmin)
do iv=1, nvp1d
  write(fnum2, '(2G18.10)') dvp1d(iv), fmax+t1
  write(fnum2, '(2G18.10)') dvp1d(iv), fmin-t1
  write(fnum2, '("     ")')
end do
do iv=1,size(plotdef%path%pointarray)
call xml_NewElement(xf,"vertex")

 write(buffer,'(5G18.10)')dvp1d(iv)
call xml_addAttribute(xf,"distance",trim(adjustl(buffer)))
 write(buffer,'(5G18.10)')fmax+t1
call xml_addAttribute(xf,"upperboundary",trim(adjustl(buffer)))
 write(buffer,'(5G18.10)')fmin-t1
call xml_addAttribute(xf,"lowerboundary",trim(adjustl(buffer)))
call xml_addAttribute(xf,"label",trim(adjustl(plotdef%path%pointarray(iv)%point%label)))
 write(buffer,'(5G18.10)')plotdef%path%pointarray(iv)%point%coord
call xml_addAttribute(xf,"coord",trim(adjustl(buffer)))
call xml_endElement(xf,"vertex")
end do
call xml_close(xf)
close(fnum1)
close(fnum2)
deallocate(fp)
deallocate(dvp1d)
deallocate(vplp1d)
deallocate(dpp1d)
return
end subroutine
!EOC
