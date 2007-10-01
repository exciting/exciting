
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepmagmt(is,wfmt1,wfmt2,wfmt3,wfmt4,zvfmt)
use modmain
implicit none
! arguments
integer,  intent(in) :: is
complex(8), intent(in) ::  wfmt1(lmmaxvr,nrcmtmax)
complex(8), intent(in) ::  wfmt2(lmmaxvr,nrcmtmax)
complex(8), intent(in) ::  wfmt3(lmmaxvr,nrcmtmax)
complex(8), intent(in) ::  wfmt4(lmmaxvr,nrcmtmax)
complex(8), intent(out) :: zvfmt(lmmaxvr,nrcmtmax,ndmag)
! local variables
integer nr
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:)
allocate(zfmt(lmmaxvr,nrcmtmax,2))
! muffin-tin part
nr=nrcmt(is)
! up-up spin density
call vnlrhomt(is,wfmt1,wfmt3,zfmt(1,1,1))
! dn-dn spin density
call vnlrhomt(is,wfmt2,wfmt4,zfmt(1,1,2))
! calculate the z-component of mangetisation: up-up - dn-dn
zvfmt(:,1:nr,ndmag)=zfmt(:,1:nr,1)-zfmt(:,1:nr,2)
! non-collinear case
if (ndmag.eq.3) then
! up-dn spin density
  call vnlrhomt(is,wfmt1,wfmt4,zfmt(1,1,1))
! dn-up spin density
  call vnlrhomt(is,wfmt2,wfmt3,zfmt(1,1,2))
! calculate the x-component: up-dn + dn-up
  zvfmt(:,1:nr,1)=zfmt(:,1:nr,1)+zfmt(:,1:nr,2)
! calculate the y-component: i*(dn-up - up-dn)
  zvfmt(:,1:nr,2)=zi*(zfmt(:,1:nr,2)-zfmt(:,1:nr,1))
end if
deallocate(zfmt)
return
end subroutine

