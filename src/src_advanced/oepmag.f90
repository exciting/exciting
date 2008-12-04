
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepmag(tsh,wfmt1,wfmt2,wfir1,wfir2,zmagmt,zmagir)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh
complex(8), intent(in) ::  wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfir1(ngrtot,nspinor)
complex(8), intent(in) ::  wfir2(ngrtot,nspinor)
complex(8), intent(out) :: zmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
complex(8), intent(out) :: zmagir(ngrtot,ndmag)
! local variables
integer is,ia,ias,nrc,ir,idm
complex(8) zt1,zt2
! allocatable arrays
complex(8), allocatable :: zvfmt(:,:,:)
allocate(zvfmt(lmmaxvr,nrcmtmax,ndmag))
! muffin-tin part
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call oepmagmt(tsh,is,wfmt1(:,:,ias,1),wfmt1(:,:,ias,2),wfmt2(:,:,ias,1), &
     wfmt2(:,:,ias,2),zvfmt)
    do idm=1,ndmag
      zmagmt(:,1:nrc,ias,idm)=zvfmt(:,1:nrc,idm)
    end do
  end do
end do
! interstitial part
do ir=1,ngrtot
! calculate the z-component of mangetisation: up-up - dn-dn
  zmagir(ir,ndmag)=conjg(wfir1(ir,1))*wfir2(ir,1)-conjg(wfir1(ir,2))*wfir2(ir,2)
  if (ncmag) then
! up-dn spin density
    zt1=conjg(wfir1(ir,1))*wfir2(ir,2)
! dn-up spin density
    zt2=conjg(wfir1(ir,2))*wfir2(ir,1)
! calculate the x-component: up-dn + dn-up
    zmagir(ir,1)=zt1+zt2
! calculate the y-component: i*(dn-up - up-dn)
    zmagir(ir,2)=zi*(zt2-zt1)
  end if
end do
deallocate(zvfmt)
return
end subroutine

