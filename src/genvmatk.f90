
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,wfmt,wfir,vmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot)
real(8), intent(in) :: vir(ngrtot)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngrtot,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer is,ia,ias,nrc,irc
integer ispn,ist,jst
real(8) t1
complex(8) zt1
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfmt(:,:)
complex(8), allocatable :: zfir(:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
! allocate local arrays
allocate(rfir(ngrtot))
allocate(zfmt(lmmaxvr,nrcmtmax))
allocate(zfir(ngrtot))
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do jst=1,nstsv
  do is=1,nspecies
    nrc=nrcmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ispn=1,nspinor
! apply potential to wavefunction
        do irc=1,nrc
          zfmt(:,irc)=vmt(:,irc,ias)*wfmt(:,irc,ias,ispn,jst)
        end do
        do ist=1,jst
! compute inner product (functions are in spherical coordinates)
          zt1=zfmtinp(.false.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
           wfmt(:,:,ias,ispn,ist),zfmt)
          vmat(ist,jst)=vmat(ist,jst)+zt1
        end do
      end do
    end do
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
rfir(:)=vir(:)*cfunir(:)
t1=omega/dble(ngrtot)
do jst=1,nstsv
  do ispn=1,nspinor
! apply potential to wavefunction
    zfir(:)=rfir(:)*wfir(:,ispn,jst)
    do ist=1,jst
      zt1=zdotc(ngrtot,wfir(:,ispn,ist),1,zfir,1)
      vmat(ist,jst)=vmat(ist,jst)+t1*zt1
    end do
  end do
end do
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
deallocate(rfir,zfmt,zfir)
return
end subroutine

