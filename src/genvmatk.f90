
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,wfmt,wfir,vmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: vir(ngrtot)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngrtot,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer is,ia,ias,nr,ir,irc
integer ispn,ist1,ist2
real(8) t1
complex(8) zt1
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfmt(:,:)
complex(8), allocatable :: zfir(:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
! allocate local arrays
allocate(rfmt(lmmaxvr,nrcmtmax,natmtot))
allocate(rfir(ngrtot))
allocate(zfmt(lmmaxvr,nrcmtmax))
allocate(zfir(ngrtot))
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
! convert muffin-tin potential to spherical coordinates
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,vmt(1,ir,ias), &
       1,0.d0,rfmt(1,irc,ias),1)
    end do
  end do
end do
do ist2=1,nstsv
  do is=1,nspecies
    nr=nrcmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ispn=1,nspinor
! apply potential to wavefunction
        do irc=1,nr
          zfmt(:,irc)=rfmt(:,irc,ias)*wfmt(:,irc,ias,ispn,ist2)
        end do
        do ist1=1,ist2
! compute inner product (functions are in spherical coordinates)
          zt1=zfmtinp(lmaxvr,nr,rcmt(1,is),lmmaxvr,wfmt(1,1,ias,ispn,ist1),zfmt)
          vmat(ist1,ist2)=vmat(ist1,ist2)+zt1*fourpi/dble(lmmaxvr)
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
do ist2=1,nstsv
  do ispn=1,nspinor
! apply potential to wavefunction
    zfir(:)=rfir(:)*wfir(:,ispn,ist2)
    do ist1=1,ist2
      zt1=zdotc(ngrtot,wfir(1,ispn,ist1),1,zfir,1)
      vmat(ist1,ist2)=vmat(ist1,ist2)+t1*zt1
    end do
  end do
end do
! lower triangular part
do ist1=1,nstsv
  do ist2=1,ist1-1
    vmat(ist1,ist2)=conjg(vmat(ist2,ist1))
  end do
end do
deallocate(rfmt,rfir,zfmt,zfir)
return
end subroutine

