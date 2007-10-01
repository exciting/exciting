
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair
! !INTERFACE:
subroutine ggair(grhoir,gupir,gdnir,g2upir,g2dnir,g3rhoir,g3upir,g3dnir)
! !INPUT/OUTPUT PARAMETERS:
!   grhoir  : |grad rho| (out,real(ngrtot))
!   gupir   : |grad rhoup| (out,real(ngrtot))
!   gdnir   : |grad rhodn| (out,real(ngrtot))
!   g2upir  : grad^2 rhoup (out,real(ngrtot))
!   g2dnir  : grad^2 rhodn (out,real(ngrtot))
!   g3rhoir : (grad rho).(grad |grad rho|) (out,real(ngrtot))
!   g3upir  : (grad rhoup).(grad |grad rhoup|) (out,real(ngrtot))
!   g3dnir  : (grad rhodn).(grad |grad rhodn|) (out,real(ngrtot))
! !DESCRIPTION:
!   Computes $|\nabla\rho|$, $|\nabla\rho^{\uparrow}|$,
!   $|\nabla\rho^{\downarrow}|$, $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho\cdot(\nabla|\nabla\rho|)$,
!   $\nabla\rho^{\uparrow}\cdot(\nabla|\nabla\rho^{\uparrow}|)$ and
!   $\nabla\rho^{\downarrow}\cdot(\nabla|\nabla\rho^{\downarrow}|)$ for the
!   interstitial charge density, as required by the generalised gradient
!   approximation functional for spin-polarised densities. In the case of spin
!   unpolarised calculations, $|\nabla\rho|$, $\nabla^2\rho$ and
!   $\nabla\rho\cdot(\nabla|\nabla\rho|)$ are returned in the arrays
!   {\tt gupir}, {\tt g2upir} and {\tt g3upir}, respectively, while
!   {\tt grhoir}, {\tt gdnir}, {\tt g2dnir}, {\tt g3rhoir} and {\tt g3dnir} are
!   not referenced. See routines {\tt potxc} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
real(8), intent(out) :: grhoir(ngrtot)
real(8), intent(out) :: gupir(ngrtot)
real(8), intent(out) :: gdnir(ngrtot)
real(8), intent(out) :: g2upir(ngrtot)
real(8), intent(out) :: g2dnir(ngrtot)
real(8), intent(out) :: g3rhoir(ngrtot)
real(8), intent(out) :: g3upir(ngrtot)
real(8), intent(out) :: g3dnir(ngrtot)
! local variables
integer i,ig,ifg,ir
! allocatable arrays
real(8), allocatable :: rfir1(:,:)
real(8), allocatable :: rfir2(:,:)
complex(8), allocatable :: zfft1(:)
complex(8), allocatable :: zfft2(:)
allocate(rfir1(ngrtot,3))
allocate(rfir2(ngrtot,3))
allocate(zfft1(ngrtot))
allocate(zfft2(ngrtot))
if (spinpol) then
! rhoup for spin-polarised case
  zfft1(:)=0.5d0*(rhoir(:)+magir(:,ndmag))
else
! rho for spin-unpolarised case
  zfft1(:)=rhoir(:)
end if
call zfftifc(3,ngrid,-1,zfft1)
! |grad rhoup|
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
  end do
  call zfftifc(3,ngrid,1,zfft2)
  rfir1(:,i)=dble(zfft2(:))
end do
do ir=1,ngrtot
  gupir(ir)=sqrt(rfir1(ir,1)**2+rfir1(ir,2)**2+rfir1(ir,3)**2)
end do
! grad^2 rhoup
zfft2(:)=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3,ngrid,1,zfft2)
g2upir(:)=dble(zfft2(:))
! (grad rhoup).(grad |grad rhoup|)
zfft1(:)=gupir(:)
call zfftifc(3,ngrid,-1,zfft1)
g3upir(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
  end do
  call zfftifc(3,ngrid,1,zfft2)
  do ir=1,ngrtot
    g3upir(ir)=g3upir(ir)+rfir1(ir,i)*dble(zfft2(ir))
  end do
end do
if (spinpol) then
! rhodn
  zfft1(:)=0.5d0*(rhoir(:)-magir(:,ndmag))
  call zfftifc(3,ngrid,-1,zfft1)
! |grad rhodn|
  do i=1,3
    zfft2(:)=0.d0
    do ig=1,ngvec
      ifg=igfft(ig)
      zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
    end do
    call zfftifc(3,ngrid,1,zfft2)
    rfir2(:,i)=dble(zfft2(:))
  end do
  do ir=1,ngrtot
    gdnir(ir)=sqrt(rfir2(ir,1)**2+rfir2(ir,2)**2+rfir2(ir,3)**2)
  end do
! grad^2 rhodn
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
  end do
  call zfftifc(3,ngrid,1,zfft2)
  g2dnir(:)=dble(zfft2(:))
! (grad rhodn).(grad |grad rhodn|)
  zfft1(:)=gdnir(:)
  call zfftifc(3,ngrid,-1,zfft1)
  g3dnir(:)=0.d0
  do i=1,3
    zfft2(:)=0.d0
    do ig=1,ngvec
      ifg=igfft(ig)
      zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
    end do
    call zfftifc(3,ngrid,1,zfft2)
    do ir=1,ngrtot
      g3dnir(ir)=g3dnir(ir)+rfir2(ir,i)*dble(zfft2(ir))
    end do
  end do
! |grad rho|
  do ir=1,ngrtot
    grhoir(ir)=sqrt((rfir1(ir,1)+rfir2(ir,1))**2 &
                   +(rfir1(ir,2)+rfir2(ir,2))**2 &
                   +(rfir1(ir,3)+rfir2(ir,3))**2)
  end do
! (grad rho).(grad |grad rho|)
  zfft1(:)=grhoir(:)
  call zfftifc(3,ngrid,-1,zfft1)
  g3rhoir(:)=0.d0
  do i=1,3
    zfft2(:)=0.d0
    do ig=1,ngvec
      ifg=igfft(ig)
      zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
    end do
    call zfftifc(3,ngrid,1,zfft2)
    do ir=1,ngrtot
      g3rhoir(ir)=g3rhoir(ir)+(rfir1(ir,i)+rfir2(ir,i))*dble(zfft2(ir))
    end do
  end do
end if
deallocate(rfir1,rfir2,zfft1,zfft2)
return
end subroutine
!EOC

