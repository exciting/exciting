
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gencfun
! !INTERFACE:
subroutine gencfun
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the smooth characteristic function. This is the function which is
!   0 within the muffin-tins and 1 in the intersitial region and is constructed
!   from radial step function form-factors with $G<G_{\rm max}$. The form
!   factors are given by
!   $$ \tilde{\Theta}_i(G)=\begin{cases}
!    \frac{4\pi R_i^3}{3 \Omega} & G=0 \\
!    \frac{4\pi R_i^3}{\Omega}\frac{j_1(GR_i)}{GR_i} & 0<G\le G_{\rm max} \\
!    0 & G>G_{\rm max}\end{cases}, $$
!   where $R_i$ is the muffin-tin radius of the $i$th species and $\Omega$ is
!   the unit cell volume. Therefore the characteristic function in $G$-space is
!   $$ \tilde{\Theta}({\bf G})=\delta_{G,0}-\sum_{ij}\exp(-i{\bf G}\cdot
!    {\bf r}_{ij})\tilde{\Theta}_i(G), $$
!   where ${\bf r}_{ij}$ is the position of the $j$th atom of the $i$th species.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ig,ifg
real(8) t1,t2,jl(0:1)
! allocatable arrays
real(8), allocatable :: ff(:,:)
complex(8), allocatable :: zfft(:)
allocate(ff(ngvec,nspecies))
allocate(zfft(ngrtot))
! allocate global characteristic function arrays
if (allocated(cfunig)) deallocate(cfunig)
allocate(cfunig(ngvec))
if (allocated(cfunir)) deallocate(cfunir)
allocate(cfunir(ngrtot))
! step function form factors for each species
do is=1,nspecies
  t1=(fourpi/omega)*rmt(is)**3
  do ig=1,ngvec
    if (gc(ig).gt.epslat) then
      t2=gc(ig)*rmt(is)
      call sbessel(1,t2,jl)
      ff(ig,is)=t1*jl(1)/t2
    else
      ff(ig,is)=t1/3.d0
    end if
  end do
end do
! set up the characteristic function in G-space with G < gmaxvr
cfunig(:)=0.d0
cfunig(1)=1.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=1,ngvec
      cfunig(ig)=cfunig(ig)-conjg(sfacg(ig,ias))*ff(ig,is)
    end do
  end do
end do
zfft(:)=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  zfft(ifg)=cfunig(ig)
end do
! Fourier transform to real-space
call zfftifc(3,ngrid,1,zfft)
cfunir(:)=dble(zfft(:))
deallocate(ff,zfft)
return
end subroutine
!EOC

