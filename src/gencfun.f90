!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gencfun
! !INTERFACE:
!
!
Subroutine gencfun
! !USES:
      Use modinput
      Use modmain
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
      Implicit None
! local variables
      Integer :: is, ia, ig, ifg
      Real (8) :: t1, t2
      Complex (8) zt1
! allocatable arrays
      Real (8), Allocatable :: ffacg (:)
      Complex (8), Allocatable :: zfft (:)
      Allocate (ffacg(ngrtot))
      Allocate (zfft(ngrtot))
! allocate global characteristic function arrays
      If (allocated(cfunig)) deallocate (cfunig)
      Allocate (cfunig(ngrtot))
      If (allocated(cfunir)) deallocate (cfunir)
      Allocate (cfunir(ngrtot))
      cfunig (:) = 0.d0
      cfunig (1) = 1.d0
      t1 = fourpi / omega
! begin loop over species
      Do is = 1, nspecies
! smooth step function form factors for each species
         Do ig = 1, ngrtot
            If (gc(ig) .Gt. input%structure%epslat) Then
               t2 = gc (ig) * rmt (is)
               ffacg (ig) = t1 * (Sin(t2)-t2*Cos(t2)) / (gc(ig)**3)
            Else
               ffacg (ig) = (t1/3.d0) * rmt (is) ** 3
            End If
         End Do
! loop over atoms
         Do ia = 1, natoms (is)
            Do ig = 1, ngrtot
! structure factor
               t2 = - dot_product (vgc(:, ig), atposc(:, ia, is))
               zt1 = cmplx (Cos(t2), Sin(t2), 8)
! add to characteristic function in G-space
               cfunig (ig) = cfunig (ig) - zt1 * ffacg (ig)
            End Do
         End Do
      End Do
      Do ig = 1, ngrtot
         ifg = igfft (ig)
         zfft (ifg) = cfunig (ig)
      End Do
! Fourier transform to real-space
      Call zfftifc (3, ngrid, 1, zfft)
      cfunir (:) = dble (zfft(:))
! damp the characteristic function in G-space if required
      If (input%groundstate%cfdamp .Gt. 0.d0) Then
         t1 = input%groundstate%cfdamp / input%groundstate%gmaxvr
         Do ig = 1, ngrtot
            t2 = Exp (-(t1*gc(ig))**2)
            cfunig (ig) = cfunig (ig) * t2
         End Do
      End If
      Deallocate (ffacg, zfft)
      Return
End Subroutine
!EOC
