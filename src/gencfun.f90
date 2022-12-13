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
      use constants, only: fourpi, zzero
      use sirius_api, only: get_step_function
      use sirius_init, only: sirius_options
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
      integer :: ig
      real (8) :: t1, t2
      integer :: ngv
! allocatable arrays
      complex (8), allocatable :: zfft (:)

      ! Exciting generates the full set of G-vectors equal to the size of FFT box in real space.
      ! The full set of G-vectors is used only here, in generating step function. For SIRIUS
      ! only G-vectors within a plane-wave cutoff are generated. This is not problem as it leads to
      ! negligible errors in total energy.
      ngv = ngrtot
      if (associated(input%groundstate%sirius)) then
        ngv = ngvec
      end if

! allocate global characteristic function arrays
      If (allocated(cfunig)) deallocate (cfunig)
      Allocate (cfunig(ngv))
      If (allocated(cfunir)) deallocate (cfunir)
      Allocate (cfunir(ngrtot))

      if (associated(input%groundstate%sirius) .and. sirius_options%use_c_function) then
        call get_step_function( cfunig, cfunir, ngrtot )
        return
      endif

      Allocate (zfft(ngrtot), source=zzero)
      call gencfunig( ngv, gc, vgc, cfunig)

      zfft( igfft(1:ngv)) = cfunig

! Fourier transform to real-space
      Call zfftifc (3, ngrid, 1, zfft)
      cfunir (:) = dble (zfft(:))
! damp the characteristic function in G-space if required
      If (input%groundstate%cfdamp .Gt. 0.d0) Then
         t1 = input%groundstate%cfdamp / input%groundstate%gmaxvr
         Do ig = 1, ngv
            t2 = Exp (-(t1*gc(ig))**2)
            cfunig (ig) = cfunig (ig) * t2
         End Do
      End If
      Deallocate (zfft)
      Return
End Subroutine
!EOC

!> Generate the Fourier coefficients of the smooth step function
!> for a given set of G-vectors given by
!> \[ \hat{\Theta}({\bf G}) = \begin{cases}
!>    \delta_{{bf G}0} - \sum\alpha {\rm e}^{-{\rm i}{\bf G}\cdot{\bf \tau}_\alpha} 
!>    \frac{4\pi R_\alpha^3}{3\Omega} & {\bf G} = 0 \\
!>    \sum\alpha {\rm e}^{-{\rm i}{\bf G}\cdot{\bf \tau}_\alpha} 
!>    \frac{4\pi R_\alpha^3}{\Omega} \frac{j_1(G\,R_\alpha)}{G\,R_\alpha} & {\bf G} \neq 0
!>    \end{cases} \;. \]
subroutine gencfunig( ng, gc, vgc, cfunig)
  use precision, only: dp
  use modinput
  use constants, only: zzero, zone, fourpi
  use mod_lattice, only: omega
  use mod_atoms, only: nspecies, natoms, atposc
  use mod_muffin_tin, only: rmt
  
  !Note: assumed-size declaration needed because routine is not in a module.
  !> number of G-vectors
  integer, intent(in) :: ng
  !> length of the G-vectors
  real(dp), intent(in) :: gc(*)
  !> G-vectors in Cartesian coordinates
  real(dp), intent(in) :: vgc(3,*)
  !> Fourier coefficients of smooth step function
  complex(dp), intent(out) :: cfunig(*)

  integer :: is, ia, i, ig
  real(dp) :: t1, t2

  integer, allocatable :: finite_g_indices(:), zero_g_indices(:)
  real(dp), allocatable :: ffacg(:)

  t1 = fourpi/omega

  allocate( ffacg(ng))
  finite_g_indices = pack( [(ig,ig=1,ng)], [(gc(ig) > input%structure%epslat, ig=1, ng)])
  zero_g_indices = pack( [(ig,ig=1,ng)], [(gc(ig) <= input%structure%epslat, ig=1, ng)])
  
  cfunig(finite_g_indices) = zzero
  cfunig(zero_g_indices) = zone

!$omp parallel default(shared) private(is,ia,i,ig,t2)
  do is = 1, nspecies
    ! smooth step function form factors for each species
!$omp do
    do i = 1, size( finite_g_indices)
      ig = finite_g_indices(i)
      t2 = gc(ig) * rmt(is)
      ffacg(ig) = t1 * (sin(t2) - t2 * cos(t2)) / (gc(ig)**3)
    end do
!$omp end do
!$omp do
    do i = 1, size( zero_g_indices)
      ig = zero_g_indices(i)
      ffacg(ig) = t1 / 3._dp * rmt(is)**3
    end do
!$omp end do
!$omp barrier

    do ia = 1, natoms(is)
!$omp do
      do ig = 1, ng
        t2 = -dot_product( vgc(:,ig), atposc(:,ia,is))
        ! add to characteristic function in G-space
        cfunig(ig) = cfunig(ig) - ffacg(ig)*cmplx(cos(t2),sin(t2),dp)
      end do
!$omp end do
    end do

  end do
!$omp end parallel
end subroutine
