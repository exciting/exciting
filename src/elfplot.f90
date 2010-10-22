!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: elfplot
! !INTERFACE:
!
!
Subroutine elfplot
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Outputs the electron localisation function (ELF) for 1D, 2D or 3D plotting.
!   The spin-averaged ELF is given by
!   $$ f_{\rm ELF}({\bf r})=\frac{1}{1+[D({\bf r})/D^0({\bf r})]^2}, $$
!   where
!   $$ D({\bf r})=\frac{1}{2}\left(\tau({\bf r})-\frac{1}{4}
!    \frac{[\nabla n({\bf r})]^2}{n({\bf r})}\right) $$
!   and
!   $$ \tau({\bf r})=\sum_{i=1}^N \left|\nabla\Psi_i({\bf r})
!    \right|^2 $$
!   is the spin-averaged kinetic energy density from the spinor wavefunctions.
!   The function $D^0$ is the kinetic energy density for the homogeneous
!   electron gas evaluated for $n({\bf r})$:
!   $$ D^0({\bf r})=\frac{3}{5}(6\pi^2)^{2/3}\left(\frac{n({\bf r})}{2}
!    \right)^{5/3}. $$
!   The ELF is useful for the topological classification of bonding. See for
!   example T. Burnus, M. A. L. Marques and E. K. U. Gross [Phys. Rev. A 71,
!   10501 (2005)].
!
! !REVISION HISTORY:
!   Created September 2003 (JKD)
!   Fixed bug found by F. Wagner (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, is, ia, ias
      Integer :: ir, i, itp, ig, ifg
      Real (8) :: t1, t2
! allocatable arrays
      Real (8), Allocatable :: rftp1 (:), rftp2 (:), rftp3 (:)
      Real (8), Allocatable :: grfmt (:, :, :)
      Real (8), Allocatable :: grfir (:)
      Real (8), Allocatable :: gwf2mt (:, :, :)
      Real (8), Allocatable :: gwf2ir (:)
      Real (8), Allocatable :: elfmt (:, :, :)
      Real (8), Allocatable :: elfir (:)
      Complex (8), Allocatable :: zfft1 (:), zfft2 (:)
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
! initialise universal variables
      Call init0
      Call init1
! allocate local arrays
      Allocate (rftp1(lmmaxvr), rftp2(lmmaxvr), rftp3(lmmaxvr))
      Allocate (grfmt(lmmaxvr, nrmtmax, 3))
      Allocate (grfir(ngrtot))
      Allocate (gwf2mt(lmmaxvr, nrmtmax, natmtot))
      Allocate (gwf2ir(ngrtot))
      Allocate (elfmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (elfir(ngrtot))
      Allocate (zfft1(ngrtot), zfft2(ngrtot))
! allocate first-variational eigenvector array
      Allocate (evecfv(nmatmax, nstfv))
! allocate second-variational eigenvector array
      Allocate (evecsv(nstsv, nstsv))
! read density and potentials from file
      Call readstate
! read Fermi energy from file
      Call readfermi
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! generate the core wavefunctions and densities
      Call gencore
! set the gradient squared to zero
      gwf2mt (:, :, :) = 0.d0
      gwf2ir (:) = 0.d0
      Do ik = 1, nkpt
! get the eigenvectors and occupancies from file
         Call getoccsv (vkl(:, ik), occsv(:, ik))
         Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(:, ik), evecsv)
! add the valence wavefunction gradient squared
         Call gwf2val (ik, evecfv, evecsv, gwf2mt, gwf2ir)
      End Do
! add core wavefunction gradient squared
      Call gwf2cr (gwf2mt)
!------------------------!
!     muffin-tin ELF     !
!------------------------!
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! compute the gradient of the density
            Call gradrfmt (input%groundstate%lmaxvr, nrmt(is), spr(:, &
           & is), lmmaxvr, nrmtmax, rhomt(:, :, ias), grfmt)
            Do ir = 1, nrmt (is)
! convert rho from spherical harmonics to spherical coordinates
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
              & lmmaxvr, rhomt(:, ir, ias), 1, 0.d0, rftp1, 1)
               rftp2 (:) = 0.d0
! compute the square of the gradient of rho
               Do i = 1, 3
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 & lmmaxvr, grfmt(:, ir, i), 1, 0.d0, rftp3, 1)
                  Do itp = 1, lmmaxvr
                     rftp2 (itp) = rftp2 (itp) + rftp3 (itp) ** 2
                  End Do
               End Do
               Do itp = 1, lmmaxvr
! D for inhomogeneous density
                  t1 = (1.d0/2.d0) * (gwf2mt(itp, ir, &
                 & ias)-(1.d0/4.d0)*rftp2(itp)/rftp1(itp))
! D0 for uniform electron gas
                  t2 = (3.d0/5.d0) * ((6.d0*pi**2)**(2.d0/3.d0)) * &
                 & (Abs(rftp1(itp))/2.d0) ** (5.d0/3.d0)
! ELF function
                  rftp3 (itp) = 1.d0 / (1.d0+(t1/t2)**2)
               End Do
! convert from spherical coordinates to spherical harmonics
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
              & lmmaxvr, rftp3, 1, 0.d0, elfmt(:, ir, ias), 1)
            End Do
         End Do
      End Do
!--------------------------!
!     interstitial ELF     !
!--------------------------!
! Fourier transform density to G-space
      zfft1 (:) = rhoir (:)
      Call zfftifc (3, ngrid,-1, zfft1)
      grfir (:) = 0.d0
      Do i = 1, 3
         zfft2 (:) = 0.d0
         Do ig = 1, ngvec
            ifg = igfft (ig)
! take the gradient
            zfft2 (ifg) = zi * vgc (i, ig) * zfft1 (ifg)
         End Do
! Fourier transform gradient to real-space
         Call zfftifc (3, ngrid, 1, zfft2)
         Do ir = 1, ngrtot
            grfir (ir) = grfir (ir) + dble (zfft2(ir)) ** 2
         End Do
      End Do
      Do ir = 1, ngrtot
! D for inhomogeneous density
         t1 = (1.d0/2.d0) * &
        & (gwf2ir(ir)-(1.d0/4.d0)*grfir(ir)/rhoir(ir))
! D0 for homogeneous electron gas
         t2 = (3.d0/5.d0) * ((6.d0*pi**2)**(2.d0/3.d0)) * &
        & (Abs(rhoir(ir))/2.d0) ** (5.d0/3.d0)
! ELF function
         elfir (ir) = 1.d0 / (1.d0+(t1/t2)**2)
      End Do
! symmetrise the ELF
      Call symrf (1, elfmt, elfir)
! plot the ELF to file
!
      If (associated(input%properties%elfplot%plot1d)) Then
!
         Call plot1d ("ELF", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & elfmt, elfir, input%properties%elfplot%plot1d)
!
         Write (*,*)
         Write (*, '("Info(elfplot):")')
         Write (*, '(" 1D ELF plot written to ELF1D.OUT")')
         Write (*, '(" vertex location lines written to ELFLINES.OUT")')
      End If
      If (associated(input%properties%elfplot%plot2d)) Then
         Call plot2d ("ELF", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & elfmt, elfir, input%properties%elfplot%plot2d)
         Write (*,*)
         Write (*, '("Info(elfplot): 2D ELF plot written to ELF2d.xml")&
        &')
      End If
      If (associated(input%properties%elfplot%plot3d)) Then
         Call plot3d ("ELF", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & elfmt, elfir, input%properties%elfplot%plot3d)
         Write (*,*)
         Write (*, '("Info(elfplot): 3D ELF plot written to ELF3d.xml")&
        &')
      End If
      Write (*,*)
      Deallocate (rftp1, rftp2, rftp3)
      Deallocate (grfmt, grfir, gwf2mt, gwf2ir, elfmt, elfir)
      Deallocate (zfft1, zfft2, evecfv, evecsv)
      Return
End Subroutine
!EOC
