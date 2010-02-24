!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phlwidth
      Use modmain
      Implicit None
! local variables
      Integer :: n, i, j, iq, iv
      Integer :: lwork, info
      Real (8) :: gmin, gmax
! allocatable arrays
      Real (8), Allocatable :: wq (:)
      Real (8), Allocatable :: gq (:, :)
      Real (8), Allocatable :: gp (:, :)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: ev (:, :), b (:, :)
      Complex (8), Allocatable :: gmq (:, :, :)
      Complex (8), Allocatable :: gmr (:, :, :)
      Complex (8), Allocatable :: gmp (:, :)
      Complex (8), Allocatable :: work (:)
! initialise universal variables
      Call init0
      Call init2
      n = 3 * natmtot
      Allocate (wq(n))
      Allocate (gq(n, nqpt))
      Allocate (gp(n, npp1d))
      Allocate (rwork(3*n))
      Allocate (dynq(n, n, nqpt))
      Allocate (ev(n, n), b(n, n))
      Allocate (gmq(n, n, nqpt))
      Allocate (gmr(n, n, ngridq(1)*ngridq(2)*ngridq(3)))
      Allocate (gmp(n, n))
      lwork = 2 * n
      Allocate (work(lwork))
! read in the dynamical matrices
      Call readdyn (.true.,dynq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! read in the phonon linewidths for each q-point
      Call readgamma (gq)
! loop over phonon q-points
      Do iq = 1, nqpt
! diagonalise the dynamical matrix
         Call dyndiag (dynq(:, :, iq), wq, ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues are the phonon linewidths
         Do i = 1, n
            Do j = 1, n
               b (i, j) = gq (i, iq) * conjg (ev(j, i))
            End Do
         End Do
         Call zgemm ('N', 'N', n, n, n, zone, ev, n, b, n, zzero, &
        & gmq(:, :, iq), n)
      End Do
! Fourier transform the gamma matrices to real-space
      Call dynqtor (gmq, gmr)
! generate a set of q-point vectors along a path in the Brillouin zone
      Call connect (bvec, nvp1d, npp1d, vvlp1d, vplp1d, dvp1d, dpp1d)
      gmin = 1.d8
      gmax = 0.d0
! compute the linewidths along the path
      Do iq = 1, npp1d
! compute the gamma matrix at this particular q-point
         Call dynrtoq (vplp1d(:, iq), gmr, gmp)
! diagonalise the gamma matrix
         Call zheev ('N', 'U', n, gmp, n, gp(:, iq), work, lwork, &
        & rwork, info)
         gmin = Min (gmin, gp(1, iq))
         gmax = Max (gmax, gp(n, iq))
      End Do
      gmax = gmax + (gmax-gmin) * 0.5d0
      gmin = gmin - (gmax-gmin) * 0.5d0
! output the vertex location lines
      Open (50, File='PHLWLINES.OUT', Action='WRITE', Form='FORMATTED')
      Do iv = 1, nvp1d
         Write (50, '(2G18.10)') dvp1d (iv), gmin
         Write (50, '(2G18.10)') dvp1d (iv), gmax
         Write (50, '("     ")')
      End Do
      Close (50)
! output the phonon linewidth dispersion
      Open (50, File='PHLWIDTH.OUT', Action='WRITE', Form='FORMATTED')
      Do i = 1, n
         Do iq = 1, npp1d
            Write (50, '(2G18.10)') dpp1d (iq), gp (i, iq)
         End Do
         Write (50, '("     ")')
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(phlwidth):")')
      Write (*, '(" phonon linewidth dispersion written to PHLWIDTH.OUT&
     &")')
      Write (*, '(" vertex location lines written to PHLWLINES.OUT")')
      Write (*,*)
      Deallocate (wq, gq, gp, rwork, dynq)
      Deallocate (ev, b, gmq, gmr, gmp, work)
      Return
End Subroutine
