!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phdisp
      Use modmain
      use modmpi
      Implicit None
! local variables
      Integer :: iq, i, n, iv
      Real (8) :: wmin, wmax
! allocatable arrays
      Real (8), Allocatable :: wp (:, :)
      Complex (8), Allocatable :: ev (:, :)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: dynp (:, :)
      Complex (8), Allocatable :: dynr (:, :, :)
! writeout only in master process
      if (rank .ne. 0) goto 10
! initialise universal variables
      Call init0
      Call init2
      nvp1d = size(input%phonons%phonondispplot%plot1d%path%pointarray)
      npp1d = input%phonons%phonondispplot%plot1d%path%steps
      n = 3 * natmtot
      Allocate (wp(n, npp1d))
      Allocate (ev(n, n))
      Allocate (dynq(n, n, nqpt))
      Allocate (dynp(n, n))
      Allocate (dynr(n, n, ngridq(1)*ngridq(2)*ngridq(3)))
! read in the dynamical matrices
      Call readdyn (.true.,dynq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! Fourier transform the dynamical matrices to real-space
      Call dynqtor (dynq, dynr)
! generate a set of q-point vectors along a path in the Brillouin zone
      If (allocated(dvp1d)) deallocate (dvp1d)
      Allocate (dvp1d(nvp1d))
      If (allocated(vplp1d)) deallocate (vplp1d)
      Allocate (vplp1d(3, npp1d))
      If (allocated(dpp1d)) deallocate (dpp1d)
      Allocate (dpp1d(npp1d))
      Call connect (bvec, input%phonons%phonondispplot%plot1d,nvp1d, npp1d, &
        vplp1d, dvp1d, dpp1d)
!      Call connect (bvec, nvp1d, npp1d, vvlp1d, vplp1d, dvp1d, dpp1d)
      wmin = 0.d0
      wmax = 0.d0
! compute the phonon frequencies along the path
      Do iq = 1, npp1d
! compute the dynamical matrix at this particular q-point
         Call dynrtoq (vplp1d(:, iq), dynr, dynp)
! find the phonon frequencies and eigenvectors
         Call dyndiag (dynp, wp(:, iq), ev)
         wmin = Min (wmin, wp(1, iq))
         wmax = Max (wmax, wp(n, iq))
      End Do
      wmax = wmax + (wmax-wmin) * 0.5d0
      wmin = wmin - (wmax-wmin) * 0.5d0
! output the vertex location lines
      Open (50, File='PHDLINES.OUT', Action='WRITE', Form='FORMATTED')
      Do iv = 1, nvp1d
         Write (50, '(2G18.10)') dvp1d (iv), wmin
         Write (50, '(2G18.10)') dvp1d (iv), wmax
         Write (50, '("     ")')
      End Do
      Close (50)
! output the phonon dispersion
      Open (50, File='PHDISP.OUT', Action='WRITE', Form='FORMATTED')
      Do i = 1, n
         Do iq = 1, npp1d
            Write (50, '(2G18.10)') dpp1d (iq), wp (i, iq)
         End Do
         Write (50, '("     ")')
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(phdisp):")')
      Write (*, '(" phonon dispersion written to PHDISP.OUT")')
      Write (*, '(" vertex location lines written to PHDLINES.OUT")')
      Write (*,*)
      Deallocate (wp, ev, dynq, dynp, dynr)
10    continue
#ifdef MPI
      call MPI_Barrier(MPI_Comm_World, ierr)
#endif
      Return
End Subroutine
