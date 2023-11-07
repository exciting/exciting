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
      use phonons_util, only: ph_util_setup_interpolation, ph_util_interpolate, ph_util_diag_dynmat
      use phonons_io_util, only: ph_io_read_dielten, ph_io_read_borncharge
      use bz_path
      use matrix_fourier_interpolation, only: mfi_type
      Implicit None
! local variables
      Integer :: iq, i, n, iv
      Real (8) :: wmin, wmax, dielten(3, 3)
      type(mfi_type) :: mfi
      type(bz_path_type) :: path
      logical :: success
! allocatable arrays
      Real (8), Allocatable :: borncharge (:, :, :)
      Real (8), Allocatable :: wp (:, :)
      Complex (8), Allocatable :: ev (:, :)
      Complex (8), Allocatable :: dynp (:, :, :)
      Complex (8), Allocatable :: dynr (:, :, :)
! writeout only in master process
      if (rank .ne. 0) goto 10
! initialise universal variables
      Call init0
      Call init2
      n = 3 * natmtot
      ! try to read dielectric tensor and Born effective charges
      if( input%phonons%polar ) then
        call ph_io_read_dielten( dielten, 'EPSINF.OUT', success )
        call terminate_if_false( success, '(phdisp) &
          Failed to read dielectric tensor from EPSINF.OUT.' )
        call ph_io_read_borncharge( borncharge, 'ZSTAR.OUT', success )
        call terminate_if_false( success, '(phdisp) &
          Failed to read Born effective charge tensors from ZSTAR.OUT.' )
      end if
      ! set up interpolation
      if( input%phonons%polar ) then
        call ph_util_setup_interpolation( bvec, ngridq, nqpt, ivq(:, 1:nqpt), vql(:, 1:nqpt), mfi, dynr, &
          dielten=dielten, borncharge=borncharge, &
          sumrule=input%phonons%sumrule )
      else
        call ph_util_setup_interpolation( bvec, ngridq, nqpt, ivq(:, 1:nqpt), vql(:, 1:nqpt), mfi, dynr, &
          sumrule=input%phonons%sumrule )
      end if
      ! generate Brillouine zone path
      if( input%phonons%polar ) then
        path = bz_path_type( bvec, input%phonons%phonondispplot%plot1d%path, gamma_offset=1.d-6 )
      else
        path = bz_path_type( bvec, input%phonons%phonondispplot%plot1d%path )
      end if
      Allocate (wp(n, path%num_points))
      Allocate (ev(n, n))
      wmin = 0.d0
      wmax = 0.d0
      ! interpolate dynamical matrix along path
      if( input%phonons%polar ) then
        call ph_util_interpolate( path%num_points, &
          reshape( [(path%points(iq)%coord_lat, iq=1, path%num_points)], [3, path%num_points] ), &
          mfi, dynr, dynp, &
          dielten=dielten, borncharge=borncharge )
      else
        call ph_util_interpolate( path%num_points, &
          reshape( [(path%points(iq)%coord_lat, iq=1, path%num_points)], [3, path%num_points] ), &
          mfi, dynr, dynp )
      end if
      ! compute the phonon frequencies along the path
      Do iq = 1, path%num_points
         ! find the phonon frequencies and eigenvectors
         call ph_util_diag_dynmat( dynp(:, :, iq), wp(:, iq), ev )
         wmin = Min (wmin, wp(1, iq))
         wmax = Max (wmax, wp(n, iq))
      End Do
      wmax = wmax + (wmax-wmin) * 0.5d0
      wmin = wmin - (wmax-wmin) * 0.5d0
! output the vertex location lines
      Open (50, File='PHDLINES.OUT', Action='WRITE', Form='FORMATTED')
      Do iv = 1, path%num_vertices
         Write (50, '(2G18.10)') path%vertices(iv)%distance, wmin
         Write (50, '(2G18.10)') path%vertices(iv)%distance, wmax
         Write (50, '("     ")')
      End Do
      Close (50)
! output the phonon dispersion
      Open (50, File='PHDISP.OUT', Action='WRITE', Form='FORMATTED')
      Do i = 1, n
         Do iq = 1, path%num_points
            Write (50, '(2G18.10)') path%points(iq)%distance, wp (i, iq)
         End Do
         Write (50, '("     ")')
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(phdisp):")')
      Write (*, '(" phonon dispersion written to PHDISP.OUT")')
      Write (*, '(" vertex location lines written to PHDLINES.OUT")')
      Write (*,*)
      Deallocate (wp, ev, dynp, dynr)
      call mfi%destroy
10    continue
#ifdef MPI
      call MPI_Barrier(MPI_Comm_World, ierr)
#endif
      Return
End Subroutine
