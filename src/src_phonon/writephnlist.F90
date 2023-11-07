
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine writephnlist(nppt,vpl,twrev,fname)
      Use modmain
      use modmpi
      use phonons_util, only: ph_util_setup_interpolation, ph_util_interpolate, ph_util_diag_dynmat
      use phonons_io_util, only: ph_io_read_dielten, ph_io_read_borncharge
      use matrix_fourier_interpolation, only: mfi_type
      Implicit None
! arguments
      integer, intent(in) :: nppt
      real(8), intent(in) :: vpl(3,nppt)
      logical, intent(in) :: twrev
      character(*), intent(in) :: fname
! local variables
      Integer :: n, iq, i, j, is, ia, ip
      Real (8) :: er, ei, eeps, dielten(3, 3)
      type(mfi_type) :: mfi
      logical :: success
! allocatable arrays
      Real (8), Allocatable :: borncharge (:, :, :)
      Real (8), Allocatable :: w (:)
      Complex (8), Allocatable :: ev (:, :)
      Complex (8), Allocatable :: dynp (:, :, :)
      Complex (8), Allocatable :: dynr (:, :, :)
      n = 3 * natmtot
      Allocate (w(n))
      Allocate (ev(n, n))
      ! set threshold for the eigenvector components
      eeps = 1E-15
      ! try to read dielectric tensor and Born effective charges
      if( input%phonons%polar ) then
        call ph_io_read_dielten( dielten, 'EPSINF.OUT', success )
        call terminate_if_false( success, '(writephnlist) &
          Failed to read dielectric tensor from EPSINF.OUT.' )
        call ph_io_read_borncharge( borncharge, 'ZSTAR.OUT', success )
        call terminate_if_false( success, '(writephnlist) &
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
      ! interpolate dynamical matrices on list
      if( input%phonons%polar ) then
        call ph_util_interpolate( nppt, vpl(:, 1:nppt), mfi, dynr, dynp, &
          dielten=dielten, borncharge=borncharge )
      else
        call ph_util_interpolate( nppt, vpl(:, 1:nppt), mfi, dynr, dynp )
      end if
      ! write phonons
      Open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
      Do iq = 1, nppt
         call ph_util_diag_dynmat( dynp(:, :, iq), w, ev )
         Write (50,*)
         Write (50, '(I6, 3G18.10, " : q-point, vpl")') iq, vpl &
        & (:, iq)
         Do j = 1, n
            if (twrev) Write (50,*)
            Write (50, '(I6, G18.10, " : mode, frequency")') j, w (j)
            if (twrev) then
                i = 0
                Do is = 1, nspecies
                   Do ia = 1, natoms (is)
                      Do ip = 1, 3
                         i = i + 1
                         er = dreal(ev(i, j))
                         if (dabs(er) .lt. eeps) then 
                             er = 0.E0
                         end if
                         ei = dimag(ev(i, j))
                         if (dabs(ei) .lt. eeps) then 
                             ei = 0.E0
                         end if
                         ev(i, j) = cmplx(er, ei)
                         If (i .Eq. 1) Then
                            Write (50, '(3I4, 2G18.10, " : species, atom, polarisation, eigenvector")')&
                           & is, ia, ip, ev (i, j)
                         Else
                            Write (50, '(3I4, 2G18.10)') is, ia, ip, ev (i, j)
                         End If
                      End Do
                   End Do
                End Do
            end if
         End Do
         Write (50,*)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writephn): phonon frequencies (and eigenvectors) &
     &written to ",a)') trim(fname)
      Write (*, '(" for all q-vectors in the phwrite list")')
      Write (*,*)
      ! clean up
      Deallocate (w, ev )
      if( allocated( dynr ) ) deallocate( dynr )
      if( allocated( dynp ) ) deallocate( dynp )
      if( allocated( borncharge ) ) deallocate( borncharge )
      call mfi%destroy
      Return
End Subroutine
