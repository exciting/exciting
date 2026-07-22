
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine readbroydsteps_and_init_SY (noldsteps, n, S, Y, potential, residual)
!
      use modmixermsec, Only: record_of_last_iter, noldstepsmax
      use modmain
      use mod_large_io, only: inquire_large, open_direct_unformatted_large
      use precision, only: i32, dp, long_int, str_256

      Implicit None
      integer(i32), intent (inout)     :: noldsteps
      integer(long_int), intent(inout) :: n
      real(dp), intent (out) :: S (n, noldstepsmax), Y (n, noldstepsmax)
      real(dp), intent (in) :: potential(n), residual(n)
      integer(i32) :: i, skipp, io_unit, rectoread, firstrec
      character(str_256), External :: outfilenamestring
      integer(long_int) :: reclength

      call inquire_large( reclength, potential, residual )
      call open_direct_unformatted_large( io_unit, trim( scrpath ) // "BROYDEN.OUT", "read", reclength, "old" )

      firstrec = 1
      if ( noldsteps >= noldstepsmax ) firstrec = mod( record_of_last_iter, noldstepsmax ) + 1

      S = 0
      Y = 0
      skipp = noldstepsmax - noldsteps
      Do i = 1, noldsteps
         rectoread = firstrec - 1 + i
         If (rectoread .Gt. noldstepsmax) rectoread = rectoread - &
        & noldstepsmax
         Read (io_unit, Rec=rectoread) S (:, i+skipp), Y (:, i+skipp)
      End Do
      Close (io_unit)
!
      Do i = 1, noldsteps
         S (:, i+skipp) = S (:, i+skipp) - potential
         Y (:, i+skipp) = Y (:, i+skipp) - residual
      End Do
!
!
End Subroutine
