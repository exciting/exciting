
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine write_current_to_broyden_file (n, iscl, potential, residual)
      Use modmixermsec, Only: record_of_last_iter, noldstepsmax, &
     & noldstepsin_file
      use mod_large_io, only: inquire_large, open_direct_unformatted_large
      use mod_misc, only: scrpath
      use precision, only: i32, dp, long_int, str_256

      Implicit None
      Integer(long_int), Intent (In) :: n
      Integer(i32), Intent (In) :: iscl
      Real (dp), Intent (In) :: potential (n), residual (n)
      Integer(long_int) :: reclength
      integer(i32) :: io_unit
      Character (str_256), External :: outfilenamestring
      Character (str_256) :: filetag
      filetag = "BROYDEN"
      record_of_last_iter = Mod (record_of_last_iter, noldstepsmax) + 1
      
      call inquire_large( reclength, potential, residual )
      call open_direct_unformatted_large( io_unit, trim( scrpath ) // "BROYDEN.OUT", "write", reclength, "unknown" )
      write( io_unit, rec = record_of_last_iter ) potential, residual
      close( io_unit )

      noldstepsin_file = noldstepsin_file + 1
      noldstepsin_file = Min (noldstepsin_file, noldstepsmax)
End Subroutine
