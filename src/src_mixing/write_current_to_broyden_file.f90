
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine write_current_to_broyden_file (n, iscl, potential, residual)
      Use modmixermsec, Only: record_of_last_iter, noldstepsmax, &
     & noldstepsin_file
      Implicit None
      Integer, Intent (In) :: n, iscl
      Real (8), Intent (In) :: potential (n), residual (n)
      Integer :: reclength
      Character (256), External :: outfilenamestring
      Character (256) :: filetag
      filetag = "BROYDEN"
      record_of_last_iter = Mod (record_of_last_iter, noldstepsmax) + 1
      Inquire (IoLength=reclength) potential, residual
!
      Open (23, File=outfilenamestring(filetag, 1), Access="DIRECT", &
     & Recl=reclength, Form='UNFORMATTED')
      Write (23, Rec=record_of_last_iter) potential, residual
      Close (23)
      noldstepsin_file = noldstepsin_file + 1
      noldstepsin_file = Min (noldstepsin_file, noldstepsmax)
End Subroutine
