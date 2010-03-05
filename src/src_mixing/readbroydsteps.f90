
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine readbroydsteps_and_init_SY (noldsteps, n, S, Y, potential, &
& residual)
!
      Use modmixermsec, Only: record_of_last_iter, noldstepsmax
      Implicit None
      Integer, Intent (In) :: noldsteps, n
      Real (8), Intent (Out) :: S (n, noldstepsmax), Y (n, &
     & noldstepsmax)
      Real (8), Intent (In) :: potential (n), residual (n)
      Integer :: i, skipp
      Character (256), External :: outfilenamestring
      Character (256) :: filetag
      Integer :: reclength, rectoread, firstrec
      Inquire (IoLength=reclength) potential, residual
      filetag = "BROYDEN"
      Open (23, File=outfilenamestring(filetag, 1), Access="DIRECT", &
     & Recl=reclength, Action="READ", Form='UNFORMATTED')
      If (noldsteps .Lt. noldstepsmax) Then
         firstrec = 1
      Else
         firstrec = Mod (record_of_last_iter, noldstepsmax) + 1
      End If
      S = 0
      Y = 0
      skipp = noldstepsmax - noldsteps
      Do i = 1, noldsteps
         rectoread = firstrec - 1 + i
         If (rectoread .Gt. noldstepsmax) rectoread = rectoread - &
        & noldstepsmax
         Read (23, Rec=rectoread) S (:, i+skipp), Y (:, i+skipp)
      End Do
      Close (23)
!
      Do i = 1, noldsteps
         S (:, i+skipp) = S (:, i+skipp) - potential
         Y (:, i+skipp) = Y (:, i+skipp) - residual
      End Do
!
!
End Subroutine
