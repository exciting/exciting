
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Subroutine rescaleYS (noldsteps, n, S, Y, potential, residual)
      Use modmixermsec, Only: noldstepsmax, PWHIST, FHIST, CLMHIST, &
     & yhist, FHIST, icond, scl_plane, MSECINFO
      Use mod_Gvector, Only: ngrtot
      Use mod_muffin_tin, Only: lmmaxvr, nrmtmax
      Use mod_atoms, Only: natmtot
      Implicit None
      Integer, Intent (In) :: n, noldsteps
      Real (8), Intent (Inout) :: S (n, noldstepsmax), Y (n, &
     & noldstepsmax)
      Real (8), Intent (Inout) :: potential (n), residual (n)
      Real (8) :: PWAVE, CLAVE, Rescale, T1
      Integer :: i, j, k, nmt, firstpw, lastpw
!
      nmt = lmmaxvr * nrmtmax * natmtot
      firstpw = n - ngrtot
      lastpw = n
      FHIST (noldsteps) = dot_product (residual, residual)
!
      PWHIST (noldsteps) = dot_product (residual(firstpw:lastpw), &
     & residual(firstpw:lastpw))
      CLMHIST (noldsteps) = FHIST (noldsteps) - PWHIST (noldsteps)
!Preconditioner Omega_n Pg 21
#ifdef Full_Understanding_how_this_should_work
      PWAVE = 0.
      CLAVE = 0.
      Do i = 1, noldsteps
         PWAVE = PWAVE + Sqrt (PWHIST(i)/FHIST(i))
         CLAVE = CLAVE + Sqrt (CLMHIST(i)/FHIST(i))
      End Do
!
!       For te PW's rescale so residue is comparable to that of CLMs
!       Default is icond=1, take sqrt
!       This makes sense because relative weights appear as Rescale**2 in algorithm
      Rescale = CLAVE / PWAVE
      If (icond .Gt. 0) Rescale = Sqrt (Rescale)
        !
        !rescale doesnt work yet in exciting
      Rescale = 1
        !
      MSECINFO (1) = Rescale
1002  Format (':INFO : ', a, 10D11.3)
      Write (*, 1002) ' Dynamic rescale ', Rescale
      potential (firstpw:lastpw) = potential (firstpw:lastpw) * Rescale
      residual (firstpw:lastpw) = residual (firstpw:lastpw) * Rescale
      Y (firstpw:lastpw, 1:noldstepsmax) = Y (firstpw:lastpw, &
     & 1:noldstepsmax) * Rescale
      S (firstpw:lastpw, 1:noldstepsmax) = S (firstpw:lastpw, &
     & 1:noldstepsmax) * Rescale
      PWHIST (1:noldsteps) = PWHIST (1:noldsteps) * Rescale * Rescale
      FHIST (1:noldsteps) = PWHIST (1:noldsteps) + CLMHIST &
     & (1:noldsteps)
      scl_plane = scl_plane * Rescale
!
!
#endif
!
      Do j = 1, noldstepsmax
         T1 = 0.
!
         Do k = 1, n
            T1 = T1 + Y (k, j) * Y (k, j)
         End Do
         yhist (j) = T1
!       Renormalize
         T1 = 1.D0 / Sqrt (T1)
         Do k = 1, n
            S (k, j) = S (k, j) * T1
            Y (k, j) = Y (k, j) * T1
         End Do
      End Do
End Subroutine
