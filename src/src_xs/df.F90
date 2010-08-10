!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: df
! !INTERFACE:
!
!
Subroutine df
! !USES:
      Use modinput
      Use modmain
      Use modxs
      Use modmpi
      Use m_writegqpts
      Use m_xsgauntgen
      Use m_findgntn0
      Use m_genfilname
! !DESCRIPTION:
!   Control routine for setting up the Kohn-Sham response function or the
!   microscopic dielectric function/matrix for all specified ${\bf q}$-points.
!   Can be run with MPI parallelization for ${\bf q}$-points.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'df'
      Character (256) :: filex
      Integer :: iq
      If ( .Not. tscreen) Call genfilname (setfilext=.True.)
      Call init0
  ! initialise universal variables
      Call init1
  ! save Gamma-point variables
      Call xssave0
  ! initialize q-point set
      Call init2
      If (tscreen) Then
     ! generate Gaunt coefficients
         Call xsgauntgen (Max(input%groundstate%lmaxapw, lolmax), &
        & input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
     ! find indices for non-zero Gaunt coefficients
         Call findgntn0 (Max(input%xs%lmaxapwwf, lolmax), &
        & Max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
      End If
  ! read Fermi energy
      Call readfermi
  ! w-point parallelization for dielectric function
      If (tscreen) Then
         nwdf = 1
         Call genparidxran ('q', nqpt)
      Else
         Call genparidxran ('w', nwdf)
      End If
  ! set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
      input%xs%emattype = 1
  ! write out q-points
      Call writeqpts
  ! loop over q-points
      Do iq = qpari, qparf
         Call genfilname (iq=iq, fileext=filex)
     ! call for q-point
         Call dfq (iq)
         If (tscreen) Call writegqpts (iq, filex)
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): Kohn Sahm&
        & response function finished for q - point:', iq
      End Do
  ! synchronize
      Call barrier
      If ((procs .Gt. 1) .And. (rank .Eq. 0) .And. ( .Not. tscreen)) &
     & Call dfgather
      Call barrier
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): Kohn-Sham&
     & response function finished"
      If ( .Not. tscreen) Call genfilname (setfilext=.True.)
      If (tscreen) Call findgntn0_clear
End Subroutine df
!EOC
