!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writeemat
      Use modmain
      Use modinput
      Use modxs
      Use modmpi
      Use m_xsgauntgen
      Use m_findgntn0
      Use m_filedel
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'writeemat'
      Integer :: iq
  ! initialise universal variables
      Call init0
      Call init1
      Call init2
  ! k-point parallelization for TDDFT
      If ((task .Ge. 300) .And. (task .Le. 399)) Call genparidxran ('k',&
     &  nkpt)
  ! q-point parallelization for screening
      If ((task .Ge. 400) .And. (task .Le. 499)) Call genparidxran ('q',&
     &  nqpt)
   ! write q-point set
      If (rank .Eq. 0) Call writeqpts
  ! read Fermi energy from file
      Call readfermi
  ! save variables for the Gamma q-point
      Call xssave0
  ! generate Gaunt coefficients
      Call xsgauntgen (Max(input%groundstate%lmaxapw, lolmax), &
     & input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
  ! find indices for non-zero Gaunt coefficients
      Call findgntn0 (Max(input%xs%lmaxapwwf, lolmax), &
     & Max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
      Write (unitout, '(a, 3i5)') 'Info(' // thisnam // '): Gaunt coeff&
     &icients generated within lmax values:', &
     & input%groundstate%lmaxapw, input%xs%lmaxemat, &
     & input%groundstate%lmaxapw
      Write (unitout, '(a, i5)') 'Info(' // thisnam // '): number of q-&
     &points: ', nqpt
      Call flushifc (unitout)
  ! loop over q-points
      Do iq = 1, nqpt
     ! call for q-point
         Call ematq (iq)
         Write (unitout, '(a, i5)') 'Info(' // thisnam // '): matrix el&
        &ements of the exponentials finished for q - point:', iq
         Call flushifc (unitout)
      End Do
  ! synchronize
      Call barrier
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): matrix el&
     &ements of exponential expression finished"
      Call findgntn0_clear
      Call genfilname (setfilext=.True.)
End Subroutine writeemat
