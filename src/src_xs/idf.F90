!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine idf
      Use modmain
      Use modinput
      Use modmpi
      Use modxs
      Use modfxcifc
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'idf'
      Integer :: iq
  ! initialise universal variables
      Call init0
      Call init1
  ! save Gamma-point variables
      Call xssave0
  ! initialize q-point set
      Call init2
      Call readfermi
  ! w-point parallelization for dielectric function
      Call genparidxran ('w', nwdf)
      Write (unitout, '("Exchange-correlation kernel type :", i4)') &
     & input%xs%tddft%fxctypenumber
      Write (unitout, '("  ", a)') trim (fxcdescr)
  ! loop over q-points
      Do iq = 1, nqpt
     ! call for q-point
         Call idfq (iq)
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): inverse d&
        &ielectric function finished for q - point:', iq
      End Do
      Call barrier
      If ((procs .Gt. 1) .And. (rank .Eq. 0)) Then
         Call idfgather
         Write (unitout, '(a)') 'Info(' // thisnam // '): inverse diele&
        &ctric function gathered for q - point:'
      End If
      Call barrier
      If (rank .Eq. 0) Then
         Do iq = 1, nqpt
        ! call for q-point
            Call xslinopt (iq)
            Write (unitout, '(a, i8)') 'Info(' // thisnam // '): &
            &TDDFT linear optics finished for q - point:', iq
         End Do
      End If
      Call barrier
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): TDDFT lin&
     &ear optics finished"
      Call genfilname (setfilext=.True.)
End Subroutine idf
