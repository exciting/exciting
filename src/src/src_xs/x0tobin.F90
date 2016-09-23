!
!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine x0tobin
      Use modmain
      Use modinput
      Use modxs
      Use m_getx0
      Use m_putx0
      Use m_getunit
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'x0tobin'
      Character (256) :: filnam, filnama
      Integer :: n, iq, iw, un
      Complex (8), Allocatable :: chi0 (:, :), chi0wg (:, :, :), chi0hd &
     & (:, :)
      Logical :: tq0
      Logical, External :: tqgamma
      Call init0
  ! initialise universal variables
      Call init1
  ! save Gamma-point variables
      Call xssave0
  ! initialize q-point set
      Call init2
  ! loop over q-points
      Do iq = 1, nqpt
         tq0 = tqgamma (iq)
     ! calculate k+q and G+k+q related variables
         Call init1offs (qvkloff(1, iq))
     ! size of local field effects
         n = ngq (iq)
     ! allocate
         Allocate (chi0(n, n), chi0wg(n, 2, 3), chi0hd(3, 3))
     ! filenames
         Call genfilname (asc=.True., basename='X0', bzsampl=bzsampl, &
        & acont=input%xs%tddft%acont, nar= .Not. input%xs%tddft%aresdf, &
        & tord=input%xs%tddft%torddf, markfxcbse=tfxcbse, iqmt=iq, &
        & filnam=filnama)
         Call genfilname (basename='X0', bzsampl=bzsampl, &
        & acont=input%xs%tddft%acont, nar= .Not. input%xs%tddft%aresdf, &
        & tord=input%xs%tddft%torddf, markfxcbse=tfxcbse, iqmt=iq, &
        & filnam=filnam)
     ! open file to read ASCI
         Call getunit (un)
         Open (Unit=un, File=trim(filnama), Form='formatted', Action='r&
        &ead', Status='old')
         Do iw = 1, nwdf
        ! read from ASCII file
            If (tq0) Then
               Read (un,*) ngq (iq), vql (:, iq), chi0, chi0wg, chi0hd
            Else
               Read (un,*) ngq (iq), vql (:, iq), chi0
            End If
        ! write to binary file
            Call putx0 (tq0, iq, iw, trim(filnam), '', chi0, chi0wg, &
           & chi0hd)
         End Do
     ! close file
         Close (un)
         Deallocate (chi0, chi0wg, chi0hd)
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): Kohn Sham&
        & response function converted to binary file for q - point:', &
        & iq
      End Do
      Call genfilname (setfilext=.True.)
End Subroutine x0tobin
