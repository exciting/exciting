!
!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine x0toasc
      Use modmain
      Use modinput
      Use modxs
      Use m_getx0
      Use m_putx0
      Use m_getunit
      Use m_genfilname
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'x0toasc'
      Character (256) :: filnam, filnama
      Integer :: n, iq, igq, igqp, iw, oct1, oct2, noct, un
      Complex (8) :: zt
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
     ! open file to write ASCI
         Call getunit (un)
         Open (Unit=un, File=trim(filnama), Form='formatted', Action='w&
        &rite', Status='replace')
         noct = 1
         If (tq0) noct = 3
         Do iw = 1, nwdf
        ! read from binary file
            Call getx0 (tq0, iq, iw, trim(filnam), '', chi0, chi0wg, &
           & chi0hd)
        ! write to ASCII file
            Do igq = 1, ngq (iq)
               Do igqp = 1, ngq (iq)
                  If (tq0) Then
                     Do oct1 = 1, noct
                        Do oct2 = 1, noct
                       ! head
                           If ((igq .Eq. 1) .And. (igqp .Eq. 1)) chi0 &
                          & (igq, igqp) = chi0hd (oct1, oct2)
                       ! wings
                           If ((n .Gt. 1) .And. (igq .Eq. 1) .And. &
                          & (igqp .Gt. 1)) chi0 (igq, igqp) = chi0wg &
                          & (igqp, 1, oct1)
                           If ((n .Gt. 1) .And. (igq .Gt. 1) .And. &
                          & (igqp .Eq. 1)) chi0 (igq, igqp) = chi0wg &
                          & (igq, 2, oct2)
                           zt = chi0 (igq, igqp)
                           Write (un, '(6i6, 3g18.10)') iq, iw, igq, &
                          & igqp, oct1, oct2, zt, Abs (zt)
                        End Do
                     End Do
                  Else
                     zt = chi0 (igq, igqp)
                     Write (un, '(6i6, 3g18.10)') iq, iw, igq, igqp, 0, &
                    & 0, zt, Abs (zt)
                  End If
               End Do
            End Do
         End Do
     ! close file
         Close (un)
         Deallocate (chi0, chi0wg, chi0hd)
         Write (unitout, '(a, i8)') 'Info(' // thisnam // '): Kohn Sham&
        & response function converted to ASCII file for q - point:', iq
      End Do
      Call genfilname (setfilext=.True.)
End Subroutine x0toasc
