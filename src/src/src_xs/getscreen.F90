!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getscreen (iqr, ngq, scrh, scrw, scrb)
      Use modxs, Only:
      Use modinput
      Use m_genfilname
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iqr, ngq
      Complex (8), Intent (Out) :: scrb (ngq, ngq), scrw (ngq, 2, 3), &
     & scrh (3, 3)
  ! local variables
      Character (256) :: fname
      Real (8) :: rm (3, 3, 3)
      Integer :: igq1, igq2, i, j, it1, it2, un, bzsampl
  ! sampling of Brillouin zone
      bzsampl = 0
      If (input%xs%tetra%tetradf) bzsampl = 1
  ! read in screening
      Call genfilname (basename='SCREEN', iq=iqr, bzsampl=bzsampl, &
     & filnam=fname)
      Call getunit (un)
      Open (un, File=trim(fname), Form='formatted', Action='read', &
     & Status='old')
      Do igq1 = 1, ngq
         Do igq2 = 1, ngq
            If (iqr .Eq. 1) Then
               If ((igq1 .Eq. 1) .And. (igq2 .Eq. 1)) Then
                  Read (un,*) ((it1, it2, rm(1, i, j), rm(2, i, j), &
                 & rm(3, i, j), j=1, 3), i=1, 3)
                  scrh (:, :) = cmplx (rm(1, :, :), rm(2, :, :), 8)
               End If
               If ((igq1 .Eq. 1) .And. (igq2 .Ne. 1)) Then
                  Read (un,*) (it1, it2, rm(1, 1, j), rm(2, 1, j), &
                 & rm(3, 1, j), j=1, 3)
                  scrw (igq2, 1, :) = cmplx (rm(1, 1, :), rm(2, 1, :), &
                 & 8)
               End If
               If ((igq1 .Ne. 1) .And. (igq2 .Eq. 1)) Then
                  Read (un,*) (it1, it2, rm(1, 1, j), rm(2, 1, j), &
                 & rm(3, 1, j), j=1, 3)
                  scrw (igq1, 2, :) = cmplx (rm(1, 1, :), rm(2, 1, :), &
                 & 8)
               End If
               If ((igq1 .Ne. 1) .And. (igq2 .Ne. 1)) read (un,*) it1, &
              & it2, rm (1, 1, 1), rm (2, 1, 1), rm (3, 1, 1)
               scrb (igq1, igq2) = cmplx (rm(1, 1, 1), rm(2, 1, 1), 8)
            Else
               Read (un,*) it1, it2, rm (1, 1, 1), rm (2, 1, 1), rm (3, &
              & 1, 1)
               scrb (igq1, igq2) = cmplx (rm(1, 1, 1), rm(2, 1, 1), 8)
            End If
         End Do
      End Do
      Close (un)
End Subroutine getscreen
