!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writeematasc
      Use modmain
      Use modinput
      Use modxs
      Use m_getunit
      Use m_getemat
      Use m_genfilname
      Implicit None
      Character (256) :: filnam
      Integer :: un, iq, ik, i, j, ib, jb, igq
      Complex (8) :: zt
      Call init0
      Call init1
      Call xssave0
      Call init2
      Call readfermi
      Call getunit (un)
  ! loop over q-points
      Do iq = 1, nqpt
         Call genfilname (iqmt=iq, setfilext=.True.)
     ! calculate k+q and G+k+q related variables
         Call init1offs (qvkloff(1, iq))
     ! find highest (partially) occupied and lowest (partially) unoccupied
     ! states
        call findocclims(iq, ikmapikq(:,iq), istocc0, istunocc0, isto0, isto, istu0, istu)
        istunocc = istunocc0
        istocc = istocc0
     ! set limits for band combinations
         Call ematbdcmbs (input%xs%emattype)
         If (allocated(xiou)) deallocate (xiou)
         Allocate (xiou(nst1, nst2, ngq(iq)))
         If (input%xs%emattype .Ne. 0) Then
            If (allocated(xiuo)) deallocate (xiuo)
            Allocate (xiuo(nst3, nst4, ngq(iq)))
         End If
     ! filename for matrix elements file
         Call genfilname (basename='EMAT', asc=.True., iqmt=iq, &
        & etype=input%xs%emattype, filnam=filnam)
         Open (un, File=trim(filnam), Action='write')
     ! read matrix elements of exponential expression
         Call genfilname (basename='EMAT', iqmt=iq, &
        & etype=input%xs%emattype, filnam=fnemat)
     ! loop over k-points
         Do ik = 1, nkpt
            If (input%xs%emattype .Eq. 0) Then
               Call getemat (iq, ik, .True., trim(fnemat), ngq(iq), &
              & istl1, istu1, istl2, istu2, xiou)
            Else
               Call getemat (iq, ik, .True., trim(fnemat), ngq(iq), &
              & istl1, istu1, istl2, istu2, xiou, istl3, istu3, istl4, &
              & istu4, xiuo)
            End If
            Do igq = 1, ngq (iq)
               Do i = 1, nst1
                  ib = i + istl1 - 1
                  Do j = 1, nst2
                     jb = j + istl2 - 1
                     zt = xiou (i, j, igq)
                     Write (un, '(5i8, 3g18.10)') iq, ik, igq, ib, jb, &
                    & zt, Abs (zt) ** 2
                  End Do
               End Do
            End Do
            Do igq = 1, ngq (iq)
               Do i = 1, nst3
                  ib = i + istl3 - 1
                  Do j = 1, nst4
                     jb = j + istl4 - 1
                     zt = xiuo (i, j, igq)
                     Write (un, '(5i8, 3g18.10)') iq, ik, igq, ib, jb, &
                    & zt, Abs (zt) ** 2
                  End Do
               End Do
            End Do
         End Do ! ik
         Close (un)
         Deallocate (xiou)
         If (input%xs%emattype .Ne. 0) deallocate (xiuo)
     ! end loop over q-points
      End Do
      Call genfilname (setfilext=.True.)
End Subroutine writeematasc
