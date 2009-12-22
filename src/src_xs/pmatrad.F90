!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine pmatrad
      Use modmain
      Use modinput
      Use modxs
      Implicit None
  ! local variables
      Integer :: is, ia, ias, nr, ir
      Integer :: l1, m1, lm1, l3, m3, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, j
  ! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
  ! allocatable arrays
      Real (8), Allocatable :: fapw (:), flo (:), dapwfr (:, :, :, :, &
     & :), dlofr (:, :, :, :, :)
  ! allocate local arrays for radial derivatives
      Allocate (fapw(nrmtmax))
      Allocate (dapwfr(lmmaxapw, nrmtmax, 3, apwordmax, lmmaxapw))
      dapwfr (:, :, :, :, :) = 0.d0
      If (nlotot .Gt. 0) Then
         Allocate (flo(nrmtmax))
         Allocate (dlofr(lmmaxapw, nrmtmax, 3, nlomax,-lolmax:lolmax))
         dlofr (:, :, :, :, :) = 0.d0
      End If
      ripaa (:, :, :, :, :, :) = 0.d0
      If (nlotot .Gt. 0) Then
         ripalo (:, :, :, :, :, :) = 0.d0
         riploa (:, :, :, :, :, :) = 0.d0
         riplolo (:, :, :, :, :, :) = 0.d0
      End If
  ! begin loop over species
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
        ! calculate r^2
            r2 (ir) = spr (ir, is) ** 2
         End Do
     ! begin loop over atoms
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
        !--------------------!
        !     derivatives    !
        !--------------------!
        ! APW functions
            Do l1 = 0, input%groundstate%lmaxapw
               Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  Do io = 1, apword (l1, is)
                     fapw (:) = apwfr (:, 1, io, l1, ias)
                     Call gradzfmtr (input%groundstate%lmaxapw, nr, &
                    & spr(1, is), l1, m1, lmmaxapw, nrmtmax, fapw, &
                    & dapwfr(1, 1, 1, io, lm1))
                  End Do
               End Do
            End Do
            If (nlotot .Gt. 0) Then
           ! local orbital functions
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     flo (:) = lofr (:, 1, ilo, ias)
                     Call gradzfmtr (input%groundstate%lmaxapw, nr, &
                    & spr(1, is), l1, m1, lmmaxapw, nrmtmax, flo, &
                    & dlofr(1, 1, 1, ilo, m1))
                  End Do
               End Do
            End If
        !----------------!
        !     APW-APW    !
        !----------------!
            Do l1 = 0, input%groundstate%lmaxapw
               Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  Do io1 = 1, apword (l1, is)
                     Do l3 = 0, input%groundstate%lmaxapw
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io2 = 1, apword (l3, is)
                              Do j = 1, 3
                                 fr (:) = r2 (1:nr) * apwfr (1:nr, 1, &
                                & io1, l1, ias) * dapwfr (lm1, 1:nr, j, &
                                & io2, lm3)
                                 Call fderiv (-1, nr, spr(1, is), fr, &
                                & gr, cf)
                                 ripaa (io1, lm1, io2, lm3, ias, j) = &
                                & gr (nr)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
            If (nlotot .Gt. 0) Then
           !----------------------------!
           !     APW-local-orbital      !
           !----------------------------!
               Do l1 = 0, input%groundstate%lmaxapw
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     Do io1 = 1, apword (l1, is)
                        Do ilo = 1, nlorb (is)
                           l3 = lorbl (ilo, is)
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Do j = 1, 3
                                 fr (:) = r2 (1:nr) * apwfr (1:nr, 1, &
                                & io1, l1, ias) * dlofr (lm1, 1:nr, j, &
                                & ilo, m3)
                                 Call fderiv (-1, nr, spr(1, is), fr, &
                                & gr, cf)
                                 ripalo (io1, lm1, ilo, m3, ias, j) = &
                                & gr (nr)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
           !----------------------------!
           !     local-orbital-APW      !
           !----------------------------!
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     Do l3 = 0, input%groundstate%lmaxapw
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io2 = 1, apword (l3, is)
                              Do j = 1, 3
                                 fr (:) = r2 (1:nr) * lofr (:, 1, ilo, &
                                & ias) * dapwfr (lm1, 1:nr, j, io2, &
                                & lm3)
                                 Call fderiv (-1, nr, spr(1, is), fr, &
                                & gr, cf)
                                 riploa (ilo, m1, io2, lm3, ias, j) = &
                                & gr (nr)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
           !------------------------------------!
           !     local-orbital-local-orbital    !
           !------------------------------------!
               Do ilo1 = 1, nlorb (is)
                  l1 = lorbl (ilo1, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     Do ilo2 = 1, nlorb (is)
                        l3 = lorbl (ilo2, is)
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do j = 1, 3
                              fr (:) = r2 (1:nr) * lofr (:, 1, ilo1, &
                             & ias) * dlofr (lm1, 1:nr, j, ilo2, m3)
                              Call fderiv (-1, nr, spr(1, is), fr, gr, &
                             & cf)
                              riplolo (ilo1, m1, ilo2, m3, ias, j) = gr &
                             & (nr)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End If
        ! end loops over atoms and species
         End Do
      End Do
  ! deallocate
      Deallocate (fapw, dapwfr)
      If (nlotot .Gt. 0) deallocate (flo, dlofr)
End Subroutine pmatrad
