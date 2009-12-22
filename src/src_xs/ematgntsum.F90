!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematgntsum (iq, igq)
      Use modmain
      Use modinput
      Use modxs
      Use m_findgntn0
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, igq
  ! local variables
      Integer :: is, ia, ias
      Integer :: l1, l2, l3, m2, lm2
      Integer :: ilo, ilo1, ilo2, io, io1, io2
      Integer :: lmax1, lmax2, lmax3, lmmax1, lmmax2, lmmax3
      Integer :: u1, u2, u3, u4
      Integer :: m1, m3, lm1, lm3, cl1, cm1, cl2, cm2, cl3, cm3
      lmax1 = Max (input%xs%lmaxapwwf, lolmax)
      lmax2 = input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same
      lmax3 = lmax1
      lmmax1 = (lmax1+1) ** 2
      lmmax2 = (lmax2+1) ** 2
      lmmax3 = (lmax3+1) ** 2
  ! allocate arrays for radial integrals and Bessel functions
      If (allocated(intrgaa)) deallocate (intrgaa)
      If (allocated(intrgloa)) deallocate (intrgloa)
      If (allocated(intrgalo)) deallocate (intrgalo)
      If (allocated(intrglolo)) deallocate (intrglolo)
      Allocate (intrgaa(lmmax1, apwordmax, lmmax3, apwordmax, natmtot))
      Allocate (intrgloa(-lolmax:lolmax, nlomax, lmmax3, apwordmax, &
     & natmtot))
      Allocate (intrgalo(-lolmax:lolmax, nlomax, lmmax3, apwordmax, &
     & natmtot))
      Allocate (intrglolo(-lolmax:lolmax, nlomax,-lolmax:lolmax, &
     & nlomax, natmtot))
  ! allocate temporary arrays
      intrgaa (:, :, :, :, :) = zzero
      intrgloa (:, :, :, :, :) = zzero
      intrgalo (:, :, :, :, :) = zzero
      intrglolo (:, :, :, :, :) = zzero
      If (input%xs%dbglev .Gt. 2) Then
     ! APW-APW
         Call getunit (u1)
         Open (Unit=u1, File='IRADGAUNTaa'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u1, '(a)') 'igq, ias, lm1, io1, lm3, io2,   intrgaa'
         Write (u1, '(a)') '-------------------------------------------&
        &-----------'
     ! APW-lo
         Call getunit (u2)
         Open (Unit=u2, File='IRADGAUNTalo'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u2, '(a)') 'igq, ias, m3, ilo, lm1, io,     intrgalo'
         Write (u2, '(a)') '-------------------------------------------&
        &-----------'
     ! lo-APW
         Call getunit (u3)
         Open (Unit=u3, File='IRADGAUNTloa'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u3, '(a)') 'igq, ias, m1, ilo, lm3, io,     intrgloa'
         Write (u3, '(a)') '-------------------------------------------&
        &-----------'
     ! lo-lo
         Call getunit (u4)
         Open (Unit=u4, File='IRADGAUNTlolo'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u4, '(a)') 'igq, ias, m1, ilo1, m3, ilo2,   intrglolo'
         Write (u4, '(a)') '-------------------------------------------&
        &-----------'
      End If
  ! begin loop over species
      Do is = 1, nspecies
     ! begin loop over atoms
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
        !---------------------------!
        !     APW-APW integrals     !
        !---------------------------!
            Do cl1 = 1, l1shape
               l1 = l1map (cl1)
               Do cm1 = 1, m1shape (l1)
                  m1 = m1map (l1, cm1)
                  lm1 = idxlm (l1, m1)
                  Do io1 = 1, apword (l1, is)
                     Do cl2 = 1, l2shape (l1, m1)
                        l3 = l2map (l1, m1, cl2)
                        Do cm2 = 1, m2shape (l1, m1, l3)
                           m3 = m2map (l1, m1, l3, cm2)
                           lm3 = idxlm (l3, m3)
                           Do io2 = 1, apword (l3, is)
                              Do cl3 = 1, l3shape (l1, m1, l3, m3)
                                 l2 = l3map (l1, m1, l3, m3, cl3)
                                 Do cm3 = 1, m3shape (l1, m1, l3, m3, &
                                & l2)
                                    m2 = m3map (l1, m1, l3, m3, l2, &
                                   & cm3)
                                    lm2 = idxlm (l2, m2)
                                    intrgaa (lm1, io1, lm3, io2, ias) = &
                                   & intrgaa (lm1, io1, lm3, io2, ias) &
                                   & + conjg (zil(l2)) * riaa (l1, io1, &
                                   & l3, io2, l2, ias, igq) * conjg &
                                   & (ylmgq(lm2, igq, iq)) * xsgnt &
                                   & (lm1, lm2, lm3)
                                 End Do
                              End Do
                              If (input%xs%dbglev .Gt. 2) Then
                                 Write (u1, '(6i5, 2g18.10)') igq, ias, &
                                & lm1, io1, lm3, io2, intrgaa (lm1, &
                                & io1, lm3, io2, ias)
                              End If
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
        !-------------------------------------!
        !     APW-local-orbital integrals     !
        !-------------------------------------!
            Do cl1 = 1, l1shape
               l1 = l1map (cl1)
               Do cm1 = 1, m1shape (l1)
                  m1 = m1map (l1, cm1)
                  lm1 = idxlm (l1, m1)
                  Do io = 1, apword (l1, is)
                     Do ilo = 1, nlorb (is)
                        l3 = lorbl (ilo, is)
                        Do cm2 = 1, m2shape (l1, m1, l3)
                           m3 = m2map (l1, m1, l3, cm2)
                           lm3 = idxlm (l3, m3)
                           Do cl3 = 1, l3shape (l1, m1, l3, m3)
                              l2 = l3map (l1, m1, l3, m3, cl3)
                              Do cm3 = 1, m3shape (l1, m1, l3, m3, l2)
                                 m2 = m3map (l1, m1, l3, m3, l2, cm3)
                                 lm2 = idxlm (l2, m2)
                                 intrgalo (m3, ilo, lm1, io, ias) = &
                                & intrgalo (m3, ilo, lm1, io, ias) + &
                                & conjg (zil(l2)) * riloa (ilo, l1, io, &
                                & l2, ias, igq) * conjg (ylmgq(lm2, &
                                & igq, iq)) * xsgnt (lm1, lm2, lm3)
                              End Do
                           End Do
                           If (input%xs%dbglev .Gt. 2) Then
                              Write (u2, '(6i5, 2g18.10)') igq, ias, &
                             & m3, ilo, lm1, io, intrgalo (m3, ilo, &
                             & lm1, io, ias)
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Do
        !-------------------------------------!
        !     local-orbital-APW integrals     !
        !-------------------------------------!
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do cm1 = 1, m1shape (l1)
                  m1 = m1map (l1, cm1)
                  lm1 = idxlm (l1, m1)
                  Do cl2 = 1, l2shape (l1, m1)
                     l3 = l2map (l1, m1, cl2)
                     Do cm2 = 1, m2shape (l1, m1, l3)
                        m3 = m2map (l1, m1, l3, cm2)
                        lm3 = idxlm (l3, m3)
                        Do io = 1, apword (l3, is)
                           Do cl3 = 1, l3shape (l1, m1, l3, m3)
                              l2 = l3map (l1, m1, l3, m3, cl3)
                              Do cm3 = 1, m3shape (l1, m1, l3, m3, l2)
                                 m2 = m3map (l1, m1, l3, m3, l2, cm3)
                                 lm2 = idxlm (l2, m2)
                                 intrgloa (m1, ilo, lm3, io, ias) = &
                                & intrgloa (m1, ilo, lm3, io, ias) + &
                                & conjg (zil(l2)) * riloa (ilo, l3, io, &
                                & l2, ias, igq) * conjg (ylmgq(lm2, &
                                & igq, iq)) * xsgnt (lm1, lm2, lm3)
                              End Do
                           End Do
                           If (input%xs%dbglev .Gt. 2) Then
                              Write (u3, '(6i5, 2g18.10)') igq, ias, &
                             & m1, ilo, lm3, io, intrgloa (m1, ilo, &
                             & lm3, io, ias)
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Do
        !-----------------------------------------------!
        !     local-orbital-local-orbital integrals     !
        !-----------------------------------------------!
            Do ilo1 = 1, nlorb (is)
               l1 = lorbl (ilo1, is)
               Do cm1 = 1, m1shape (l1)
                  m1 = m1map (l1, cm1)
                  lm1 = idxlm (l1, m1)
                  Do ilo2 = 1, nlorb (is)
                     l3 = lorbl (ilo2, is)
                     Do cm2 = 1, m2shape (l1, m1, l3)
                        m3 = m2map (l1, m1, l3, cm2)
                        lm3 = idxlm (l3, m3)
                        Do cl3 = 1, l3shape (l1, m1, l3, m3)
                           l2 = l3map (l1, m1, l3, m3, cl3)
                           Do cm3 = 1, m3shape (l1, m1, l3, m3, l2)
                              m2 = m3map (l1, m1, l3, m3, l2, cm3)
                              lm2 = idxlm (l2, m2)
                              intrglolo (m1, ilo1, m3, ilo2, ias) = &
                             & intrglolo (m1, ilo1, m3, ilo2, ias) + &
                             & conjg (zil(l2)) * rilolo (ilo1, ilo2, &
                             & l2, ias, igq) * conjg (ylmgq(lm2, igq, &
                             & iq)) * xsgnt (lm1, lm2, lm3)
                           End Do
                        End Do
                        If (input%xs%dbglev .Gt. 2) Then
                           Write (u4, '(6i5, 2g18.10)') igq, ias, m1, &
                          & ilo1, m3, ilo2, intrglolo (m1, ilo1, m3, &
                          & ilo2, ias)
                        End If
                     End Do
                  End Do
               End Do
            End Do
        ! end loops over atoms and species
         End Do
      End Do
  ! deallocate
      If (input%xs%dbglev .Gt. 2) Then
     ! close files
         Close (u1)
         Close (u2)
         Close (u3)
         Close (u4)
      End If
End Subroutine ematgntsum
