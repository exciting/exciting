!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematgntsum (iq, igq, integrals)
      Use modmain
      Use modinput
      Use modxs
      Use m_findgntn0
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, igq
      complex(8) :: integrals(apwmaxsize+lomaxsize,apwmaxsize+lomaxsize,natmtot)
  ! local variables
      Integer :: is, ia, ias,iaug1,iaug2
      Integer :: l1, l2, l3, m2, lm2
      Integer :: ilo, ilo1, ilo2, io, io1, io2
      Integer :: lmax1, lmax2, lmax3, lmmax1, lmmax2, lmmax3
      Integer :: u1, u2, u3, u4
      Integer :: m1, m3, lm1, lm3, cl1, cm1, cl2, cm2, cl3, cm3
     
      complex(8), allocatable :: intrgaa(:,:,:,:),intrgloa(:,:,:,:),intrglolo(:,:,:,:),intrgalo(:,:,:,:)

!      print *,'howdy'
      lmax1 = Max (input%xs%lmaxapwwf, lolmax)
      lmax2 = input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same
      lmax3 = lmax1
      lmmax1 = (lmax1+1) ** 2
      lmmax2 = (lmax2+1) ** 2
      lmmax3 = (lmax3+1) ** 2
  ! allocate arrays for radial integrals and Bessel functions
      Allocate (intrgaa(lmmax1, apwordmax, lmmax3, apwordmax))
      Allocate (intrgloa(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
      Allocate (intrgalo(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
      Allocate (intrglolo(-lolmax:lolmax, nlomax,-lolmax:lolmax, nlomax))


!      If (allocated(integrals%aa)) deallocate (integrals%aa)
!      If (allocated(integrals%alo)) deallocate (integrals%alo)
!      If (allocated(integrals%loa)) deallocate (integrals%loa)
!      If (allocated(integrals%lolo)) deallocate (integrals%lolo)
!      Allocate (integrals%aa(apwmaxsize,apwmaxsize, natmtot))
!      Allocate (integrals%alo(lomaxsize, apwmaxsize, natmtot))
!      Allocate (integrals%loa(apwmaxsize, lomaxsize, natmtot))
!      Allocate (integrals%lolo(lomaxsize, lomaxsize, natmtot))
  ! allocate temporary arrays
      integrals=zzero
!      integrals%aa(:, :, :) = zzero
!      integrals%alo( :, :, :) = zzero
!      integrals%loa(:, :, :) = zzero
!      integrals%lolo(:, :, :) = zzero
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
            intrgaa (:, :, :, :) = zzero
            intrgloa (:, :, :, :) = zzero
            intrgalo (:, :, :, :) = zzero
            intrglolo (:, :, :, :) = zzero
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
                                    intrgaa (lm1, io1, lm3, io2) = &
                                   & intrgaa (lm1, io1, lm3, io2) &
                                   & + conjg (zil(l2)) * riaa (l1, io1, &
                                   & l3, io2, l2, ias, igq) * conjg &
                                   & (ylmgq(lm2, igq, iq)) * xsgnt &
                                   & (lm1, lm2, lm3)
                                 End Do
                              End Do
                              If (input%xs%dbglev .Gt. 2) Then
                                 Write (u1, '(6i5, 2g18.10)') igq, ias, &
                                & lm1, io1, lm3, io2, intrgaa (lm1, &
                                & io1, lm3, io2)
                              End If
                           End Do
                        End Do
                     End Do
                  End Do


               End Do
            End Do

          iaug2=0
          do l3=0,input%xs%lmaxapwwf
            do m3=-l3,l3
              do io2=1,apword(l3,is)
                iaug2=iaug2+1
                lm3=idxlm(l3,m3)

                iaug1=0
                do l1=0,input%xs%lmaxapwwf
                  do m1=-l1,l1
                      do io1=1,apword(l1,is)
                        iaug1=iaug1+1
                        lm1=idxlm(l1,m1)
                        integrals(iaug2,iaug1,ias)=intrgaa (lm1, io1, lm3, io2)
                      enddo
                    enddo
                  enddo


              enddo
            enddo
          enddo

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
                                 intrgalo (m3, ilo, lm1, io) = &
                                & intrgalo (m3, ilo, lm1, io) + &
                                & conjg (zil(l2)) * riloa (ilo, l1, io, &
                                & l2, ias, igq) * conjg (ylmgq(lm2, &
                                & igq, iq)) * xsgnt (lm1, lm2, lm3)
                              End Do
                           End Do
                           If (input%xs%dbglev .Gt. 2) Then
                              Write (u2, '(6i5, 2g18.10)') igq, ias, &
                             & m3, ilo, lm1, io, intrgalo (m3, ilo, &
                             & lm1, io)
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Do

          iaug2=0
          Do ilo = 1, nlorb (is)
            l3 = lorbl (ilo, is)
            do m3=-l3,l3
              iaug2=iaug2+1
              lm3=idxlm(l3,m3)
             
              iaug1=0
              do l1=0,input%xs%lmaxapwwf
                do m1=-l1,l1
                  do io=1,apword(l1,is)
                    iaug1=iaug1+1
                    lm1=idxlm(l1,m1)
                    integrals(apwsize(is)+iaug2,iaug1,ias)=intrgalo (m3, ilo, lm1, io)
                  enddo
                enddo
               enddo

            enddo
          enddo
!write(*,*) sum(intrgalo2)
!write(*,*) sum(intrgalo)




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
                                 intrgloa (m1, ilo, lm3, io) = &
                                & intrgloa (m1, ilo, lm3, io) + &
                                & conjg (zil(l2)) * riloa (ilo, l3, io, &
                                & l2, ias, igq) * conjg (ylmgq(lm2, &
                                & igq, iq)) * xsgnt (lm1, lm2, lm3)
                              End Do
                           End Do
                           If (input%xs%dbglev .Gt. 2) Then
                              Write (u3, '(6i5, 2g18.10)') igq, ias, &
                             & m1, ilo, lm3, io, intrgloa (m1, ilo, &
                             & lm3, io)
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Do

          iaug2=0
          Do ilo = 1, nlorb (is)
            l3 = lorbl (ilo, is)
            do m3=-l3,l3
              iaug2=iaug2+1
              lm3=idxlm(l3,m3)



              iaug1=0
              do l1=0,input%xs%lmaxapwwf
                do m1=-l1,l1
                  do io=1,apword(l1,is)
                    iaug1=iaug1+1
                    lm1=idxlm(l1,m1)
                    integrals(iaug1,apwsize(is)+iaug2,ias)=intrgloa (m3, ilo, lm1, io)
                  enddo
                enddo
               enddo

            enddo
          enddo
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
                              intrglolo (m1, ilo1, m3, ilo2) = &
                             & intrglolo (m1, ilo1, m3, ilo2) + &
                             & conjg (zil(l2)) * rilolo (ilo1, ilo2, &
                             & l2, ias, igq) * conjg (ylmgq(lm2, igq, &
                             & iq)) * xsgnt (lm1, lm2, lm3)
                           End Do
                        End Do
                        If (input%xs%dbglev .Gt. 2) Then
                           Write (u4, '(6i5, 2g18.10)') igq, ias, m1, &
                          & ilo1, m3, ilo2, intrglolo (m1, ilo1, m3, &
                          & ilo2)
                        End If
                     End Do
                  End Do
               End Do
            End Do

          iaug2=0
          Do ilo2 = 1, nlorb (is)
            l3 = lorbl (ilo2, is)
            do m3=-l3,l3
              iaug2=iaug2+1



              iaug1=0
              Do ilo1 = 1, nlorb (is)
                l1 = lorbl (ilo1, is)
                do m1=-l1,l1
                  iaug1=iaug1+1
                  integrals(apwsize(is)+iaug2,apwsize(is)+iaug1,ias)=intrglolo (m1, ilo1, m3, ilo2)
                enddo
              enddo

            enddo
          enddo

        ! end loops over atoms and species
         End Do
      End Do

      If (allocated(intrgaa)) deallocate (intrgaa)
      If (allocated(intrgloa)) deallocate (intrgloa)
      If (allocated(intrgalo)) deallocate (intrgalo)
      If (allocated(intrglolo)) deallocate (intrglolo)
  ! deallocate
      If (input%xs%dbglev .Gt. 2) Then
     ! close files
         Close (u1)
         Close (u2)
         Close (u3)
         Close (u4)
      End If
      
End Subroutine ematgntsum
