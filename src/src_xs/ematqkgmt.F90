!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqkgmt (iq, ik, igq)
      Use modmain
      Use modinput
      Use modxs
      Use m_zaxpyc
      Use m_xszoutpr
      Use m_xszoutpr3
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik, igq
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqkgmt'
      Integer :: is, ia, ias, l1, m1, lm1, l3, m3, lm3, io, io1, io2, &
     & ilo, ilo1, ilo2
      Integer :: lmax1, lmax3, ikt, i, j
      Complex (8), Allocatable :: zv (:)
      Allocate (zv(nstsv))
      ikt = ik
      lmax1 = input%xs%lmaxapwwf
      lmax3 = lmax1
      xih (:, :) = zzero
      xiuhloa (:, :) = zzero
      xiohalo (:, :) = zzero
      xiou (:, :, igq) = zzero
  ! loop over species and atoms
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call cpu_time (cmt0)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
        ! loop over (l',m',p')
            Do l1 = 0, lmax1
               Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  Do io1 = 1, apword (l1, is)
                     zv (:) = zzero
                 ! loop over (l'',m'',p'')
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io2 = 1, apword (l3, is)
                              Call zaxpy (nstsv, intrgaa(lm1, io1, lm3, &
                             & io2, ias), apwcmt(1, io2, lm3, ias), 1, &
                             & zv, 1)
                           End Do
                        End Do ! m3
                     End Do ! l3
                     Call xszoutpr (nst1, nst2, &
                    & fourpi*conjg(sfacgq(igq, ias, iq)), &
                    & apwcmt0(istl1:istu1, io1, lm1, ias), &
                    & zv(istl2:istu2), xiou(:, :, igq))
                 ! end loop over (l',m',p')
                  End Do ! io1
               End Do ! m1
            End Do ! l1
            Call cpu_time (cmt1)
            If (input%xs%fastemat) Then
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
           ! loop over local orbitals
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     zv (:) = zzero
                 ! loop over (l'',m'',p'')
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io = 1, apword (l3, is)
                              Call zaxpy (nstsv, intrgloa(m1, ilo, lm3, &
                             & io, ias), apwcmt(1, io, lm3, ias), 1, &
                             & zv, 1)
                           End Do ! io
                        End Do ! m3
                     End Do ! l3
                     Call xszoutpr (nst1, nst2, &
                    & fourpi*conjg(sfacgq(igq, ias, iq)), &
                    & locmt0(istl1:istu1, ilo, m1, ias), &
                    & zv(istl2:istu2), xiou(:, :, igq))
                  End Do ! m1
               End Do ! ilo
               Call cpu_time (cmt2)
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
           ! loop over (l'',m'',p'')
               Do l3 = 0, lmax3
                  Do m3 = - l3, l3
                     lm3 = idxlm (l3, m3)
                     Do io = 1, apword (l3, is)
                        zv (:) = zzero
                    ! loop over local orbitals
                        Do ilo = 1, nlorb (is)
                           l1 = lorbl (ilo, is)
                           Do m1 = - l1, l1
                              lm1 = idxlm (l1, m1)
                              Call zaxpy (nstsv, intrgalo(m1, ilo, lm3, &
                             & io, ias), locmt(1, ilo, m1, ias), 1, zv, &
                             & 1)
                           End Do ! m1
                        End Do ! ilo
                        Call xszoutpr (nst1, nst2, &
                       & fourpi*conjg(sfacgq(igq, ias, iq)), &
                       & apwcmt0(istl1:istu1, io, lm3, ias), &
                       & zv(istl2:istu2), xiou(:, :, igq))
                     End Do ! io
                  End Do ! m3
               End Do ! l3
               Call cpu_time (cmt3)
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
               Do ilo1 = 1, nlorb (is)
                  l1 = lorbl (ilo1, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     zv (:) = zzero
                     Do ilo2 = 1, nlorb (is)
                        l3 = lorbl (ilo2, is)
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Call zaxpy (nstsv, intrglolo(m1, ilo1, m3, &
                          & ilo2, ias), locmt(1, ilo2, m3, ias), 1, zv, &
                          & 1)
                        End Do ! m3
                     End Do ! ilo2
                     Call xszoutpr (nst1, nst2, &
                    & fourpi*conjg(sfacgq(igq, ias, iq)), &
                    & locmt0(istl1:istu1, ilo1, m1, ias), &
                    & zv(istl2:istu2), xiou(:, :, igq))
                  End Do ! m1
               End Do ! ilo1
            Else
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
           ! loop over local orbitals
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     i = idxlo (lm1, ilo, ias)
                     zv (:) = zzero
                 ! loop over (l'',m'',p'')
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io = 1, apword (l3, is)
                              Call zaxpy (nstsv, intrgloa(m1, ilo, lm3, &
                             & io, ias), apwcmt(1, io, lm3, ias), 1, &
                             & zv, 1)
                           End Do ! io
                        End Do ! m3
                     End Do ! l3
                     xiuhloa (i, :) = xiuhloa (i, :) + fourpi * conjg &
                    & (sfacgq(igq, ias, iq)) * zv (istl2:istu2)
                  End Do ! m1
               End Do ! ilo
               Call cpu_time (cmt2)
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
           ! loop over local orbitals
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     j = idxlo (lm1, ilo, ias)
                     zv (:) = zzero
                 ! loop over (l'',m'',p'')
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io = 1, apword (l3, is)
                              Call zaxpyc (nstsv, intrgalo(m1, ilo, &
                             & lm3, io, ias), apwcmt0(:, io, lm3, ias), &
                             & 1, zv, 1)
                           End Do ! io
                        End Do ! m3
                     End Do ! l3
                     xiohalo (:, j) = xiohalo (:, j) + fourpi * conjg &
                    & (sfacgq(igq, ias, iq)) * zv (istl1:istu1)
                  End Do ! m1
               End Do ! ilo
               Call cpu_time (cmt3)
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
               Do ilo1 = 1, nlorb (is)
                  l1 = lorbl (ilo1, is)
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     i = idxlo (lm1, ilo1, ias)
                     Do ilo2 = 1, nlorb (is)
                        l3 = lorbl (ilo2, is)
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           j = idxlo (lm3, ilo2, ias)
                           xih (i, j) = xih (i, j) + fourpi * conjg &
                          & (sfacgq(igq, ias, iq)) * intrglolo (m1, &
                          & ilo1, m3, ilo2, ias)
                        End Do ! m3
                     End Do ! ilo2
                  End Do ! m1
               End Do ! ilo1
           ! strategy of emat-calculation
            End If
            Call cpu_time (cmt4)
            cpumtaa = cpumtaa + cmt1 - cmt0
            cpumtloa = cpumtloa + cmt2 - cmt1
            cpumtalo = cpumtalo + cmt3 - cmt2
            cpumtlolo = cpumtlolo + cmt4 - cmt3
        ! end loop over species and atoms
         End Do ! ia
      End Do ! is
  ! deallocate
      Deallocate (zv)
End Subroutine ematqkgmt
