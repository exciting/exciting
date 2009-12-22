!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematrad (iq)
      Use modmain
      Use modinput
      Use modxs
      Use m_getunit
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Integer :: is, ia, ias, nr, ir, igq
      Integer :: l1, l2, l3
      Integer :: ilo, ilo1, ilo2, io, io1, io2
      Real (8) :: t1
      Integer :: lmax1, lmax2, lmax3
      Integer :: u11, u22, u33
  ! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
  ! allocatable arrays
      Real (8), Allocatable :: jl (:, :), jhelp (:)
!
      lmax1 = Max (input%xs%lmaxapwwf, lolmax)
      lmax2 = input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same!
      lmax3 = lmax1
!
  ! allocate arrays for radial integrals and Bessel functions
      If (allocated(riaa)) deallocate (riaa)
      If (allocated(riloa)) deallocate (riloa)
      If (allocated(rilolo)) deallocate (rilolo)
      Allocate (riaa(0:lmax1, apwordmax, 0:lmax3, apwordmax, 0:lmax2, &
     & natmtot, ngq(iq)))
      Allocate (riloa(nlomax, 0:lmax3, apwordmax, 0:lmax2, natmtot, &
     & ngq(iq)))
      Allocate (rilolo(nlomax, nlomax, 0:lmax2, natmtot, ngq(iq)))
  ! allocate temporary arrays
      Allocate (jl(0:lmax2, nrmtmax))
      Allocate (jhelp(0:lmax2))
      jl (:, :) = 0.d0
      jhelp (:) = 0.d0
  ! zero arrays for radial integrals
      riaa (:, :, :, :, :, :, :) = 0.d0
      riloa (:, :, :, :, :, :) = 0.d0
      rilolo (:, :, :, :, :) = 0.d0
!
      If (input%xs%dbglev .Gt. 1) Then
     ! APW-APW
         Call getunit (u11)
         Open (Unit=u11, File='IRADaa'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u11, '(a)') 'igq, ias, l1, io1, l3, io2, l2	 iraa'
         Write (u11, '(a)') '------------------------------------------&
        &-----------'
     ! lo-APW
         Call getunit (u22)
         Open (Unit=u22, File='IRADalo'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u22, '(a)') 'igq, ias, ilo, l1, l3, io, l2,	 iralo'
         Write (u22, '(a)') '------------------------------------------&
        &-----------'
     ! lo-lo
         Call getunit (u33)
         Open (Unit=u33, File='IRADlolo'//filext, Form='formatted', &
        & Action='write', Status='replace')
         Write (u33, '(a)') 'igq, ias, ilo1, l1, ilo2, l3, l2,   irlolo&
        &'
         Write (u33, '(a)') '------------------------------------------&
        &-----------'
      End If
!
  ! begin loop over G+q vectors
      Do igq = 1, ngq (iq)
     ! begin loop over species
         Do is = 1, nspecies
            nr = nrmt (is)
            Do ir = 1, nr
           ! calculate r^2
               r2 (ir) = spr (ir, is) ** 2
           ! calculate spherical Bessel functions of first kind j_l(|G+q|r_a)
               Call sbessel (lmax2, gqc(igq, iq)*spr(ir, is), jhelp)
               jl (:, ir) = jhelp (:)
            End Do
        ! begin loop over atoms
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
           !----------------!
           !     APW-APW    !
           !----------------!
               Do l1 = 0, lmax1
                  Do io1 = 1, apword (l1, is)
                     Do l3 = 0, lmax3
                        Do io2 = 1, apword (l3, is)
                           Do l2 = 0, lmax2
                              Do ir = 1, nr
                                 t1 = apwfr (ir, 1, io1, l1, ias) * &
                                & apwfr (ir, 1, io2, l3, ias) * r2 (ir)
                                 fr (ir) = t1 * jl (l2, ir)
                              End Do
                              Call fderiv (-1, nr, spr(1, is), fr, gr, &
                             & cf)
                              riaa (l1, io1, l3, io2, l2, ias, igq) = &
                             & gr (nr)
                              If (input%xs%dbglev .Gt. 1) Then
                                 Write (u11, '(7i5, g18.10)') igq, ias, &
                                & l1, io1, l3, io2, l2, gr (nr)
                              End If
                           End Do
                        End Do ! io2
                     End Do ! l3
                  End Do ! io1
               End Do ! l1
           !----------------------------!
           !     local-orbital-APW      !
           !----------------------------!
               Do ilo = 1, nlorb (is)
                  l1 = lorbl (ilo, is)
                  Do l3 = 0, lmax3
                     Do io = 1, apword (l3, is)
                        Do l2 = 0, lmax2
                           Do ir = 1, nr
                              t1 = lofr (ir, 1, ilo, ias) * apwfr (ir, &
                             & 1, io, l3, ias) * r2 (ir)
                              fr (ir) = t1 * jl (l2, ir)
                           End Do
                           Call fderiv (-1, nr, spr(1, is), fr, gr, cf)
                           riloa (ilo, l3, io, l2, ias, igq) = gr (nr)
                           If (input%xs%dbglev .Gt. 1) Then
                              Write (u22, '(7i5, g18.10)') igq, ias, &
                             & ilo, l1, l3, io, l2, gr (nr)
                           End If
                        End Do ! l2
                     End Do ! io
                  End Do ! l3
               End Do ! ilo
           !------------------------------------!
           !     local-orbital-local-orbital    !
           !------------------------------------!
               Do ilo1 = 1, nlorb (is)
                  l1 = lorbl (ilo1, is)
                  Do ilo2 = 1, nlorb (is)
                     l3 = lorbl (ilo2, is)
                     Do l2 = 0, lmax2
                        Do ir = 1, nr
                           t1 = lofr (ir, 1, ilo1, ias) * lofr (ir, 1, &
                          & ilo2, ias) * r2 (ir)
                           fr (ir) = t1 * jl (l2, ir)
                        End Do
                        Call fderiv (-1, nr, spr(1, is), fr, gr, cf)
                        rilolo (ilo1, ilo2, l2, ias, igq) = gr (nr)
                        If (input%xs%dbglev .Gt. 1) Then
                           Write (u33, '(7i5, g18.10)') igq, ias, ilo1, &
                          & l1, ilo2, l3, l2, gr (nr)
                        End If
                     End Do ! l2
                  End Do ! ilo2
               End Do ! ilo1
           ! end loops over atoms and species
            End Do
         End Do
     ! end loop over G+q vectors
      End Do
!
  ! deallocate
      Deallocate (jl, jhelp)
      If (input%xs%dbglev .Gt. 1) Then
     ! close files
         Close (u11)
         Close (u22)
         Close (u33)
      End If
End Subroutine ematrad
