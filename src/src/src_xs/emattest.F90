!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine emattest
      Use modmain
      Use modxs
      Use m_getunit
      Use m_getpmat
      Use m_getemat
      Use m_genfilname
      Implicit None
      Complex (8), Allocatable :: pmat (:, :, :, :), x (:, :, :, :)
      Real (8), Allocatable :: d (:, :, :), scis12 (:, :), scis21 (:, &
     & :)
  ! compare first Q-point
      Integer, Parameter :: iq = 1
      Complex (8) :: x_sc, p_sc
      Real (8) :: d1, d2, d3, a, p
      Integer :: n, ik, ikq, ist1, ist2
      Character (256) :: filename
      Logical, External :: tqgamma
      Call init0
      Call init1
      Call init2
      Call xssave0
      n = ngq (iq)
      If (tqgamma(iq)) Then
         Write (*,*)
         Write (*, '("Info(emattest): Q-point is Gamma point - no compa&
        &rison possible")')
         Write (*,*)
         Call terminate
      End If
  ! file extension for q-point
      Call genfilname (iqmt=iq, setfilext=.True.)
  ! calculate k+q and G+k+q related variables
      Call init1offs (vql(1, iq))
      Call findocclims (iq, istocc0, istocc, istunocc0, istunocc, &
     & isto0, isto, istu0, istu)
  ! allocate arrays
      If (allocated(deou)) deallocate (deou)
      If (allocated(deuo)) deallocate (deuo)
      Allocate (deou(nst1, nst2))
      Allocate (deuo(nst2, nst1))
      If (allocated(docc12)) deallocate (docc12)
      If (allocated(docc21)) deallocate (docc21)
      Allocate (docc12(nst1, nst2))
      Allocate (docc21(nst2, nst1))
      If (allocated(xiou)) deallocate (xiou)
      If (allocated(xiuo)) deallocate (xiuo)
      Allocate (xiou(nst1, nst2, n))
      Allocate (xiuo(nst2, nst1, n))
  ! allocate local arrays
      Allocate (x(nst1, nst2, n, nkpt))
      Allocate (d(nst1, nst2, nkpt))
      Allocate (pmat(3, nstsv, nstsv, nkpt))
      Allocate (scis12(nst1, nst2))
      Allocate (scis21(nst2, nst1))
      Call getunit (unit1)
      Call genfilname (basename='emat_pmat', iqmt=iq, filnam=filename)
      Open (unit1, File=trim(filename), Action='write', Status='replace&
     &')
  ! annotate magnitude of q-vector
      Write (*,*) 'Info(emattest): length of q-vector:', gqc (1, iq)
  ! test matrix elements
      Do ik = 1, nkpt
     ! read matrix elemets of exponential expression
         Call getpmat (ik, vkl0, 1, nstsv, 1, nstsv, .True., &
        & trim(fnpmat), pmat(:, :, :, ik))
     ! read matrix elemets of exponential expression
         Call getemat (iq, ik, .True., trim(fnpmat), ngq(iq), istl1, &
        & istu1, istl2, istu2, xiou, istl3, istu3, istl4, istu4, xiuo)
         ikq = ikmapikq (ik, iq)
         Call getdevaldoccsv (iq, ik, ikq, istl1, istu1, istl2, istu2, &
        & deou, docc12, scis12)
         Call getdevaldoccsv (iq, ik, ikq, istl2, istu2, istl1, istu1, &
        & deuo, docc21, scis21)
         x (:, :, :, ik) = xiou (:, :, :)
         d (:, :, ik) = deou
         Do ist1 = 1, nst1
            Do ist2 = 1, nst2
               x_sc = x (ist1, ist2, 1, ik) / gqc (1, iq)
               p_sc = dot_product (vgqc(:, 1, iq)/gqc(1, iq), pmat(:, &
              & istl1-1+ist1, istl2-1+ist2, ik)) / (-d(ist1, ist2, ik))
               a = dble (x_sc) ** 2 + aimag (x_sc) ** 2
               p = dble (p_sc) ** 2 + aimag (p_sc) ** 2
               d1 = Abs (x_sc-p_sc)
               d2 = Min (d1, Abs(x_sc+p_sc))
               d3 = Abs (x_sc) - Abs (p_sc)
               Write (unit1, '(100g18.10)') ik, ist1, ist2, x_sc, p_sc, &
              & d1, d2, d3, a, p
            End Do
         End Do
      End Do
      Deallocate (deou, deuo, docc12, docc21, scis12, scis21, xiou, &
     & xiuo, pmat, x, d)
      Close (unit1)
      Close (unit3)
      Close (unit4)
End Subroutine emattest
