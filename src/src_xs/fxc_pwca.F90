!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
Module m_fxc_alda
      Implicit None
Contains
!
!
      Subroutine fxc_alda (iq, msiz, fxcg)
         Use modmain, Only: ivg, gc, fourpi, lmmaxvr, nrmtmax, natmtot, &
        & ngrtot, ngvec, ivgig
         Use modinput
         Use modxs, Only: unitout, fxcmt, gqc, fxcir, igqig
         Use m_ftfun, Only: ftfun
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, msiz
         Complex (8), Intent (Out) :: fxcg (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'fxc_alda'
         Complex (8), Allocatable :: fxcft1 (:)
         Integer :: sh (2), ig, ngmax, igq1, igq2, iv1 (3), iv (3)
         sh = shape (fxcg)
         If ((sh(1) .Lt. msiz) .Or. (sh(2) .Lt. msiz)) Then
            Write (unitout, '(a,2i9,a,i9,a)') 'Error(' // trim &
           & (thisnam) // '): size of fxc is to small (required)', sh, &
           & '(', msiz, ')'
            Stop
         End If
         If (allocated(fxcmt)) deallocate (fxcmt)
         If (allocated(fxcir)) deallocate (fxcir)
         Allocate (fxcmt(lmmaxvr, nrmtmax, natmtot))
         Allocate (fxcir(ngrtot))
    ! calculate exchange-correlation kernel in real space
         Call kernxc
    ! determine G-vector cutoff for 2*|G+q|_max
         Do ngmax = 1, ngvec
            If (2.d0*input%xs%gqmax .Lt. gc(ngmax)) Exit
         End Do
    ! Fourier transform of muffin-tin and interstitial kernel
         Allocate (fxcft1(ngmax))
         Call ftfun (ngmax, input%xs%tddft%lmaxalda, .True., .True., &
        & fxcir, fxcmt, fxcft1)
    ! transform G''=G-G' to G and G'
         Do igq1 = 1, msiz
            iv1 (:) = ivg (:, igqig(igq1, iq))
            Do igq2 = 1, msiz
               iv (:) = iv1 (:) - ivg (:, igqig(igq2, iq))
               ig = ivgig (iv(1), iv(2), iv(3))
               fxcg (igq1, igq2) = fxcft1 (ig)
          ! renormalization to symmetrized quantity wrt. G-space
               fxcg (igq1, igq2) = fxcg (igq1, igq2) * (gqc(igq1, &
              & iq)*gqc(igq2, iq)) / fourpi
            End Do
         End Do
    ! deallocate
         Deallocate (fxcmt, fxcir, fxcft1)
      End Subroutine fxc_alda
!
End Module m_fxc_alda
