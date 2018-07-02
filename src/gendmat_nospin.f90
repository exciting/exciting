!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gendmat_nospin ( fst, lst, lmin, lmax, is, ia, ngp, apwalm, evecfv, ld, dmat)
      Use modinput
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: lmin, fst, lst
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv)
      Complex (8), Intent (In) :: evecfv( nmatmax, fst:lst, nspnfv)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: dmat( ld, ld, fst:lst)
! local variables
      Integer :: lmmax
      Integer :: l, m1, m2, lm1, lm2
      Integer :: i, j, n, ist, irc
      Real (8) :: t1, t2
      Complex (8) zt1
! automatic arrays
      Real (8) :: fr1 (nrcmtmax), fr2 (nrcmtmax)
      Real (8) :: gr (nrcmtmax), cf (3, nrcmtmax)
! allocatable arrays
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      If (lmin .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(gendmat): lmin < 0 : ", I8)') lmin
         Write (*,*)
         Stop
      End If
      If (lmax .Gt. input%groundstate%lmaxapw) Then
         Write (*,*)
         Write (*, '("Error(gendmat): lmax > lmaxapw : ", 2I8)') lmax, &
        & input%groundstate%lmaxapw
         Write (*,*)
         Stop
      End If
      lmmax = (lmax+1) ** 2
! allocate local arrays
      Allocate( wfmt2( lmmax, nrcmtmax, nspinor))
! zero the density matrix
      dmat = zzero
      n = lmmax * nrcmt (is)
! begin loop over second-variational states
      Do j = fst, lst
! spin-unpolarised wavefunction
         Call wavefmt( input%groundstate%lradstep, lmax, is, ia, ngp, apwalm, evecfv(:, j, 1), lmmax, wfmt2)
         Do l = lmin, lmax
            Do m1 = - l, l
               lm1 = idxlm (l, m1)
               Do m2 = - l, l
                  lm2 = idxlm (l, m2)
                  Do irc = 1, nrcmt (is)
                     zt1 = wfmt2( lm1, irc, 1)*conjg( wfmt2( lm2, irc, 1))
                     t1 = rcmt( irc, is)**2
                     fr1( irc) = dble( zt1)*t1
                     fr2( irc) = aimag( zt1)*t1
                  End Do
                  Call fderiv (-1, nrcmt(is), rcmt(:, is), fr1, gr, cf)
                  t1 = gr( nrcmt( is))
                  Call fderiv (-1, nrcmt(is), rcmt(:, is), fr2, gr, cf)
                  t2 = gr( nrcmt( is))
                  dmat( lm1, lm2, j) = cmplx( t1, t2, 8)
               End Do
            End Do
         End Do
! end loop over second-variational states
      End Do
      Deallocate( wfmt2)
      Return
End Subroutine
