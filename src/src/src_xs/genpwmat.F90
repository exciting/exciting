!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genpwmat
! !INTERFACE:
!
!
Subroutine genpwmat (vpl, ngpmax, ngp, vgpc, gpc, igpig, ylmgp, sfacgp, &
& vklk, ngkk, igkigk, apwalmk, evecfvk, evecsvk, vklkp, ngkkp, igkigkp, &
& apwalmkp, evecfvkp, evecsvkp, pwmat)
! !USES:
      Use modmain
      Use modxs
      Use modinput
! !INPUT/OUTPUT PARAMETERS:
! !DESCRIPTION:
!   Calculates the matrix elements of the plane wave
!   $$ M_{ij{\bf k}}=\langle\Psi_{i,{\bf k}}|e^{-i({\bf G}+{\bf q}){\bf r}}|
!      \Psi_{j,{\bf k}}\rangle. $$
!  Straightforward implementation for checking.
!
! !REVISION HISTORY:
!   Created November 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vpl (3)
      Integer, Intent (In) :: ngpmax
      Integer, Intent (In) :: ngp
      Real (8), Intent (In) :: vgpc (3, ngpmax)
      Real (8), Intent (In) :: gpc (ngpmax)
      Integer, Intent (In) :: igpig (ngpmax)
      Complex (8), Intent (In) :: ylmgp (lmmaxapw, ngpmax)
      Complex (8), Intent (In) :: sfacgp (ngpmax, natmtot)
      Real (8), Intent (In) :: vklk (3)
      Integer, Intent (In) :: ngkk
      Integer, Intent (In) :: igkigk (ngkmax)
      Complex (8), Intent (In) :: apwalmk (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: evecfvk (nmatmax, nstfv)
      Complex (8), Intent (In) :: evecsvk (nstsv, nstsv)
      Real (8), Intent (In) :: vklkp (3)
      Integer, Intent (In) :: ngkkp
      Integer, Intent (In) :: igkigkp (ngkmax)
      Complex (8), Intent (In) :: apwalmkp (ngkmax, apwordmax, &
     & lmmaxapw, natmtot)
      Complex (8), Intent (In) :: evecfvkp (nmatmax, nstfv)
      Complex (8), Intent (In) :: evecsvkp (nstsv, nstsv)
      Complex (8), Intent (Out) :: pwmat (ngp, nstsv, nstsv)
  ! local variables
      Integer :: i, j, k, l, m, lm, irc, igp, ig, ir, iv (3), ivu (3), &
     & ispn, is, ia, ias, ist, jst
      Integer :: igp1, igp2, ig1, ig2, iv1 (3)
      Real (8) :: v1 (3), t1
      Complex (8) zt1, zt2
  ! allocatable arrays
      Real (8), Allocatable :: jlgpr (:, :)
      Complex (8), Allocatable :: wfmtk (:, :, :)
      Complex (8), Allocatable :: wfmtkp (:, :, :)
      Complex (8), Allocatable :: wfmt1 (:, :)
      Complex (8), Allocatable :: wfmt2 (:, :)
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: pwfmt (:, :)
      Complex (8), Allocatable :: wfirk (:, :)
      Complex (8), Allocatable :: wfirkp (:, :)
      Complex (8), Allocatable :: pm (:, :, :)
      Complex (8), Allocatable :: cfunt (:, :), h (:, :), pmt (:, :)
      Complex (8), Allocatable :: evecfvt1 (:, :), evecfvt2 (:, :)
  ! external functions
      Real (8), External :: r3taxi
      Complex (8), External :: zfmtinp
      t1 = vgpc (1, 1)! for convenience:)
  ! allocate arrays
      Allocate (cfunt(ngkk, ngkkp))
      Allocate (h(ngkk, nstfv))
      Allocate (pmt(nstfv, nstfv))
      Allocate (evecfvt1(nstfv, ngkk), evecfvt2(ngkkp, nstfv))
  ! check if q-point is commensurate with k-mesh
      If (any(Abs(vpl*input%groundstate%ngridk-&
     & Nint(vpl*input%groundstate%ngridk)) .Gt. &
     & input%structure%epslat)) Then
         Write (*,*)
         Write (*, '("Error(genpwmat): q-point not commensurate with k-&
        &mesh : ",3g18.10)') vpl
         Write (*,*)
         Stop
      End If
  ! kp-point umklapp G-vector
      v1 (:) = vpl (:) + vklk (:)
      Call r3frac (input%structure%epslat, v1, ivu)
  ! check k+q=kp+Gw
      If (r3taxi(v1, vklkp) .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(genpwmat): q-point not commensurate with k-&
        &point and kp-point")')
         Write (*, '(" k-point         : ",3g18.10)') vklk
         Write (*, '(" q-point         : ",3g18.10)') vpl
         Write (*, '(" kp-point        : ",3g18.10)') vklkp
         Write (*, '(" umklapp G-vector: ",3i8)') ivu
         Write (*,*)
         Stop
      End If
      Allocate (wfmtk(lmmaxapw, nrcmtmax, nstfv))
      Allocate (wfmtkp(lmmaxapw, nrcmtmax, nstfv))
      Allocate (wfmt1(lmmaxapw, nrcmtmax))
      Allocate (wfmt2(lmmaxapw, nrcmtmax))
      Allocate (zfmt(lmmaxapw, nrcmtmax))
      Allocate (pwfmt(lmmaxapw, nrcmtmax))
      Allocate (wfirk(ngrtot, nstfv))
      Allocate (wfirkp(ngrtot, nstfv))
      Allocate (pm(ngp, nstfv, nstfv))
      Allocate (jlgpr(0:input%groundstate%lmaxapw, nrcmtmax))
  ! zero arrays
      wfmtk (:, :, :) = zzero
      wfmtkp (:, :, :) = zzero
      wfmt1 (:, :) = zzero
      wfmt2 (:, :) = zzero
      zfmt (:, :) = zzero
      pwfmt (:, :) = zzero
      pm (:, :, :) = zzero
  ! loop over G+p vectors
      Do igp = 1, ngp
         Write (*,*) 'Info(genpwmat): igp: ', igp
     ! calculate matrix elements of the plane wave in the muffin-tin
         Do is = 1, nspecies
        ! calculate bessel functions j_l(|G+q||r|)
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Call sbessel (input%groundstate%lmaxapw, &
              & gpc(igp)*spr(ir, is), jlgpr(0, irc))
           ! set up plane wave factor from Rayleigh formula
               Do l = 0, input%xs%lmaxemat
                  Do m = - l, l
                     lm = idxlm (l, m)
                     pwfmt (lm, irc) = fourpi * conjg (zil(l)) * jlgpr &
                    & (l, irc) * conjg (ylmgp(lm, igp))
                  End Do
               End Do
            End Do
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ist = 1, nstfv
              ! calculate the wavefunction for k-point
                  Call wavefmt (input%groundstate%lradstep, &
                 & input%groundstate%lmaxapw, is, ia, ngkk, apwalmk, &
                 & evecfvk(1, ist), lmmaxapw, wfmtk(1, 1, ist))
              ! calculate the wavefunction for kp-point
                  Call wavefmt (input%groundstate%lradstep, &
                 & input%groundstate%lmaxapw, is, ia, ngkkp, apwalmkp, &
                 & evecfvkp(1, ist), lmmaxapw, wfmtkp(1, 1, ist))
              ! convert wavefunction and plane wave to spherical coordinates
                  Call zgemm ('N', 'N', lmmaxapw, nrcmt(is), lmmaxapw, &
                 & zone, zbshtapw, lmmaxapw, wfmtkp(1, 1, ist), &
                 & lmmaxapw, zzero, wfmt1, lmmaxapw)
                  Call zgemm ('N', 'N', lmmaxapw, nrcmt(is), lmmaxapw, &
                 & zone, zbshtapw, lmmaxapw, pwfmt, lmmaxapw, zzero, &
                 & wfmt2, lmmaxapw)
              ! direct evaluation of exp(-i(G+q)r) on spherical grid
              ! note that the Rayleigh expansion seems to be more precise
              !do l=0,lmaxapw
              !   do m=-l,l
              !      lm=idxlm(l,m)
              !      irc=0
              !      do ir=1,nrmt(is),lradstp
              !         irc=irc+1
              !         vr(:)=spr(ir,is)*sphcov(:,lm) *** sphcover must be
              !         ! invoked to obtain "sphcov(:,lm)" ***
              !         t1=dot_product(vgpc(:,igp),vr)
              !         wfmt2(lm,irc)=cmplx(cos(t1),-sin(t1),8)
              !      end do
              !   end do
              !end do
              ! calculate product in muffin-tin in real space
                  Do irc = 1, nrcmt (is)
                     zfmt (:, irc) = wfmt1 (:, irc) * wfmt2 (:, irc)
                  End Do
              ! convert to spherical harmonics
                  Call zgemm ('N', 'N', lmmaxapw, nrcmt(is), lmmaxapw, &
                 & zone, zfshtapw, lmmaxapw, zfmt, lmmaxapw, zzero, &
                 & wfmtkp(1, 1, ist), lmmaxapw)
               End Do
           ! structure factor for G+p vector
               zt2 = conjg (sfacgp(igp, ias))
               Do ist = 1, nstfv
                  Do jst = 1, nstfv
                     zt1 = zfmtinp (.True., input%groundstate%lmaxapw, &
                    & nrcmt(is), rcmt(1, is), lmmaxapw, wfmtk(1, 1, &
                    & ist), wfmtkp(1, 1, jst))
                     pm (igp, ist, jst) = pm (igp, ist, jst) + zt1 * &
                    & zt2
                  End Do
               End Do
           ! end loops over atoms and species
            End Do
         End Do
     ! calculate matrix elements of the plane wave in the interstitial region
         Forall (ist=1:nstfv)
            evecfvt1 (ist, :) = conjg (evecfvk(1:ngkk, ist))
         End Forall
         evecfvt2 (:, :) = evecfvkp (1:ngkkp, :)
         Do igp1 = 1, ngkk
            ig1 = igkigk (igp1)
            iv1 (:) = ivg (:, ig1)
            Do igp2 = 1, ngkkp
               ig2 = igkigkp (igp2)
               iv (:) = iv1 (:) - ivg (:, ig2) + ivg (:, igpig(igp)) + &
              & ivu (:)
               ig = ivgig (iv(1), iv(2), iv(3))
               cfunt (igp1, igp2) = cfunig (ig)
            End Do
         End Do
	 ! h = cfunt*evecfvt2
         Call zgemm ('n', 'n', ngkk, nstfv, ngkkp, zone, cfunt, ngkk, &
        & evecfvt2, ngkkp, zzero, h, ngkk)
     ! pmt = evecfvt2*h
         Call zgemm ('n', 'n', nstfv, nstfv, ngkk, zone, evecfvt1, &
        & nstfv, h, ngkk, zzero, pmt, nstfv)
         pm (igp, :, :) = pm (igp, :, :) + pmt (:, :)
     ! compute the second-variational matrix elements of the plane wave
         If (input%groundstate%tevecsv) Then
            Do i = 1, nstsv
               Do j = 1, nstsv
                  zt1 = 0.d0
                  k = 0
                  Do ispn = 1, nspinor
                     Do ist = 1, nstfv
                        k = k + 1
                        l = (ispn-1) * nstfv
                        Do jst = 1, nstfv
                           l = l + 1
                           zt1 = conjg (evecsvk(k, i)) * evecsvkp (l, &
                          & j)
                           zt1 = zt1 + zt1 * pm (igp, ist, jst)
                        End Do
                     End Do
                  End Do
                  pwmat (igp, i, j) = zt1
               End Do
            End Do
         Else
            pwmat (:, :, :) = pm (:, :, :)
         End If
     ! end loop over G+p vectors
      End Do
      Deallocate (wfmtk, wfmtkp, wfmt1, wfmt2, zfmt, pwfmt, wfirk, &
     & wfirkp, pm, jlgpr)
      Deallocate (cfunt, h, pmt, evecfvt1, evecfvt2)
End Subroutine genpwmat
!EOC
