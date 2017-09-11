!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhovalk
! !INTERFACE:
!
!
Subroutine genrhoir (ik, evecfv, evecsv)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the partial valence charge density from the eigenvectors at
!   $k$-point {\tt ik}. In the muffin-tin region, the wavefunction is obtained
!   in terms of its $(l,m)$-components from both the APW and local-orbital
!   functions. Using a backward spherical harmonic transform (SHT), the
!   wavefunction is converted to real-space and the density obtained from its
!   modulus squared. This density is then transformed with a forward SHT and
!   accumulated in the global variable {\tt rhomt}. A similar proccess is used
!   for the intersitial density in which the wavefunction in real-space is
!   obtained from a Fourier transform of the sum of APW functions. The
!   interstitial density is added to the global array {\tt rhoir}. See routines
!   {\tt wavefmt}, {\tt genshtmat} and {\tt seceqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: nsd, ispn, jspn, is, ia, ias, ist
      Integer :: ir, irc, itp, igk, ifg, i, j, n
      Real (8) :: t1, t2, t3, t4
      Real (8) :: ts0, ts1
      Complex (8) zt1, zt2, zt3
! allocatable arrays
      Complex (8), Allocatable :: zfft (:, :)
!addidional arrays to make truncationerrors predictable
      Real (8) :: rhoir_k (ngrtot)
      Real (8) :: magir_k (ngrtot, ndmag)

      Call timesec (ts0)
      rhoir_k (:) = 0
      magir_k (:, :) = 0
      If (associated(input%groundstate%spin)) Then
         If (ncmag) Then
            nsd = 4
         Else
            nsd = 2
         End If
      Else
         nsd = 1
      End If
      Allocate (zfft(ngrtot, nspinor))
!------------------------------!
!     interstitial density     !
!------------------------------!
      Do j = 1, nstsv
         t1 = wkpt (ik) * occsv (j, ik)
         If (Abs(t1) .Gt. input%groundstate%epsocc) Then
            t2 = t1 / omega
            zfft (:, :) = 0.d0
            If (input%groundstate%tevecsv) Then
! generate spinor wavefunction from second-variational eigenvectors
               i = 0
               Do ispn = 1, nspinor
                  If (isspinspiral()) Then
                     jspn = ispn
                  Else
                     jspn = 1
                  End If
                  Do ist = 1, nstfv
                     i = i + 1
                     zt1 = evecsv (i, j)
                     If (Abs(dble(zt1))+Abs(aimag(zt1)) .Gt. &
                    & input%groundstate%epsocc) Then
                        Do igk = 1, ngk (jspn, ik)
                           ifg = igfft (igkig(igk, jspn, ik))
                           zfft (ifg, ispn) = zfft (ifg, ispn) + zt1 * &
                          & evecfv (igk, ist, jspn)
                        End Do
                     End If
                  End Do
               End Do
            Else
! spin-unpolarised wavefunction
               Do igk = 1, ngk (1, ik)
                  ifg = igfft (igkig(igk, 1, ik))
                  zfft (ifg, 1) = evecfv (igk, j, 1)
               End Do
            End If
            
! Fourier transform wavefunction to real-space
            Do ispn = 1, nspinor
               Call zfftifc (3, ngrid, 1, zfft(:, ispn))
            End Do
!
            If (associated(input%groundstate%spin)) Then
! spin-polarised
               Do ir = 1, ngrtot
                  zt1 = zfft (ir, 1)
                  zt2 = zfft (ir, 2)
                  zt3 = zt1 * conjg (zt2)
                  t3 = dble (zt1) ** 2 + aimag (zt1) ** 2
                  t4 = dble (zt2) ** 2 + aimag (zt2) ** 2
                  rhoir_k (ir) = rhoir_k (ir) + t2 * (t3+t4)
                  If (ncmag) Then
                     magir_k (ir, 1) = magir_k (ir, 1) + 2.d0 * t2 * &
                    & dble (zt3)
                     magir_k (ir, 2) = magir_k (ir, 2) - 2.d0 * t2 * &
                    & aimag (zt3)
                     magir_k (ir, 3) = magir_k (ir, 3) + t2 * (t3-t4)
                  Else
                     magir_k (ir, 1) = magir_k (ir, 1) + t2 * (t3-t4)
                  End If
               End Do
            Else
! spin-unpolarised
               Do ir = 1, ngrtot
                  zt1 = zfft (ir, 1)
                  rhoir_k (ir) = rhoir_k (ir) + t2 * &
                 & (dble(zt1)**2+aimag(zt1)**2)
               End Do
            End If
         End If
      End Do
      !write(*,*) ik
      !write(*,*) ik, ik
      !write(*,'(3F13.6)') vkl( :, ik)
      !write(*,'(SP,5F13.6)') rhoir_k( 1:100)
      !write(*,*)
      Deallocate (zfft)
!      Call timesec (ts1)
!$OMP CRITICAL
      rhoir (:) = rhoir (:) + rhoir_k (:)
!
      If (associated(input%groundstate%spin)) Then
         magir (:, :) = magir (:, :) + magir_k (:, :)
      End If
!      timerho = timerho + ts1 - ts0
!$OMP END CRITICAL
      Call timesec (ts1)
      timerho = timerho + ts1 - ts0
      Return
End Subroutine
!EOC
