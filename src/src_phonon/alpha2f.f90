!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine alpha2f
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: n, ik, iq, i, j
      Integer :: i1, i2, i3, iw
      Integer :: lwork, info
      Real (8) :: wmin, wmax, wd, dw, wlog
      Real (8) :: v (3), lambda, tc, t1
! allocatable arrays
      Real (8), Allocatable :: wq (:, :)
      Real (8), Allocatable :: wp (:)
      Real (8), Allocatable :: gq (:, :)
      Real (8), Allocatable :: a2fp (:)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: a2f (:)
      Real (8), Allocatable :: f (:), g (:), cf (:, :)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: dynp (:, :)
      Complex (8), Allocatable :: dynr (:, :, :)
      Complex (8), Allocatable :: ev (:, :), b (:, :)
      Complex (8), Allocatable :: a2fmq (:, :, :)
      Complex (8), Allocatable :: a2fmr (:, :, :)
      Complex (8), Allocatable :: a2fmp (:, :)
      Complex (8), Allocatable :: work (:)
! initialise universal variables
      Call init0
      Call init1
      Call init2
      n = 3 * natmtot
      Allocate (wq(n, nqpt))
      Allocate (wp(n))
      Allocate (gq(n, nqpt))
      Allocate (a2fp(n))
      Allocate (w(input%properties%dos%nwdos))
      Allocate (a2f(input%properties%dos%nwdos))
      Allocate (f(input%properties%dos%nwdos), &
     & g(input%properties%dos%nwdos), cf(3, &
     & input%properties%dos%nwdos))
      Allocate (rwork(3*n))
      Allocate (dynq(n, n, nqpt))
      Allocate (dynp(n, n))
      Allocate (dynr(n, n, ngridq(1)*ngridq(2)*ngridq(3)))
      Allocate (ev(n, n), b(n, n))
      Allocate (a2fmq(n, n, nqpt))
      Allocate (a2fmr(n, n, ngridq(1)*ngridq(2)*ngridq(3)))
      Allocate (a2fmp(n, n))
      lwork = 2 * n
      Allocate (work(lwork))
! get the eigenvalues and occupancies from file
      Do ik = 1, nkpt
         Call getevalsv (vkl(:, ik), evalsv(:, ik))
         Call getoccsv (vkl(:, ik), occsv(:, ik))
      End Do
! compute the density of states at the Fermi energy
      Call occupy
! read in the dynamical matrices
      Call readdyn (.true.,dynq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! Fourier transform the dynamical matrices to real-space
      Call dynqtor (dynq, dynr)
! read in the phonon linewidths for each q-point
      Call readgamma (gq)
! loop over phonon q-points
      Do iq = 1, nqpt
! diagonalise the dynamical matrix
         Call dyndiag (dynq(:, :, iq), wq(:, iq), ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues are the phonon linewidths divided by the frequency
         Do i = 1, n
            If (wq(i, iq) .Gt. 1.d-8) Then
               t1 = gq (i, iq) / wq (i, iq)
            Else
               t1 = 0.d0
            End If
            Do j = 1, n
               b (i, j) = t1 * conjg (ev(j, i))
            End Do
         End Do
         Call zgemm ('N', 'N', n, n, n, zone, ev, n, b, n, zzero, &
        & a2fmq(:, :, iq), n)
      End Do
! Fourier transform the matrices to real-space
      Call dynqtor (a2fmq, a2fmr)
! find the minimum and maximum frequencies
      wmin = 0.d0
      wmax = 0.d0
      Do iq = 1, nqpt
         wmin = Min (wmin, wq(1, iq))
         wmax = Max (wmax, wq(n, iq))
      End Do
      wmax = wmax + (wmax-wmin) * 0.1d0
      wmin = wmin - (wmax-wmin) * 0.1d0
      wd = wmax - wmin
      If (wd .Lt. 1.d-8) wd = 1.d0
      dw = wd / dble (input%properties%dos%nwdos)
! generate energy grid
      Do iw = 1, input%properties%dos%nwdos
         w (iw) = dw * dble (iw-1) + wmin
      End Do
      a2f (:) = 0.d0
      Do i1 = 0, input%properties%dos%ngrdos - 1
         v (1) = dble (i1) / dble (input%properties%dos%ngrdos)
         Do i2 = 0, input%properties%dos%ngrdos - 1
            v (2) = dble (i2) / dble (input%properties%dos%ngrdos)
            Do i3 = 0, input%properties%dos%ngrdos - 1
               v (3) = dble (i3) / dble (input%properties%dos%ngrdos)
! compute the dynamical matrix at this particular q-point
               Call dynrtoq (v, dynr, dynp)
! find the phonon frequencies
               Call dyndiag (dynp, wp, ev)
! compute the alpha^2F matrix at this particular q-point
               Call dynrtoq (v, a2fmr, a2fmp)
! diagonlise the alpha^2F matrix
               Call zheev ('N', 'U', n, a2fmp, n, a2fp, work, lwork, &
              & rwork, info)
               Do i = 1, n
                  t1 = (wp(i)-wmin) / dw + 1.d0
                  iw = Nint (t1)
                  If ((iw .Ge. 1) .And. (iw .Le. &
                 & input%properties%dos%nwdos)) Then
                     a2f (iw) = a2f (iw) + a2fp (i)
                  End If
               End Do
            End Do
         End Do
      End Do
      t1 = twopi * (fermidos/2.d0) * dw * dble &
     & (input%properties%dos%ngrdos) ** 3
      If (t1 .Gt. 1.d-8) Then
         t1 = 1.d0 / t1
      Else
         t1 = 0.d0
      End If
      a2f (:) = t1 * a2f (:)
! smooth Eliashberg function if required
      If (input%properties%dos%nsmdos .Gt. 0) Call fsmooth &
     & (input%properties%dos%nsmdos, input%properties%dos%nwdos, 1, &
     & a2f)
! write Eliashberg function to file
      Open (50, File='ALPHA2F.OUT', Action='WRITE', Form='FORMATTED')
      Do iw = 1, input%properties%dos%nwdos
         Write (50, '(2G18.10)') w (iw), a2f (iw)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(alpha2f):")')
      Write (*, '(" Eliashberg function written to ALPHA2F.OUT")')
! compute the total lambda
      Do iw = 1, input%properties%dos%nwdos
         If (w(iw) .Gt. 1.d-8) Then
            f (iw) = a2f (iw) / w (iw)
         Else
            f (iw) = 0.d0
         End If
      End Do
      Call fderiv (-1, input%properties%dos%nwdos, w, f, g, cf)
      lambda = 2.d0 * g (input%properties%dos%nwdos)
      Open (50, File='LAMBDA.OUT', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("Electron-phonon mass enhancement parameter, lambda &
     &: ", G18.10)') lambda
      Close (50)
      Write (*,*)
      Write (*, '("Info(alpha2f):")')
      Write (*, '(" Electron-phonon mass enhancement parameter, lambda,&
     & written to LAMBDA.OUT")')
! compute the logarithmic average frequency
      Do iw = 1, input%properties%dos%nwdos
         If (w(iw) .Gt. 1.d-8) Then
            f (iw) = a2f (iw) * Log (w(iw)) / w (iw)
         Else
            f (iw) = 0.d0
         End If
      End Do
      Call fderiv (-1, input%properties%dos%nwdos, w, f, g, cf)
      t1 = (2.d0/lambda) * g (input%properties%dos%nwdos)
      wlog = Exp (t1)
! compute McMillan-Allen-Dynes superconducting critical temperature
      t1 = (-1.04d0*(1.d0+lambda)) / (lambda-&
     & input%properties%eliashberg%mustar-&
     & 0.62d0*lambda*input%properties%eliashberg%mustar)
      tc = (wlog/(1.2d0*kboltz)) * Exp (t1)
      Open (50, File='TC_MCMILLAN.OUT', Action='WRITE', Form='FORMATTED&
     &')
      Write (50,*)
      Write (50, '("Logarithmic average frequency : ", G18.10)') wlog
      Write (50,*)
      Write (50, '("Coulomb pseudopotential, mu* : ", G18.10)') &
     & input%properties%eliashberg%mustar
      Write (50,*)
      Write (50, '("McMillan-Allen-Dynes superconducting critical tempe&
     &rature")')
      Write (50, '("[Eq. 34, Phys. Rev. B 12, 905 (1975)] (kelvin) : ",&
     & G18.10)') tc
      Write (50,*)
      Close (50)
      Write (*,*)
      Write (*, '("Info(alpha2f):")')
      Write (*, '(" McMillan-Allen-Dynes superconducting critical tempe&
     &rature, T_c,  written to TC_MCMILLAN.OUT")')
      Write (*,*)
      Deallocate (wq, wp, gq, a2fp, w, a2f, f, g, cf)
      Deallocate (rwork, dynq, dynp, dynr, ev, b)
      Deallocate (a2fmq, a2fmr, a2fmp, work)
      Return
End Subroutine
