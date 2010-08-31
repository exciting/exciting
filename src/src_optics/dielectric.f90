!
!
!
! Copyright (C) 2002-2008 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dielectric
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: ik, jk, isym
      Integer :: ist, jst, iw, i, j, l
      Integer :: recl, nsk (3)
      Real (8) :: eji, wd (2), wplas, t1, t2
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Complex (8) zv (3), eta, zt1
      Character (256) :: fname
! allocatable arrays
      Integer, Allocatable :: lspl (:)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: f (:, :)
      Real (8), Allocatable :: delta (:, :, :)
      Complex (8), Allocatable :: pmat (:, :, :)
      Complex (8), Allocatable :: sigma (:)
! initialise universal variables
      Call init0
      Call init1
! read Fermi energy from file
      Call readfermi
      Do ik = 1, nkpt
! get the eigenvalues and occupancies from file
         Call getevalsv (vkl(:, ik), evalsv(:, ik))
         Call getoccsv (vkl(:, ik), occsv(:, ik))
      End Do
! allocate local arrays
      Allocate (lspl(nkptnr))
      Allocate (w(input%properties%dos%nwdos))
      If (input%properties%dielectric%intraband) allocate (f(nstsv, nkpt))
      If (input%properties%dielectric%usegdft) allocate (delta(nstsv, nstsv, nkpt))
      Allocate (pmat(3, nstsv, nstsv))
      Allocate (sigma(input%properties%dos%nwdos))
! compute generalised DFT correction
      If (input%properties%dielectric%usegdft) Then
         Call readstate
         Call poteff
         Call linengy
         Call genapwfr
         Call genlofr
         Do ik = 1, nkpt
            Call gdft (ik, delta(:, :, ik))
         End Do
      End If
! energy interval should start from zero
      wdos (1) = 0.d0
! generate energy grid
      t1 = (wdos(2)-wdos(1)) / dble (input%properties%dos%nwdos)
      Do iw = 1, input%properties%dos%nwdos
         w (iw) = t1 * dble (iw-1) + wdos (1)
      End Do
! find crystal symmetries which map non-reduced k-points to reduced equivalents
      Do ik = 1, nkptnr
         Call findkpt (vklnr(:, ik), isym, jk)
         lspl (ik) = lsplsymc (isym)
      End Do
! find the record length for momentum matrix element file
      Inquire (IoLength=Recl) pmat
      Open (50, File='PMAT.OUT', Action='READ', Form='UNFORMATTED', &
     & Access='DIRECT', Recl=Recl)
! i divided by the complex relaxation time
      eta = cmplx (0.d0, input%groundstate%swidth)
! loop over dielectric tensor components
      Do l = 1, size(input%properties%dielectric%optcomp,2)
         i = input%properties%dielectric%optcomp(1, l)
         j = input%properties%dielectric%optcomp(2, l)
         sigma (:) = 0.d0
         If (input%properties%dielectric%intraband) f (:, :) = 0.d0
! loop over non-reduced k-points
         Do ik = 1, nkptnr
! equivalent reduced k-point
            jk = ikmap (ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
! read momentum matrix elements from direct-access file
            Read (50, Rec=jk) pmat
! valance states
            Do ist = 1, nstsv
               If (evalsv(ist, jk) .Lt. efermi) Then
! conduction states
                  Do jst = 1, nstsv
                     If (evalsv(jst, jk) .Gt. efermi) Then
! rotate the matrix elements from the reduced to non-reduced k-point
! (note that the inverse operation is used)
                        v1 (:) = dble (pmat(:, ist, jst))
                        Call r3mv (symlatc(:, :, lspl(ik)), v1, v2)
                        v1 (:) = aimag (pmat(:, ist, jst))
                        Call r3mv (symlatc(:, :, lspl(ik)), v1, v3)
                        zv (:) = cmplx (v2(:), v3(:), 8)
                        zt1 = occmax * zv (i) * conjg (zv(j))
                        eji = evalsv (jst, jk) - evalsv (ist, jk) + &
                       & input%properties%dielectric%scissor
                        If (input%properties%dielectric%usegdft) eji = eji + delta (jst, &
                       & ist, jk)
                        t1 = 1.d0 / (eji+input%groundstate%swidth)
                        Do iw = 1, input%properties%dos%nwdos
                           sigma (iw) = sigma (iw) + t1 * (zt1/(w(iw)-&
                          & eji+eta)+conjg(zt1)/(w(iw)+eji+eta))
                        End Do
                     End If
                  End Do
               End If
            End Do
         End Do
         zt1 = zi / (omega*dble(nkptnr))
         sigma (:) = zt1 * sigma (:)
! intraband contribution
         If (input%properties%dielectric%intraband) Then
            If (i .Eq. j) Then
! compute plasma frequency
               Do ik = 1, nkpt
                  Do ist = 1, nstsv
                     zt1 = pmat (i, ist, ist)
                     f (ist, ik) = dble (zt1) ** 2 + aimag (zt1**2)
                  End Do
               End Do
               wd (1) = efermi - input%groundstate%swidth
               wd (2) = efermi + input%groundstate%swidth
! number of subdivisions used for interpolation
               nsk (:) = Max &
              & (input%properties%dos%ngrdos/input%groundstate%ngridk(:&
              & ), 1)
               Call brzint (0, input%groundstate%ngridk, nsk, ikmap, 1, &
              & wd, nstsv, nstsv, evalsv, f, wplas)
               wplas = Abs (wplas) * occmax * 4.d0 * pi / omega
               wplas = Sqrt (wplas)
! write the plasma frequency to file
               Write (fname, '("PLASMA_", 2I1, ".OUT")') i, j
               Open (60, File=trim(fname), Action='WRITE', Form='FORMAT&
              &TED')
               Write (60, '(G18.10, " : plasma frequency")') wplas
               Close (60)
! add the intraband contribution to sigma
               t1 = wplas ** 2 / fourpi
               Do iw = 1, input%properties%dos%nwdos
                  sigma (iw) = sigma (iw) + t1 / &
                 & (input%groundstate%swidth-zi*w(iw))
               End Do
            End If
         End If
! write the optical conductivity to file
         Write (fname, '("SIGMA_", 2I1, ".OUT")') i, j
         Open (60, File=trim(fname), Action='WRITE', Form='FORMATTED')
         Do iw = 1, input%properties%dos%nwdos
            Write (60, '(2G18.10)') w (iw), dble (sigma(iw))
         End Do
         Write (60, '("     ")')
         Do iw = 1, input%properties%dos%nwdos
            Write (60, '(2G18.10)') w (iw), aimag (sigma(iw))
         End Do
         Close (60)
! write the dielectric function to file
         Write (fname, '("EPSILON_", 2I1, ".OUT")') i, j
         Open (60, File=trim(fname), Action='WRITE', Form='FORMATTED')
         t1 = 0.d0
         If (i .Eq. j) t1 = 1.d0
         Do iw = 1, input%properties%dos%nwdos
            If (w(iw) .Gt. 1.d-8) Then
               t2 = t1 - fourpi * aimag (sigma(iw)/(w(iw)+eta))
               Write (60, '(2G18.10)') w (iw), t2
            End If
         End Do
         Write (60, '("     ")')
         Do iw = 1, input%properties%dos%nwdos
            If (w(iw) .Gt. 1.d-8) Then
               t2 = fourpi * dble (sigma(iw)/(w(iw)+eta))
               Write (60, '(2G18.10)') w (iw), t2
            End If
         End Do
         Close (60)
! end loop over tensor components
      End Do
      Write (*,*)
      Write (*, '("Info(dielectric):")')
      Write (*, '(" dielectric tensor written to EPSILON_ij.OUT")')
      Write (*, '(" optical conductivity written to SIGMA_ij.OUT")')
      If (input%properties%dielectric%intraband) Then
         Write (*, '(" plasma frequency written to PLASMA_ij.OUT")')
      End If
      Write (*, '(" for components")')
      Do l = 1, size(input%properties%dielectric%optcomp,2)
         Write (*, '("  i = ", I1, ", j = ", I1)') &
        & input%properties%dielectric%optcomp(1:2, l)
      End Do
      Write (*,*)
      Deallocate (lspl, w, pmat, sigma)
      If (input%properties%dielectric%intraband) deallocate (f)
      If (input%properties%dielectric%usegdft) deallocate (delta)
      Return
End Subroutine
