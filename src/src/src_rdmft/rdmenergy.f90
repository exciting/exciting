!
!
!
! Copyright (C) 2005-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmenergy
! calculate the total energy for RDMFT
      Use modinput
      Use modmain
      Implicit None
! local variables
      Integer :: is, ia, ias
      Integer :: ik, ist, ir
      Real (8) :: vn, t1
      Complex (8) zt1
! allocatable arrays
      Real (8), Allocatable :: rfmt (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: c (:, :)
! external functions
      Real (8) :: rfmtinp
      Complex (8) zdotc
      External rfmtinp, zdotc
      Allocate (rfmt(lmmaxvr, nrmtmax))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (c(nstsv, nstsv))
! Coulomb energy from core states
      engyvcl = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               rfmt (1, ir) = rhocr (ir, ias) / y00
            End Do
            engyvcl = engyvcl + rfmtinp (1, 0, nrmt(is), spr(:, is), &
           & lmmaxvr, rfmt, vclmt(:, :, ias))
         End Do
      End Do
      engykn = engykncr
      Do ik = 1, nkpt
         Call getevecsv (vkl(:, ik), evecsv)
         Do ist = 1, nstsv
            t1 = wkpt (ik) * occsv (ist, ik)
! Coulomb energy from valence states
            engyvcl = engyvcl + t1 * dble (vclmat(ist, ist, ik))
! kinetic energy from valence states
            zt1 = zdotc (nstsv, evecsv(:, ist), 1, dkdc(:, ist, ik), 1)
            engykn = engykn + t1 * dble (zt1)
         End Do
      End Do
! Madelung term
      engymad = 0.d0
      Do is = 1, nspecies
! compute the bare nucleus potential at the origin
         Call potnucl (input%groundstate%ptnucl, 1, spr(:, is), &
        & spzn(is), vn)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            engymad = engymad + 0.5d0 * spzn (is) * (vclmt(1, 1, &
           & ias)*y00-vn)
         End Do
      End Do
! exchange-correlation energy
      Call rdmengyxc
! total energy
      engytot = 0.5d0 * engyvcl + engymad + engykn + engyx
      If (input%groundstate%RDMFT%rdmtemp .Gt. 0.d0) Then
         Call rdmentropy
         engytot = engytot - input%groundstate%RDMFT%rdmtemp * &
        & rdmentrpy
      End If
      Write (*,*)
      Write (*, '("Info(rdmengy): Energies :")')
      Write (*, '(" electronic kinetic", T30, ": ", G18.10)') engykn
      Write (*, '(" core electron kinetic", T30, ": ", G18.10)') &
     & engykncr
      Write (*, '(" Coulomb", T30, ": ", G18.10)') engyvcl
      Write (*, '(" Madelung", T30, ": ", G18.10)') engymad
      Write (*, '(" exchange-correlation", T30, ": ", G18.10)') engyx
      If (input%groundstate%RDMFT%rdmtemp .Gt. 0.d0) Then
         Write (*, '(" entropy", T30, ": ", G18.10)') &
        & input%groundstate%RDMFT%rdmtemp * rdmentrpy
      End If
      Write (*, '(" total", T30, ": ", G18.10)') engytot
      Write (*,*)
      Deallocate (evecsv, rfmt, c)
      Return
End Subroutine
