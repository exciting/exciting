!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: energy
! !INTERFACE:
!
!
Subroutine energy
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Computes the total energy and its individual contributions. The kinetic
!   energy is given by
!   $$ T_{\rm s}=\sum_i n_i\epsilon_i-\int\rho({\bf r})[v_{\rm C}({\bf r})
!    +v_{\rm xc}({\bf r})]d{\bf r}-\int {\bf m}({\bf r})\cdot
!    ({\bf B}_{\rm xc}({\bf r})+{\bf B}_{\rm ext}({\bf r}))d{\bf r}, $$
!   where $n_i$ are the occupancies and $\epsilon_i$ are the eigenvalues of both
!   the core and valence states; $\rho$ is the density; ${\bf m}$ is the
!   magnetisation density; $v_{\rm C}$ is the Coulomb potential; $v_{\rm xc}$
!   and ${\bf B}_{\rm xc}$ are the exchange-correlation potential and effective
!   magnetic field, respectively; and ${\bf B}_{\rm ext}$ is the external
!   magnetic field. The Hartree, electron-nuclear and nuclear-nuclear
!   electrostatic energies are combined into the Coulomb energy:
!   \begin{align*}
!    E_{\rm C}&=E_{\rm H}+E_{\rm en}+E_{\rm nn} \\
!             &=\frac{1}{2}V_{\rm C}+E_{\rm Mad},
!   \end{align*}
!   where
!   $$ V_{\rm C}=\int\rho({\bf r})v_{\rm C}({\bf r})d{\bf r} $$
!   is the Coulomb potential energy. The Madelung energy is given by
!   $$ E_{\rm Mad}=\frac{1}{2}\sum_{\alpha}z_{\alpha}R_{\alpha}, $$
!   where
!   $$ R_{\alpha}=\lim_{r\rightarrow 0}\left(v^{\rm C}_{\alpha;00}(r)Y_{00}
!    +\frac{z_{\alpha}}{r}\right) $$
!   for atom $\alpha$, with $v^{\rm C}_{\alpha;00}$ being the $l=0$ component of
!   the spherical harmonic expansion of $v_{\rm C}$ in the muffin-tin, and
!   $z_{\alpha}$ is the nuclear charge. Using the nuclear-nuclear energy
!   determined at the start of the calculation, the electron-nuclear and Hartree
!   energies can be isolated with
!   $$ E_{\rm en}=2\left(E_{\rm Mad}-E_{\rm nn}\right) $$
!   and
!   $$ E_{\rm H}=\frac{1}{2}(E_{\rm C}-E_{\rm en}). $$
!   Finally, the total energy is
!   $$ E=T_{\rm s}+E_{\rm C}+E_{\rm xc}, $$
!   where $E_{\rm xc}$ is obtained either by integrating the
!   exchange-correlation energy density, or in the case of exact exchange, the
!   explicit calculation of the Fock exchange integral. The energy from the
!   external magnetic fields in the muffin-tins, {\tt bfcmt}, is always removed
!   from the total since these fields are non-physical: their field lines do not
!   close. The energy of the physical external field, {\tt bfieldc}, is also not
!   included in the total because this field, like those in the muffin-tins,
!   is used for breaking spin symmetry and taken to be infintesimal. If this
!   field is intended to be finite, then the associated energy, {\tt engybext},
!   should be added to the total by hand. See {\tt potxc}, {\tt exxengy} and
!   related subroutines.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ik, ist, idm, jdm
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0
      Real (8) :: vn
      Complex (8) zt1
! allocatable arrays
      Complex (8), Allocatable :: evecsv (:, :), c (:, :)
! external functions
      Real (8) :: rfmtinp, rfinp, rfint
      Complex (8) zdotc
      External rfmtinp, rfinp, rfint, zdotc
!-----------------------------------------------!
!     exchange-correlation potential energy     !
!-----------------------------------------------!
      engyvxc = rfinp (1, rhomt, vxcmt, rhoir, vxcir)
!-----------------------------------------------------!
!     exchange-correlation effective field energy     !
!-----------------------------------------------------!
      engybxc = 0.d0
      Do idm = 1, ndmag
         engybxc = engybxc + rfinp (1, magmt(:, :, :, idm), bxcmt(:, :, &
        & :, idm), magir(:, idm), bxcir(:, idm))
      End Do
!------------------------------------------!
!     external magnetic field energies     !
!------------------------------------------!
      engybext = 0.d0
      engybmt = 0.d0
      Do idm = 1, ndmag
         If (ncmag) Then
            jdm = idm
         Else
            jdm = 3
         End If
! energy of physical global field
         engybext = engybext + ga4 * momtot (idm) * &
        & input%groundstate%spin%bfieldc(jdm)
! energy of non-physical muffin-tin fields
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               engybmt = engybmt + ga4 * mommt (idm, ias) * input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(jdm)
            End Do
         End Do
      End Do
!----------------------------------!
!     Coulomb potential energy     !
!----------------------------------!
      engyvcl = rfinp (1, rhomt, vclmt, rhoir, vclir)
!-----------------------!
!     Madelung term     !
!-----------------------!
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
!---------------------------------------------!
!     electron-nuclear interaction energy     !
!---------------------------------------------!
      engyen = 2.d0 * (engymad-engynn)
!------------------------!
!     Hartree energy     !
!------------------------!
      engyhar = 0.5d0 * (engyvcl-engyen)
!------------------------!
!     Coulomb energy     !
!------------------------!
      engycl = engynn + engyen + engyhar
!-------------------------!
!     exchange energy     !
!-------------------------!
! exchange energy from the density
      engyx = rfinp (1, rhomt, exmt, rhoir, exir)
! exact exchange for OEP-EXX or Hartree-Fock on last iteration
      If ((input%groundstate%xctypenumber .Lt. 0) .Or. (task .Eq. 5) &
     & .Or. (task .Eq. 6)) Then
         If (tlast) Call exxengy
      End If
!----------------------------!
!     correlation energy     !
!----------------------------!
      engyc = rfinp (1, rhomt, ecmt, rhoir, ecir)
! zero correlation energy for Hartree-Fock
      If ((task .Eq. 5) .Or. (task .Eq. 6)) engyc = 0.d0
!----------------------!
!     LDA+U energy     !
!----------------------!
      engylu = 0.d0
      If (ldapu .Ne. 0) Then
         Do ias = 1, natmtot
            engylu = engylu + engyalu (ias)
         End Do
      End If
!-----------------------------------------------!
!     compensating background charge energy     !
!-----------------------------------------------!
      If (input%groundstate%chgexs .Ne. 0.d0) Then
         engycbc = input%groundstate%chgexs * rfint (vclmt, vclir)
      Else
         engycbc = 0.d0
      End If
!----------------------------!
!     sum of eigenvalues     !
!----------------------------!
! core eigenvalues
      evalsum = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
               If (spcore(ist, is)) evalsum = evalsum + spocc (ist, is) &
              & * evalcr (ist, ias)
            End Do
         End Do
      End Do
! valence eigenvalues
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            evalsum = evalsum + wkpt (ik) * occsv (ist, ik) * evalsv &
           & (ist, ik)
         End Do
      End Do
!------------------------!
!     kinetic energy     !
!------------------------!
! core electron kinetic energy
      Call energykncr
! total electron kinetic energy
      If ((task .Eq. 5) .Or. (task .Eq. 6)) Then
! Hartree-Fock case
         engykn = engykncr
! kinetic energy from valence states
         Allocate (evecsv(nstsv, nstsv))
         Allocate (c(nstsv, nstsv))
         Do ik = 1, nkpt
            Call getevecsv (vkl(:, ik), evecsv)
            Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, kinmatc(:, &
           & :, ik), nstsv, evecsv, nstsv, zzero, c, nstsv)
            Do ist = 1, nstsv
               zt1 = zdotc (nstsv, evecsv(:, ist), 1, c(:, ist), 1)
               engykn = engykn + wkpt (ik) * occsv (ist, ik) * dble &
              & (zt1)
            End Do
         End Do
         Deallocate (evecsv, c)
      Else
! Kohn-Sham case
         engykn = evalsum - engyvcl - engyvxc - engybxc - engybext - &
        & engybmt
      End If
!----------------------!
!     total energy     !
!----------------------!
      engytot = engykn + 0.5d0 * engyvcl + engymad + engyx + engyc + &
     & engycbc
! add the LDA+U correction if required
      If (ldapu .Ne. 0) engytot = engytot + engylu
      Return
End Subroutine
!EOC
