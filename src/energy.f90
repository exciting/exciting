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
      Use mod_hybrids, only: ihyb, exnl
! !DESCRIPTION:
!   The {\tt energy} subroutine computes the total energy and its individual contributions. 
!   The total energy is composed of kinetic, Coulomb, and exchange-correlation energy,
!   %
!   \begin{equation}
!    E_{\rm tot}\,=\,T_{\rm s}\,+\,E_{C}\,+\,E_{\rm xc}.
!   \end{equation}
!   %
!   The kinetic energy of the non-interacting system is given by
!   %
!   \begin{equation}
!    T_{\rm s} = \sum_i n_i\epsilon_i \, - \, V_{\rm eff},
!   \label{kinetic}
!   \end{equation}
!   %
!   where $n_i$ are the occupancies and $\epsilon_i$ are the eigenvalues of both the core and 
!   valence states. The effective potential energy, $V_{\rm eff}$, can be expressed as
!   %
!   \begin{eqnarray}
!    V_{\rm eff}\,&=&\,\int\rho({\bf r}) \, v_{\rm C}({\bf r}) \, d{\bf r} + \int\rho({\bf r}) \, v_{\rm xc}({\bf r})\,d{\bf r} \nonumber \\
!                 &+&\int {\bf m}({\bf r})\cdot\left[{\bf B}_{\rm xc}({\bf r})+{\bf B}_{\rm ext}({\bf r})\right]\,d{\bf r}.
!   \label{Eeff}
!   \end{eqnarray}
!   %
!   The first and second term of Eq.~(\ref{Eeff}) are the Coulomb potential energy, $V_{C}$, and 
!   exchange-correlation potential energy, $V_{\rm xc}$, respectively. ${\bf m}({\bf r})$ is the 
!   magnetization density, and ${\bf B}_{\rm xc}$ and ${\bf B}_{\rm ext}$ are the 
!   exchange-correlation effective magnetic and the external magnetic fields, respectively.
!
!   The Coulomb energy consists of the Hartree energy, $E_{\rm H}$, the electron-nuclear energy, $
!   E_{\rm en}$, and the nuclear-nuclear energy, $E_{\rm nn}$,
!   %
!   \begin{eqnarray}
!    E_{\rm C}\,&=&\,E_{\rm H}\,+\,E_{\rm en}\,+\,E_{\rm nn} \nonumber \\
!               &=&\,(\underbrace{E_{\rm H}\,+\,\frac{1}{2}E_{\rm en}})\,+\,(\underbrace{\frac{1}{2}E_{\rm en}\,+\,E_{\rm nn}}) \nonumber \\
!               &=&\, \hspace{8mm} \frac{1}{2}V_{\rm C} \hspace{9.5mm} + \hspace{9mm} E_{\rm Madelung}.
!   \label{Eq4}
!   \end{eqnarray}
!   %
!   The Madelung energy is given by 
!   %
!   \begin{eqnarray}
!   E_{\rm Madelung}=\frac{1}{2}\sum_{\alpha}z_{\alpha}R_{\alpha}, 
!   \end{eqnarray}
!   where for each atom $\alpha$ with nuclear charge $z_{\alpha}$
!   %
!   \begin{eqnarray}
!    R_{\alpha}=\lim_{r\rightarrow 0}\left(v^{\rm C}_{\alpha,00}(r)Y_{00} +\frac{z_{\alpha}}{r}\right)
!   \end{eqnarray}
!   %
!   with $v^{\rm C}_{\alpha,00}$ being the $l=0$ component of the spherical harmonic expansion of 
!   $v_{\rm C}$ in the muffin-tin region. Using Eq.~(\ref{Eq4}), the electron-nuclear and Hartree 
!   energies can be expressed as
!   %
!   \begin{eqnarray}
!    E_{\rm en}=2\left(E_{\rm Madelung}-E_{\rm nn}\right)
!   \end{eqnarray}
!   %
!   and 
!   %
!   \begin{eqnarray}
!    E_{\rm H}=\frac{1}{2}(V_{\rm C}-E_{\rm en}).
!   \end{eqnarray}
!   %
!   $E_{\rm xc}$ is obtained either by integrating the exchange-correlation energy density,
!   %
!   \begin{eqnarray}
!    E_{\rm xc}\,=\, \int \rho({\bf r})\,\epsilon_{\rm xc}({\bf r})\,d{\bf r},
!   \end{eqnarray}
!   %
!   or in the case of exact exchange, the explicit calculation of the Fock exchange integral.
!
!   The energy from the external magnetic fields in the muffin-tins, {\tt bfcmt}, is always 
!   removed from the total since these fields are non-physical: their field lines do not close. 
!   The energy of the physical external field, {\tt bfieldc}, is also not included in the total 
!   because this field, like those in the muffin-tins, is used for breaking spin symmetry and 
!   taken to be infinitesimal. If this field is intended to be finite, then the associated energy, 
!   {\tt engybext}, should be added to the total by hand.
!   See {\tt potxc}, {\tt exxengy} and related subroutines.
!
!   !REVISION HISTORY:
!   Created Jun 2013
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ik, ist, idm, jdm, ir, lm
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0
      Real (8) :: vn
      Real (8) :: v2(50)
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
!      rhomt(2:lmmaxvr,:,1)=0d0
!      vclmt(2:lmmaxvr,:,1)=0d0
      engyvcl = rfinp (1, rhomt, vclmt, rhoir, vclir)
!      rhomt(2:lmmaxvr,:,1)=0d0
!      vclmt(2:lmmaxvr,:,1)=0d0
!     do lm=1,1
!      do ir=1,nrmt(1)
!        write(*,*) spr(ir,1),rhomt(lm,ir,1),vclmt(lm,ir,1)
!      enddo
!      read(*,*)
!     enddo
!      write(*,*) engyvcl
!      stop
!-----------------------!
!     Madelung term     !
!-----------------------!
if (.false.) then 
      engymad = 0.d0
      Do is = 1, nspecies
! compute the bare nucleus potential at the origin
         Call potnucl (input%groundstate%ptnucl, 1, spr(:, is), spzn(is), vn)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            write(*,*) vclmt(1,1,ias)*y00,vn
            engymad = engymad + 0.5d0 * spzn (is) * (vclmt(1, 1, &
           & ias)*y00-vn)
         End Do
      End Do
else
      engymad = 0.d0
      Do is = 1, nspecies
! compute the bare nucleus potential at the origin
         Call potnucl (input%groundstate%ptnucl, 50, spr(:, is), spzn(is), v2)
!        do ia=1,49
!          write(*,*) spr(ia,is),0.5d0 * spzn (is) * ((vclmt(1,ia,1)*y00-v2(ia))*spr(ia+1,is)-(vclmt(1,ia+1,1)*y00-v2(ia+1))*spr(ia,is))/(spr(ia+1,is)-spr(ia,is)),0.5d0 * spzn (is) * (vclmt(1,ia,1)*y00-v2(ia))
!        enddo
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
!            write(*,*) vclmt(1,1,ias)*y00,v2(1)
!            write(*,*) vclmt(1,2,ias)*y00,v2(2)
!            write(*,*) 0.5d0 * spzn (is) * ((vclmt(1,1,ias)*y00-v2(1))*spr(2,is)-(vclmt(1,2,ias)*y00-v2(2))*spr(1,is))/(spr(2,is)-spr(1,is))
!            write(*,*) 0.5d0 * spzn (is) * (vclmt(1,1,ias)*y00-v2(1)),0.5d0 * spzn (is) * vmad(ias)
            engymad = engymad + 0.5d0 * spzn (is) * vmad(ias) !((vclmt(1,1,ias)*y00-v2(1))*spr(2,is)-(vclmt(1,2,ias)*y00-v2(2))*spr(1,is))/(spr(2,is)-spr(1,is))
         End Do
      End Do
endif

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
!write(*,*) 
!write(*,*) 'engymad',engymad
!write(*,*) 'engynn',engynn
!write(*,*) 'engyvcl',engyvcl
!write(*,*) 'engyen',engyen

!-------------------------!
!     exchange energy     !
!-------------------------!
! exchange energy from the density
      engyx = rfinp (1, rhomt, exmt, rhoir, exir)
! Hartree-Fock energy
      If ((task .Eq. 5).Or.(task .Eq. 6)) Then
         If (tlast) Call exxengy
      End If
! calculate exchange energy for OEP-EXX/Hybrids on last iteration
       If (associated(input%groundstate%OEP)) Then
          If (input%groundstate%xctypenumber .Lt. 0) engyx = 0.d0
          If (tlast) Call exxengy
       End If 
! Hybrids
      if (associated(input%groundstate%Hybrid)) then
        if ((input%groundstate%Hybrid%exchangetypenumber==1).and.(ihyb>0)) then
           engyx = engyx + ex_coef*exnl
        end if
      end if
      
!----------------------------!
!     correlation energy     !
!----------------------------!
      engyc = rfinp (1, rhomt, ecmt, rhoir, ecir)
! zero correlation energy for Hartree-Fock
      If ((task .Eq. 5) .Or. (task .Eq. 6)) engyc = 0.d0
! Hybrids
      if (associated(input%groundstate%Hybrid)) then
        if ((input%groundstate%Hybrid%exchangetypenumber==1).and.(ihyb>0)) then
           engyc = ec_coef*engyc
        end if
      end if      
      
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

      if ((task.Eq.5).Or.(task.Eq.6)) then
        !-------------------
        ! Hartree-Fock case
        !-------------------
        ! core electron kinetic energy
        Call energykncr
        engykn = engykncr
        ! kinetic energy from valence states
        Allocate (evecsv(nstsv, nstsv))
        Allocate (c(nstsv, nstsv))
        Do ik = 1, nkpt
          Call getevecsv (vkl(:, ik), evecsv)
          Call zgemm ('N', 'N', nstsv, nstsv, nstsv, zone, &
          &           kinmatc(:,:,ik), nstsv, evecsv, nstsv, zzero, c, nstsv)
          Do ist = 1, nstsv
            zt1 = zdotc (nstsv, evecsv(:, ist), 1, c(:, ist), 1)
            engykn = engykn + wkpt(ik) * occsv(ist,ik) * dble(zt1)
          End Do
        End Do
        Deallocate (evecsv, c)
        
      else if (associated(input%groundstate%Hybrid)) then
        !------------------
        ! HF-based hybrids
        !------------------
        if (input%groundstate%Hybrid%exchangetypenumber == 1) then
          if (ihyb>0) then
            engykn = engykncr
            do ik = 1, nkpt
              do ist = 1, nstfv
                engykn = engykn + wkpt(ik)*occsv(ist,ik)*engyknst(ist,ik)
              end do
            end do
          else
            ! Default way
            engykn =  evalsum - engyvcl - engyvxc - engybxc - engybext - engybmt
            call energykncr
          end if
        Else
          ! OEP-Hybrids: Default way
           engykn =  evalsum - engyvcl - engyvxc - engybxc - engybext - engybmt
        end if
        
      else
        ! Default way
        engykn =  evalsum - engyvcl - engyvxc - engybxc - engybext - engybmt
      end if
!------------------------------!
!     DFT-1/2 contribution     !
!------------------------------!
      if (associated(input%groundstate%dfthalf)) then
        engyhalf = rfinp (1, rhomt, vhalfmt, rhoir, vhalfir)
      endif      
!----------------------!
!     total energy     !
!----------------------!
      engytot = engykn + 0.5d0 * engyvcl + engymad + engyx + engyc + engycbc
!dispersion correction
      If ( tlast .And. input%groundstate%vdWcorrection .Ne. "none" ) Then
         If ( input%groundstate%vdWcorrection .Eq. "DFTD2" ) Then
            Call DFT_D2_energy
         Else If ( input%groundstate%vdWcorrection .Eq. "TSvdW" ) Then 
            Call TS_vdW_energy
         End If
         engytot = engytot + e_disp
      End If
! dipole correction      
      if ((iscl>0).and.(input%groundstate%dipolecorrection)) engytot = engytot+0.5*endipc      
! add the LDA+U correction if required
      If (ldapu .Ne. 0) engytot = engytot + engylu
! WRITE(*,*) "end energy"
      Return
End Subroutine
!EOC
