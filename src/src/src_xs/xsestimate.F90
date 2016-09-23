!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xsestimate
      Use modmain
      Use modinput
      Use modxs
      Use m_genfilname
      Use m_getunit
      Use m_gndstateq
      Implicit None
      Real (8) :: d_ev (2), d_apwcmt (2), d_pmat (2), d_emat (2), d_x0 &
     & (2), d_tetw (2)
      Real (8) :: d_fxcbse, d_sci, d_eci, d_tot (2)
      Real (8) :: m_bseham, m_bseeigvec
      Real (8) :: mb, gb
      Integer :: sreal, scmplx, nkkpt, bsehamsiz, un
      Call init0
      Call init1
      Call init2
      nst2 = input%groundstate%nempty + 1
      nst1 = nstsv - nst2
  ! estimate disk space usage
      sreal = 8
      scmplx = 16
      mb = 1024 ** 2
      gb = 1024 * mb
      nkkpt = (nkpt*(nkpt+1)) / 2
      bsehamsiz = nkpt * nbfbse * nafbse
  ! eigenvectors (spin-unpolarized)
      d_ev (1) = dble (scmplx*nkpt) * dble (nmatmax*nstfv)
      d_ev (2) = d_ev (1) * dble (nqpt+1)
  ! tetrahedron weights
      d_tetw (1) = dble (sreal*nkpt*nwdf) * dble (nst1*nst2*3)
      d_tetw (2) = d_tetw (1)
  ! expansion coefficients of muffin-tin wavefunctions
      d_apwcmt (1) = dble (scmplx*nkpt) * dble &
     & (nstfv*apwordmax*lmmaxapw*natmtot)
      d_apwcmt (2) = d_apwcmt (1) * dble (nqpt+1)
  ! matrix elements of the momentum operator
      d_pmat (1) = dble (scmplx*nkpt) * dble (3*nstsv*nstsv)
      d_pmat (2) = d_pmat (1)
  ! matrix elements of the plane wave
      d_emat (1) = dble (scmplx*nkpt) * dble (2*nst1*nst2*ngqmax)
      d_emat (2) = d_emat (1) * dble (nqpt)
  ! Kohn-Sham response function
      d_x0 (1) = dble (scmplx*nwdf) * dble (ngqmax**2+3*2*ngqmax+3)
      d_x0 (2) = d_x0 (1) * dble (nqpt)
  ! BSE fxc-kernel
      d_fxcbse = d_x0 (1)
  ! screened and exchange Coulomb interaction
      d_sci = scmplx * nkkpt * (nst1*nst2) ** 2
      d_eci = d_sci
  ! totals
      d_tot (1) = d_ev (1) + d_tetw (1) + d_apwcmt (1) + d_pmat (1) + &
     & d_emat (1) + d_x0 (1) + d_sci + d_eci + d_fxcbse
      d_tot (2) = d_ev (2) + d_tetw (2) + d_apwcmt (2) + d_pmat (2) + &
     & d_emat (2) + d_x0 (2) + d_sci + d_eci + d_fxcbse
  ! memory consumption
      m_bseham = scmplx * bsehamsiz ** 2
      m_bseeigvec = scmplx * bsehamsiz * input%xs%BSE%nexcitmax
  !write information to file
      Call getunit (un)
      Open (un, File='SIZES.OUT', Form='formatted', Action='write', &
     & Status='replace')
      Write (un,*)
      Write (un, '(a)') 'Relevant parameters:'
      Write (un, '(a, i6)') ' nwdf	   :', nwdf
      Write (un, '(a, i6)') ' nqpt	   :', nqpt
      Write (un, '(a, i6)') ' ngqmax    :', ngqmax
      Write (un, '(a, i6)') ' nkpt	   :', nkpt
      Write (un, '(a, i6)') ' nkkpt	   :', nkkpt
      Write (un, '(a, i6)') ' nstsv	   :', nstsv
      Write (un, '(a, i6)') ' nst1	   :', nst1
      Write (un, '(a, i6)') ' nst2	   :', nst2
      Write (un, '(a, i6)') ' nbfbse    :', nbfbse
      Write (un, '(a, i6)') ' nafbse    :', nafbse
      Write (un, '(a, i6)') ' nmatmax   :', nmatmax
      Write (un, '(a, i6)') ' ngkmax    :', ngkmax
      Write (un, '(a, i6)') ' nlotot    :', nlotot
      Write (un, '(a, i6)') ' apwordmax :', apwordmax
      Write (un, '(a, i6)') ' lmmaxapw  :', lmmaxapw
      Write (un, '(a, i6)') ' natmtot   :', natmtot
      Write (un,*)
      Write (un, '(a)') 'Other parameters:'
      Write (un, '(a, i8)') ' associated(input%xs%BSE) Hamiltonian size&
     &  :', bsehamsiz
      Write (un,*)
      Write (un, '(a)') 'Estimated memory consumption in GB (one q-poin&
     &t/total):'
      Write (un, '(a, 2f12.3)') ' associated(input%xs%BSE) Hamiltonian	&
     &  :', m_bseham / gb
      Write (un, '(a, 2f12.3)') ' associated(input%xs%BSE) eigenvectors&
     &	  :', m_bseeigvec / gb
      Write (un, '(a, 2f12.3)') ' associated(input%xs%BSE) totals *		  &
     &:', (m_bseham+m_bseeigvec) / gb
      Write (un,*)
      Write (un, '(a)') 'Estimated disk usage for large files in GB (on&
     &e q - point/total):'
      Write (un, '(a, 2f12.3)') ' eigenvectors		  :', d_ev (1) / gb, &
     & d_ev (2) / gb
      Write (un, '(a, 2f12.3)') ' tetrahedron weights	  :', d_tetw (1) &
     & / gb, d_tetw (2) / gb
      Write (un, '(a, 2f12.3)') ' MT expansion coeffs	  :', d_apwcmt &
     & (1) / gb, d_apwcmt (2) / gb
      Write (un, '(a, 2f12.3)') ' matr. el. of mom. op.	  :', d_pmat &
     & (1) / gb, d_pmat (2) / gb
      Write (un, '(a, 2f12.3)') ' matr. el. of plane wave	  :', d_emat &
     & (1) / gb, d_emat (2) / gb
      Write (un, '(a, 2f12.3)') ' KS response function	  :', d_x0 (1) / &
     & gb, d_x0 (2) / gb
      Write (un, '(a, 2f12.3)') ' associated(input%xs%BSE)-kernel		    &
     &   :', d_fxcbse / gb
      Write (un, '(a, 2f12.3)') ' screened Coulomb interaction :', &
     & d_sci / gb
      Write (un, '(a, 2f12.3)') ' exchange Coulomb interaction :', &
     & d_eci / gb
      Write (un, '(a, 2f12.3)') ' - - - - - - - - - - - - - - - - - - -&
     & - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - -&
     & -  - - - - - -- - - - - - '
      Write (un, '(a, 2f12.3)') ' total			  :', d_tot (1) / gb, d_tot &
     & (2) / gb
      Write (un,*)
      Close (un)
End Subroutine xsestimate
