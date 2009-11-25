!
!
! Copyright (C) 2006-2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine epsconv
      Use modmain
      Use modinput
      Use modxs
      Use modtetra
      Use m_getunit
      Use m_genfilname
      Implicit None
  ! local variables
      Logical, Parameter :: tsmooth = .False.
      Character (*), Parameter :: thisnam = 'epsconv'
      Character (256) :: filnam
      Integer :: iq, iw, iwp, m, n, oct1, oct2, nc, un
      Logical :: exis, tq0
      Real (8), Parameter :: epsc = 1.d-8
      Real (8), Allocatable :: w (:), epst (:, :), lor (:), f (:), f1 &
     & (:), g (:), g1 (:), cf (:, :)
      Complex (8), Allocatable :: eps (:)
      Integer, External :: l2int
      Logical, External :: tqgamma
  ! initialize universal variables
      Call init0
      Call init1
      Call init2
  ! original sampling method fo Brillouine zone
      bzsampl = l2int (input%xs%tetra%tetradf)
      Allocate (w(nwdf), epst(nwdf, 2), eps(nwdf), lor(nwdf))
      Allocate (f(nwdf), f1(nwdf), g(nwdf), g1(nwdf), cf(3, nwdf))
  ! loop over q-points
      Do iq = 1, nqpt
     ! matrix size for local field effects
         n = ngq (iq)
         tq0 = tqgamma (iq)
         nc = 1
         If (tq0) nc = 3
     ! neglect/include local field effects
         Do m = 1, n, Max (n-1, 1)
        ! loop over longitudinal components for optics
            Do oct1 = 1, nc
               Do oct2 = 1, nc
           ! generate filename for Tetrahedron method
                  Call genfilname (basename='EPSILON', bzsampl=bzsampl, &
                 & nar= .Not. input%xs%tddft%aresdf, nlf=(m == 1), &
                 & fxctype=input%xs%tddft%fxctypenumber, tq0=tq0, &
                 & oc1=oct1, oc2=oct2, iqmt=iq, filnam=filnam)
           ! check for file to read
                  Inquire (File=trim(filnam), Exist=exis)
                  If ( .Not. exis) Then
                     Write (*,*) 'Error(' // trim (thisnam) // '): file&
                    & does not exist: ' // trim (filnam)
                     Call terminate
                  End If
                  Call getunit (un)
                  Open (Unit=un, File=trim(filnam), Form='formatted', &
                 & Action='read', Status='old')
           ! read energies and Re and Im of eps
                  Do iw = 1, nwdf
                     Read (un,*) w (iw), epst (iw, 1), epst (iw, 2)
                  End Do
                  Close (un)
           ! generate filename for output (_s_ymmetric _l_orentzian)
                  Call genfilname (basename='EPSILON_sl', &
                 & bzsampl=bzsampl, nar= .Not. input%xs%tddft%aresdf, &
                 & nlf=(m == 1), fxctype=input%xs%tddft%fxctypenumber, &
                 & tq0=tq0, oc1=oct1, oc2=oct2, iqmt=iq, filnam=filnam)
                  Open (Unit=un, File=trim(filnam), Form='formatted', &
                 & Action='write', Status='replace')
           ! rescale energies read from file to atomic units
                  w (:) = w (:) / escale
           ! complex dielectric function
                  eps (:) = epst (:, 1) + zi * epst (:, 2)
           ! convolution with Lorentzian
                  Do iw = 1, nwdf
                     Do iwp = 1, nwdf
                 ! standard Lorentzian with peak at w(iw)
                 !lor(iwp)=(1/pi)*broad/((w(iw)-w(iwp))**2+broad**2)
                 ! antisymmetric Lorentzian at w(iw) and -w(iw)
                 ! with norm arctan(w/broad) to assure zero crossing
                        lor (iwp) = &
                       & (1.d0/(2.d0*Atan(w(iw)/input%xs%broad))) * &
                       & (input%xs%broad/((w(iw)-w(iwp))**2+&
                       & input%xs%broad**2)-input%xs%broad/((-w(iw)-&
                       & w(iwp))**2+input%xs%broad**2))
                        If (w(iw) < epsc) lor (iwp) = 0.d0
                        f (iwp) = lor (iwp) * aimag (eps(iwp))
                        f1 (iwp) = lor (iwp) * dble (eps(iwp))
                     End Do
              ! do the convolution
                     Call fderiv (-1, nwdf, w, f, g, cf)
                     Call fderiv (-1, nwdf, w, f1, g1, cf)
                     Write (un, '(4g18.10)') w (iw) * escale, g1 &
                    & (nwdf), g (nwdf), (pi/input%xs%broad) * &
                    & input%xs%broad ** 2 / &
                    & ((w(iw)-w(iwp))**2+input%xs%broad**2)
                  End Do ! iw
                  If (tsmooth) Then
                     Call fsmooth (input%xs%doswindow%nsmdos, nwdf, 1, &
                    & g)
                     Call fsmooth (input%xs%doswindow%nsmdos, nwdf, 1, &
                    & g1)
                  End If
                  Close (un)
               End Do ! oct
            End Do
         End Do
      End Do
      Deallocate (w, epst, eps, lor, f, f1, g, g1, cf)
End Subroutine epsconv
