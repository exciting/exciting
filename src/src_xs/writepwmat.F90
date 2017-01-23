!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writepwmat
! !INTERFACE:
!
!
Subroutine writepwmat
! !USES:
      Use modmain
      Use modxs
      Use m_genfilname
! !DESCRIPTION:
!   Calculates the matrix elements of the plane wave
!   $e^{-i({\bf G}+{\bf q}){\bf r}}$
!   using routine {\tt genpwmat} and writes them to
!   direct access file {\tt PWMAT.OUT}.
!
! !REVISION HISTORY:
!   Created November 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Integer :: iq
      Integer :: ik, ikp, recl, isymkp, igq, un, un2
      Real (8) :: vpl (3), vkpl (3)
      Complex (8), Allocatable :: apwalmk (:, :, :, :), apwalmkp (:, :, &
     & :, :)
      Complex (8), Allocatable :: evecfvk (:, :), evecfvkp (:, :)
      Complex (8), Allocatable :: evecsvk (:, :), evecsvkp (:, :)
      Complex (8), Allocatable :: pwmat (:, :, :), pwmatf (:, :, :)
      Integer :: j, iknr, lspl, isym, jsym, s (3, 3), vg (3), ig, ist, &
     & jst, reclfull
      Real (8) :: c (3, 3), t1, t2, t3
      Complex (8) :: zt1, zt2
      Complex (8), Allocatable :: yiou (:, :, :), yiuo (:, :, :)
      Real (8), Parameter :: epsrot = 1.d-12
      Character (256) :: fname
      Character(256) :: string
  ! initialise universal variables
      Call init0
      Call init1
      Call readfermi
      Call init2
      call findocclims(0, ikmapikq(:,1), istocc0, istunocc0, isto0, isto, istu0, istu)
      istunocc = istunocc0
      istocc = istocc0
!
      Do iq = 1, nqpt
         Allocate (apwalmk(ngkmax, apwordmax, lmmaxapw, natmtot))
         Allocate (apwalmkp(ngkmax, apwordmax, lmmaxapw, natmtot))
         Allocate (evecfvk(nmatmax, nstfv))
         Allocate (evecfvkp(nmatmax, nstfv))
         Allocate (evecsvk(nstsv, nstsv))
         Allocate (evecsvkp(nstsv, nstsv))
     ! allocate the momentum matrix elements array
         Allocate (pwmat(ngq(iq), nstsv, nstsv))
         Allocate (pwmatf(ngq(iq), nstsv, nstsv))
!
     ! allocate matrix elements array
         If (allocated(xiou)) deallocate (xiou)
         Allocate (xiou(nstocc0, nstunocc0, ngq(iq)))
         If (allocated(xiuo)) deallocate (xiuo)
         Allocate (xiuo(nstunocc0, nstocc0, ngq(iq)))
!
         If (allocated(yiou)) deallocate (yiou)
         Allocate (yiou(nstocc0, nstunocc0, ngq(iq)))
         If (allocated(yiuo)) deallocate (yiuo)
         Allocate (yiuo(nstunocc0, nstocc0, ngq(iq)))
!
! read density and potentials from file
        If (hybridhf) Then
! in case of HF hybrids use PBE potential
            string=filext
            filext='_PBE.OUT'
            Call readstate
            filext=string
        Else
           Call readstate
        End If
     ! find the new linearisation energies
         Call linengy
     ! generate the APW radial functions
         Call genapwfr
     ! generate the local-orbital radial functions
         Call genlofr
! update potential in case if HF Hybrids
        If (hybridhf) Call readstate
     ! find the record length
         Inquire (IoLength=Recl) pwmat
         Call genfilname (basename='PWMAT', iqmt=iq, asc=.True., &
        & filnam=fname)
         Open (50, File=trim(fname), Action='WRITE', Form='FORMATTED', &
        & Status='REPLACE')
         Do ik = 1, nkpt
            Write (*,*) 'Info(writepwmat): ik', ik
        ! get the eigenvectors from file for k-point
            Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfvk)
            Call getevecsv (vkl(1, ik), evecsvk)
        ! find the matching coefficients for k-point
            Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
           & sfacgk(1, 1, 1, ik), apwalmk)
        ! find k-point equivalent to k+q
            vpl (:) = vql (:, iq) + vkl (:, ik)
            Call findkpt (vpl, isymkp, ikp)
            Write (*,*) 'ik,ikp', ik, ikp
            vkpl (:) = vkl (:, ikp)
        ! get the eigenvectors from file for kp-point
            Call getevecfv (vkl(1, ikp), vgkl(1, 1, 1, ikp), evecfvkp)
            Call getevecsv (vkl(1, ikp), evecsvkp)
        ! find the matching coefficients for kp-point
            Call match (ngk(1, ikp), gkc(1, 1, ikp), tpgkc(1, 1, 1, &
           & ikp), sfacgk(1, 1, 1, ikp), apwalmkp)
        ! calculate the matrix elements of the plane wave
            Call genpwmat (vql(1, iq), ngqmax, ngq(iq), vgqc(1, 1, iq), &
           & gqc(1, iq), igqig(1, iq), ylmgq(1, 1, iq), sfacgq(1, 1, &
           & iq), vkl(1, ik), ngk(1, ik), igkig(1, 1, ik), apwalmk, &
           & evecfvk, evecsvk, vkl(1, ikp), ngk(1, ikp), igkig(1, 1, &
           & ikp), apwalmkp, evecfvkp, evecsvkp, pwmat)
        ! write to ASCII file
            Do igq = 1, ngq (iq)
               Do ist = 1, nstsv
                  Do jst = 1, nstsv
                     Write (50, '(4i5,3g18.10)') ik, igq, ist, jst, &
                    & pwmat (igq, ist, jst), Abs (pwmat(igq, ist, jst)) &
                    & ** 2
                  End Do
               End Do
            End Do
            Do igq = 1, ngq (iq)
               xiou (:, :, igq) = pwmat (igq, 1:nstocc0, &
              & nstocc0+1:nstsv)
               xiuo (:, :, igq) = pwmat (igq, nstocc0+1:nstsv, &
              & 1:nstocc0)
            End Do
        ! write to direct access file
            Inquire (IoLength=Recl) nstocc0, nstunocc0, nkpt, ngq (iq), &
           & vql (:, iq), vkl (:, ik), xiou, xiuo
            Inquire (IoLength=reclfull) pwmat
            un = 51
            un2 = 52
!
            Call genfilname (basename='EMAT', iqmt=iq, filnam=fnemat)
            Open (Unit=un, File=trim(fnemat), Form='unformatted', &
           & Action='write', Access='direct', Recl=Recl)
            Write (un, Rec=ik) nstocc0, nstunocc0, nkpt, ngq (iq), vql &
           & (:, iq), vkl (:, ik), xiou, xiuo
            Close (un)
!
        ! rotate matrix element if k-point set is reduced
            If (nkpt .Ne. nkptnr) Then
               Do j = 1, nsymcrysstr (ik)
                  iknr = ikstrmapiknr (j, ik)
                  isym = scmapstr (j, ik)
                  jsym = scimap (isym)
                  lspl = lsplsymc (jsym)
              ! rotation in Cartesian coordinates
                  c (:, :) = symlatc (:, :, lspl)
              ! rotation in lattice coordinates
                  s (:, :) = symlat (:, :, lspl)
                  Do igq = 1, ngq (iq)
                 ! let c=a^(-1) and note that in exciting the symmetry operations
                 ! are defined as y = as a (x + T_a)
                 ! (G+q).T_c
                     t1 = twopi * dot_product (vgql(:, igq, iq), &
                    & matmul(s, vtlsymc(:, jsym)))
                 ! exp(-i(G+q)T_c)
                     t2 = Cos (t1)
                     t3 = - Sin (t1)
                     If (Abs(t2) .Lt. epsrot) t2 = 0.d0
                     If (Abs(t3) .Lt. epsrot) t3 = 0.d0
                     zt1 = cmplx (t2, t3, 8)
                 ! G-vector of G+q
                     ig = igqig (igq, iq)
                     vg (:) = ivg (:, ig)
                 ! G*c
                     vg = matmul (vg, s)
                 ! index of G+q-vector of new G-vector
                     ig = ivgigq (vg(1), vg(2), vg(3), iq)
!
                     Write (80, '(4i6,3x,3i5,3x,i6)') iknr, ik, isym, &
                    & igq, vg, ig
!
                 ! write to ASCII file
                     Do ist = 1, nstsv
                        Do jst = 1, nstsv
                           zt2 = zt1 * pwmat (ig, ist, jst)
                       ! matrix elements for full k-point set
                           pwmatf (igq, ist, jst) = zt2
                           Write (90, '(7i5,5g18.10)') iknr, ik, isym, &
                          & igq, ist, jst, ig, zt2, Abs (zt2) ** 2, zt1
!
                           If ((ist .Le. nstocc0) .And. (jst .Gt. &
                          & nstocc0)) Then
                              yiou (ist, jst-nstocc0, igq) = zt2
                              Write (300+iq, '(a,4i6,3g18.10)') 'ik,igq&
                             &,i1,i2', iknr, igq, ist, jst - nstocc0, &
                             & zt2, Abs (zt2) ** 2
                           End If
                           If ((ist .Gt. nstocc0) .And. (jst .Le. &
                          & nstocc0)) Then
                              yiuo (ist-nstocc0, jst, igq) = zt2
                              Write (400+iq, '(a,4i6,3g18.10)') 'ik,igq&
                             &,i1,i2', iknr, igq, ist - nstocc0, jst, &
                             & zt2, Abs (zt2) ** 2
                           End If
                        End Do
                     End Do
                 ! end loop over G+q vectors
                  End Do
!
              ! write v-c and c-v matrix elements
                  Open (Unit=un, File='EMAT_NR_Q00001.OUT', Form='unfor&
                 &matted', Action='write', Access='direct', Recl=Recl)
                  Write (un, Rec=iknr) nstocc0, nstunocc0, nkptnr, ngq &
                 & (iq), vql (:, iq), vklnr (:, iknr), yiou, yiuo
                  Close (un)
!
              ! write full matrix elements
                  Open (Unit=un, File='EMAT_FULL_NR_Q00001.OUT', Form='&
                 &unformatted', Action='write', Access='direct', &
                 & Recl=reclfull)
                  Write (un, Rec=iknr) pwmatf
                  Close (un)
!
              ! end loop over elements of star
               End Do
            End If
!
        ! end loop over k
         End Do
         Close (50)
!
     ! read and write matrix elements of non-reduced k-points to ASCII file
         If (nkpt .Ne. nkptnr) Then
            Open (Unit=un, File='EMAT_FULL_NR_Q00001.OUT', Form='unform&
           &atted', Action='read', Access='direct', Recl=reclfull)
            Open (Unit=un2, File='PWMAT_NR_ASC.OUT', Form='formatted', &
           & Action='write', Status='replace')
            Do iknr = 1, nkptnr
               Read (un, Rec=iknr) pwmat
               Do igq = 1, ngq (iq)
                  Do ist = 1, nstsv
                     Do jst = 1, nstsv
                        Write (un2, '(4i5,3g18.10)') iknr, igq, ist, &
                       & jst, pwmat (igq, ist, jst), Abs (pwmat(igq, &
                       & ist, jst)) ** 2
                     End Do
                  End Do
               End Do
            End Do
            Close (un2)
            Close (un)
         End If
         Deallocate (apwalmk, evecfvk, evecsvk, apwalmkp, evecfvkp, &
        & evecsvkp, pwmat, pwmatf)
     ! end loop over q-points
      End Do
      Write (*,*)
      Write (*, '("Info(writepwmat):")')
      Write (*, '(" matrix elements of the plane wave written to file P&
     &WMAT.OUT")')
      Write (*,*)
End Subroutine writepwmat
!EOC
