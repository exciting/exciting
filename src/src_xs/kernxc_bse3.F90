!
!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kernxc_bse3
! !INTERFACE:
!
!
Subroutine kernxc_bse3
! !USES:
      Use modinput
      Use modmain
      Use modxs
      Use m_getevalsvr
      Use m_getpemat
      Use m_xsgauntgen
      Use m_findgntn0
      Use m_genwgrid
      Use m_xszoutpr3
      Use m_genfilname
      Use m_getunit
      Use m_xszoutpr3
! !INPUT/OUTPUT PARAMETERS:
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
!
! !REVISION HISTORY:
!   Created March 2009 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'kernxc_bse3'
      Real (8), Parameter :: eps = 1.d-5
      Integer, Parameter :: iqmt = 1, noptc = 3
      Character (256) :: filnam
      Integer :: iknr, jknr, iv, ic, jv, jc, si, sj, nv, nc, wsiz, n, &
     & iop, iw, un, recl
      Integer :: nvl, ncu
      Complex (8), Allocatable :: w (:), mat (:, :), wmat (:, :), wm &
     & (:, :, :, :)
      Complex (8), Allocatable :: l0mat (:, :), l0mata (:, :)
      Complex (8), Allocatable :: hmat (:, :), hmat2 (:, :)
      Complex (8), Allocatable :: fxc (:, :, :)
      Complex (8), Allocatable :: xiout (:, :, :), pmout (:, :, :)
      Complex (8), Allocatable :: xiuot (:, :, :), pmuot (:, :, :)
      Complex (8), Allocatable :: me (:, :), mea (:, :)
      Real (8), Allocatable :: ev (:), de (:), scisk (:, :)
      Integer, Allocatable :: widx (:, :, :)
      Integer, External :: idxkkp
!
      Write (*,*) 'initializing...'
!
      input%xs%emattype = 2
      Call init0
      Call init1
      Call init2
!
      Write (*,*) 'preparing...'
!
      Call readfermi
  ! initialize states below and above the Fermi energy
      Call initocc (nbfbse, nafbse)
      Call xssave0
      Call xsgauntgen (Max(input%groundstate%lmaxapw, lolmax), &
     & input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
      Call findgntn0 (Max(input%xs%lmaxapwwf, lolmax), &
     & Max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      call findocclims(0, ikmapikq(:,1), istocc0, istunocc0, isto0, isto, istu0, istu)
      istunocc = istunocc0
      istocc = istocc0
      Call ematbdcmbs (input%xs%emattype)
      Allocate (w(nwdf))
      Call genwgrid (nwdf, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_cmplx=w)
!
      Write (*,*) 'done.'
!
      n = ngq (iqmt)
      nv = nst1
      nc = nst3
      wsiz = nv * nc * nkptnr
  ! same window for bands as bse-routine
      nvl = nv - nbfbse + 1
      ncu = nafbse
!
      Write (*,*) 'n, nv, nc, nvl, ncu, wsiz', n, nv, nc, nvl, ncu, &
     & wsiz
!
      Allocate (wm(nv, nc, nv, nc), widx(nv, nc, nkptnr))
      Allocate (de(wsiz), wmat(wsiz, wsiz), mat(wsiz, wsiz))
      Allocate (l0mat(wsiz, wsiz), me(-3:n, wsiz))
      Allocate (l0mata(wsiz, wsiz), mea(-3:n, wsiz))
      Allocate (hmat(-3:n, wsiz), hmat2(wsiz,-3:n))
      Allocate (ev(nstsv))
      Allocate (fxc(-3:n,-3:n, nwdf))
      Allocate (scisk(nst1, nst3))
      Allocate (xiout(nv, nc, n), pmout(3, nv, nc))
      Allocate (xiuot(nc, nv, n), pmuot(3, nc, nv))
      If (allocated(pmou)) deallocate (pmou)
      Allocate (pmou(3, nv, nc))
      If (allocated(pmuo)) deallocate (pmuo)
      Allocate (pmuo(3, nc, nv))
      If (allocated(deou)) deallocate (deou)
      Allocate (deou(nst1, nst3))
      If (allocated(deuo)) deallocate (deuo)
      Allocate (deuo(nst3, nst1))
      If (allocated(docc12)) deallocate (docc12)
      Allocate (docc12(nst1, nst3))
      If (allocated(docc21)) deallocate (docc21)
      Allocate (docc21(nst3, nst1))
!
!
  ! set up indices
      si = 0
      Do iknr = 1, nkptnr
         Do iv = 1, nv
            Do ic = 1, nc
               si = si + 1
               widx (iv, ic, iknr) = si
            End Do
         End Do
      End Do
!
!
!
  ! set up energies and their differences
      Do iknr = 1, nkptnr
         Call getevalsvr ('EVALSV_SCR.OUT', 1, nstsv, vkl(:, iknr), ev)
         Do iv = 1, nv
            Do ic = 1, nc
               si = widx (iv, ic, iknr)
               de (si) = ev (istocc+ic) - ev (iv)
            End Do
         End Do
      End Do
!
!
      Write (*,*) 'calculating matrix elements....'
!
  ! calculate matrix elements
      input%xs%emattype = 1
      Call ematbdcmbs (input%xs%emattype)
      Call ematrad (iqmt)
      Call ematqalloc
      Do iknr = 1, nkptnr
         Write (unitout,*) 'matrix elements iknr=', iknr
         Call flushifc (unitout)
         Call ematqk1 (iqmt, iknr)
         Call getdevaldoccsv (iqmt, iknr, iknr, istl1, istu1, istl2, &
        & istu2, deou, docc12, scisk)
         Call getpemat (iqmt, iknr, 'PMAT_SCR.OUT', '', m12=xiout, &
        & p12=pmout)
         Do iv = 1, nv
            Do ic = 1, nc
               si = widx (iv, ic, iknr)
               Do iop = 1, noptc
                  me (-iop, si) = pmout (iop, iv, ic)
                  mea (-iop, si) = pmuot (iop, ic, iv)
               End Do
               me (1:, si) = xiout (iv, ic, :)
               mea (1:, si) = xiuot (ic, iv, :)
            End Do
         End Do
      End Do
      input%xs%emattype = 2
      Call ematbdcmbs (input%xs%emattype)
!
!
      Write (*,*) 'done.'
!
      If ((input%xs%tddft%fxctypenumber .Eq. 7) .Or. &
     & (input%xs%tddft%fxctypenumber .Eq. 8)) Then
         Call getbsediag
         Write (unitout, '("Info(", a, "): read diagonal of BSE kernel"&
        &)') trim (thisnam)
         Write (unitout, '(" mean value : ", 2g18.10)') bsed
      End If
!
  ! set up W-matrix and L^0 matrix
      wmat (:, :) = zzero
      Do iknr = 1, nkptnr
         Do jknr = iknr, nkptnr
            Call getbsemat ('SCCLI.OUT', idxkkp(iknr, jknr, nkptnr), &
           & nv, nc, wm)
            Do iv = nvl, nv
               Do ic = 1, ncu
                  si = widx (iv, ic, iknr)
                  Do jv = nvl, nv
                     Do jc = 1, ncu
                        sj = widx (jv, jc, jknr)
                        If (si .Ne. sj) Then
                           wmat (si, sj) = - wm (iv, ic, jv, jc)
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
!
      Do si = 1, wsiz
         Do sj = si + 1, wsiz
            wmat (sj, si) = conjg (wmat(si, sj))
         End Do
      End Do
!
  ! loop over w-points
      Do iw = 1, nwdf
!
         Write (*,*) 'w-loop, at: ', iw
!
    ! set up L0-matrix
         l0mat (:, :) = zzero
         l0mata (:, :) = zzero
         Do si = 1, wsiz
            l0mat (si, si) = 1.d0 / &
           & (w(iw)+bsed-de(si)+zi*input%xs%broad)
            l0mata (si, si) = - 1.d0 / &
           & (w(iw)+bsed+de(si)+zi*input%xs%broad)
         End Do
!
!    ! do the 4 matrix multiplications
!    hmat=matmul(me,l0mat)
!    hmat=matmul(hmat,wmat)
!    hmat=matmul(hmat,l0mat)
!    fxc(:,:,iw)=2.d0*matmul(hmat,conjg(transpose(me)))
!
    ! resonant contribution
         Do si = 1, wsiz
            hmat (:, si) = me (:, si) * l0mat (si, si)
         End Do
         hmat = matmul (hmat, wmat)
         Do si = 1, wsiz
            hmat2 (si, :) = l0mat (si, si) * conjg (me(:, si))
         End Do
         fxc (:, :, iw) = 2.d0 * matmul (hmat, hmat2) / (nkptnr*omega)
!
         If (input%xs%tddft%aresfxc) Then
    ! anti-resonant contribution
            Do si = 1, wsiz
               hmat (:, si) = mea (:, si) * l0mata (si, si)
            End Do
            hmat = matmul (hmat,-conjg(wmat))
            Do si = 1, wsiz
               hmat2 (si, :) = l0mata (si, si) * conjg (mea(:, si))
            End Do
            fxc (:, :, iw) = fxc (:, :, iw) + 2.d0 * matmul (hmat, &
           & hmat2) / (nkptnr*omega)
         End If
!
      End Do
!
!
   ! deallocate the wmat arrays
      Deallocate (wmat, mat)
!
      Write (*,*) 'writing out kernel...'
!
!
  ! write out kernel
      Inquire (IoLength=Recl) n, fxc (-3:-1,-3:-1, 1), fxc (-3:-1, 1:, &
     & 1), fxc (1:,-3:-1, 1), fxc (1:, 1:, 1)
      Call genfilname (basename='FXC_BSE', asc=.False., bzsampl=0, &
     & acont=input%xs%tddft%acont, nar= .Not. input%xs%tddft%aresfxc, &
     & tord=input%xs%tddft%tordfxc, iqmt=iqmt, filnam=filnam)
      Call getunit (un)
      Open (un, File=trim(filnam), Form='unformatted', Action='write', &
     & Status='replace', Access='direct', Recl=Recl)
      Do iw = 1, nwdf
         Write (un, Rec=iw) n, fxc (-3:-1,-3:-1, iw), fxc (-3:-1, 1:, &
        & iw), fxc (1:,-3:-1, iw), fxc (1:, 1:, iw)
         Write (8888, '(i6, 6g18.10)') iw, fxc (-1,-1, iw), fxc (-2,-2, &
        & iw), fxc (-3,-3, iw)
      End Do
      Close (un)
!
      Deallocate (wm, widx)
      Deallocate (de, me)
      Deallocate (ev, hmat, hmat2)
      Deallocate (xiout, pmout, xiuot, pmuot, scisk, l0mat, l0mata)
!
End Subroutine kernxc_bse3
!EOC
