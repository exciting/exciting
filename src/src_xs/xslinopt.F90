!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xslinopt (iq)
      Use modmain
      Use modinput
      Use modxs
      Use modtetra
      Use modmpi
      Use m_genwgrid
      Use m_pade
      Use m_genloss
      Use m_gensigma
      Use m_genkerr
      Use m_gensumrls
      Use m_writeeps
      Use m_writeloss
      Use m_writesigma
      Use m_writekerr
      Use m_writesumrls
      Use m_getunit
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Character (*), Parameter :: thisnam = 'xslinopt'
      Character (256) :: filnam
      Complex (8), Allocatable :: mdf (:), mdf1 (:), mdf2 (:, :, :), w &
     & (:), wr (:), sigma (:), kerr(:)
      Real (8), Allocatable :: wplot (:), loss (:, :, :)
      Real (8), Allocatable :: eps1 (:), eps2 (:), cf (:, :)
      Real (8) :: sumrls (3), brd
      Integer :: n, m, recl, iw, nc, oct1, oct2, octl, &
     & octu, optcompt (2)
      Logical :: tq0
      Logical, External :: tqgamma
      tq0 = tqgamma (iq)
  ! number of components (3 for q=0)
      nc = 1
      If (tq0) nc = 3
  ! matrix size for local field effects
      n = ngq (iq)
      Allocate (mdf1(nwdf), mdf2(3, 3, input%xs%energywindow%points), w(nwdf), &
     & wr(input%xs%energywindow%points), wplot(input%xs%energywindow%points), &
     & mdf(input%xs%energywindow%points), loss(3, 3, input%xs%energywindow%points), &
     & sigma(input%xs%energywindow%points), cf(3, input%xs%energywindow%points),&
     & kerr(input%xs%energywindow%points))
      Allocate (eps1(input%xs%energywindow%points), &
     & eps2(input%xs%energywindow%points))
      mdf2 (:, :, :) = zzero
  ! generate energy grids
      brd = 0.d0
      If (input%xs%tddft%acont) brd = input%xs%broad
      Call genwgrid (nwdf, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_cmplx=w)
      Call genwgrid (input%xs%energywindow%points, &
     & input%xs%energywindow%intv, .False., brd, w_cmplx=wr)
      wplot = dble (wr)
  ! record length
      Inquire (IoLength=Recl) mdf1 (1)
      Call getunit (unit1)
  ! neglect/include local field effects
      Do m = 1, n, Max (n-1, 1)
     ! loop over longitudinal components for optics
         Do oct1 = 1, nc
            If (input%xs%dfoffdiag) Then
               octl = 1
               octu = nc
            Else
               octl = oct1
               octu = oct1
            End If
            Do oct2 = octl, octu
           ! file name for inverse of dielectric function
               Call genfilname (basename='IDF', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctype=input%xs%tddft%fxctypenumber, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=filnam)
           ! read macroscopic dielectric function (original frequencies)
               Open (unit1, File=trim(filnam), Form='unformatted', &
              & Action='read', Status='old', Access='direct', &
              & Recl=Recl)
               Do iw = 1, nwdf
                  Read (unit1, Rec=iw) mdf1 (iw)
               End Do
               Close (unit1)
           ! analytic continuation
               If (input%xs%tddft%acont) Then
                  Call pade (input%xs%energywindow%points, wr, nwdf, w, &
                 & mdf1, mdf)
               Else
                  mdf (:) = mdf1 (:)
               End If
               mdf2 (oct1, oct2, :) = mdf (:)
            End Do
         End Do
! STK
         call genloss (mdf2, loss, nc)
!
         Do oct1 = 1, nc
            If (input%xs%dfoffdiag) Then
               octl = 1
               octu = nc
            Else
               octl = oct1
               octu = oct1
            End If
            Do oct2 = octl, octu
               optcompt (:) = (/ oct1, oct2 /)
           ! symmetrize the macroscopic dielectric function tensor
               if (tq0) Call symt2app (oct1, oct2, nwdf, symt2, mdf2, mdf)
           ! file names for spectra
               Call genfilname (basename='EPSILON', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fneps)
               Call genfilname (basename='LOSS', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fnloss)
               Call genfilname (basename='SIGMA', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fnsigma)
               if (tq0) Call genfilname (basename='SUMRULES', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fnsumrules)
           ! generate optical functions
! STK
!              Call genloss (mdf, loss)
               Call gensigma (wplot, mdf, optcompt, sigma)
               Call gensumrls (wplot, mdf, sumrls)
           ! write optical functions to file
               Call writeeps (iq, oct1, oct2, wplot, mdf, trim(fneps))
               Call writeloss (iq, wplot, loss(oct1, oct2, :), trim(fnloss))
               Call writesigma (iq, wplot, sigma, trim(fnsigma))
               if (tq0) Call writesumrls (iq, sumrls, trim(fnsumrules))
           ! end loop over optical components
            End Do
         End Do

         Call genfilname (basename='MOKE', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & iqmt=iq, filnam=fnmoke)
         
         If (input%xs%dfoffdiag.And.(nc.Eq.3).And.(m.Eq.1)) Then
            Call genkerr (wplot, mdf2, kerr)
            Call writekerr (iq, wplot, kerr, trim(fnmoke))
         End If      

      End Do ! m

  ! deallocate
      Deallocate (mdf, mdf1, mdf2, w, wr, wplot, loss, sigma)
      Deallocate (eps1, eps2, cf, kerr)
End Subroutine xslinopt
