!
!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writevars
      Implicit None
Contains
!
!
      Subroutine writevars (un, iq, iqmt)
         Use modmain
         Use modinput
         Use modxs
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, iqmt, un
    ! local variables
         Character (10) :: dat, tim
         Integer :: iqt, iqmtt
         iqt = iq
         If (iqt .Eq. 0) iqt = iqmt
         iqmtt = iqmt
         If (iqmtt .Eq. 0) iqmtt = 1
    ! write prologue to file
         Call date_and_time (date=dat, time=tim)
         Write (un,*)
         Write (un, '("## Date (YYYY-MM-DD): ", a4, "-", a2, "-", a2)') &
        & dat (1:4), dat (5:6), dat (7:8)
         Write (un, '("## Time (hh:mm:ss)  : ", a2, ":", a2, ":", a2)') &
        & tim (1:2), tim (3:4), tim (5:6)
         Write (un, '("# version	    : ", i2.2, ".", i2.2, ".", i2.2)') version
         Write (un, '(a, 2f12.6)') '# efermi (H, eV)	   :', efermi, &
        & h2ev * efermi
         Write (un, '(a, 3f12.6)') '# vgqlmt		  :', &
        & input%xs%qpointset%qpoint(:, iqmtt)
         Write (un, '(a, 3i8)') '# ivgmt	       :', ivgmt (:, iqmtt)
         Write (un, '(a, 3f12.6)') '# vql		  :', vql (:, iqmtt)
         Write (un, '(a, a)') '# input%xs%tddft%fxctype	   :', &
        & trim(input%xs%tddft%fxctype)
         Write (un, '(a, f12.6)') '# input%xs%tddft%alphalrc		:', &
        & input%xs%tddft%alphalrc
         Write (un, '(a, f12.6)') '# input%xs%tddft%alphalrcdyn	:', &
        & input%xs%tddft%alphalrcdyn
         Write (un, '(a, f12.6)') '# input%xs%tddft%betalrcdyn	:', &
        & input%xs%tddft%betalrcdyn
         Write (un, '(a, l8)') '# input%xs%tddft%intraband	     :', &
        & input%xs%tddft%intraband
         Write (un, '(a, l8)') '# input%xs%tddft%aresdf	     :', &
        & input%xs%tddft%aresdf
         Write (un, '(a, l8)') '# input%xs%bse%aresbse       :', &
        & input%xs%bse%aresbse
         Write (un, '(a, l8)') '# input%xs%tddft%aresfxc	     :', &
        & input%xs%tddft%aresfxc
         Write (un, '(a, l8)') '# input%xs%tddft%torddf	     :', &
        & input%xs%tddft%torddf
         Write (un, '(a, l8)') '# input%xs%tddft%tordfxc	     :', &
        & input%xs%tddft%tordfxc
         Write (un, '(a, l8)') '# input%xs%tddft%acont	     :', &
        & input%xs%tddft%acont
         Write (un, '(a, i8)') '# input%xs%tddft%nwacont	     :', &
        & input%xs%tddft%nwacont
         Write (un, '(a, 2f12.6)') '# input%xs%broad (H, eV)	    :', &
        & input%xs%broad, h2ev * input%xs%broad
         Write (un, '(a, 2f12.6)') '# scissor (H, eV)	  :', &
        & input%xs%scissor, h2ev * input%xs%scissor
         Write (un, '(a, i8)') '# input%xs%energywindow%points		   :', &
        & input%xs%energywindow%points
         Write (un, '(a, i8)') '# ngq 	      :', ngq (iqmtt)
         Write (un, '(a, f12.6)') '# input%xs%gqmax		  :', &
        & input%xs%gqmax
         Write (un, '(a, f12.6)') '# input%groundstate%gmaxvr 	   :', &
        & input%groundstate%gmaxvr
         Write (un, '(a, f12.6)') '# input%groundstate%rgkmax 	   :', &
        & input%groundstate%rgkmax
         Write (un, '(a, f12.6)') '# gkmax		 :', gkmax
         Write (un, '(a, 3i8)') '# input%groundstate%ngridk		  :', &
        & input%groundstate%ngridk
         Write (un, '(a, 3f12.6)') '# input%groundstate%vkloff	    :', &
        & input%groundstate%vkloff
         Write (un, '(a, l8)') '# input%groundstate%reducek		:', &
        & input%groundstate%reducek
         Write (un, '(a, i8)') '# nmatmax	      :', nmatmax
         Write (un, '(a, i8)') '# ngkmax	      :', ngkmax
         Write (un, '(a, i8)') '# nlotot	      :', nlotot
         Write (un, '(a, i8)') '# nlomax	      :', nlomax
         Write (un, '(a, i8)') '# nst1	      :', nst1
         Write (un, '(a, i8)') '# nst2	      :', nst2
         Write (un, '(a, 4i8)') '# nstlbse         :', input%xs%BSE%nstlbse
         Write (un, '(a, i8)') '# nstsv	      :', nstsv
         Write (un, '(a, 2f12.6)') '# evlmincut (H, eV)  :', evlmincut, &
        & h2ev * evlmincut
         Write (un, '(a, 2f12.6)') '# evlmaxcut (H, eV)  :', evlmaxcut, &
        & h2ev * evlmaxcut
         Write (un, '(a, 2f12.6)') '# evlmin (H, eV)	   :', evlmin, &
        & h2ev * evlmin
         Write (un, '(a, 2f12.6)') '# evlmax (H, eV)	   :', evlmax, &
        & h2ev * evlmax
         Write (un, '(a, i5, 2f12.6)') '# evlhpo (H, eV)     :', &
        & istocc0, evlhpo, h2ev * evlhpo
         Write (un, '(a, i5, 2f12.6)') '# evllpu (H, eV)     :', &
        & istunocc0, evllpu, h2ev * evllpu
         Write (un, '(a, l8)') '# ksgap	      :', ksgap
         Write (un, '(a, i8)') '# input%groundstate%lmaxapw		:', &
        & input%groundstate%lmaxapw
         Write (un, '(a, i8)') '# input%xs%lmaxapwwf	       :', &
        & input%xs%lmaxapwwf
         Write (un, '(a, i8)') '# input%groundstate%lmaxmat		:', &
        & input%groundstate%lmaxmat
         Write (un, '(a, i8)') '# input%groundstate%lmaxvr		:', &
        & input%groundstate%lmaxvr
         Write (un, '(a, i8)') '# input%groundstate%lmaxinr		:', &
        & input%groundstate%lmaxinr
         Write (un, '(a, i8)') '# lolmax	      :', lolmax
         Write (un, '(a, i8)') '# input%xs%lmaxemat	       :', &
        & input%xs%lmaxemat
         Write (un, '(a, i8)') '# input%xs%tddft%lmaxalda	     :', &
        & input%xs%tddft%lmaxalda
         Write (un, '(a, i8)') '# input%groundstate%lradstep		 :', &
        & input%groundstate%lradstep
         Write (un, '(a, l8)') '# input%xs%tevout	       :', &
        & input%xs%tevout
!         Write (un, '(a, l8)') '# input%xs%fastemat	       :', &
!        & input%xs%fastemat
         Write (un, '(a, l8)') '# input%groundstate%nosym		:', &
        & input%groundstate%nosym
         Write (un, '(a, l8)') '# input%xs%dfoffdiag	       :', &
        & input%xs%dfoffdiag
         Write (un,*)
      End Subroutine writevars
!
End Module m_writevars
