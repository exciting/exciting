



! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writevars
  implicit none
contains


subroutine writevars(un, iq, iqmt)
    use modmain
    use modinput
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: iq, iqmt, un
    ! local variables
    character(10) :: dat, tim
    integer :: iqt, iqmtt
    iqt=iq
    if (iqt.eq.0) iqt=iqmt
    iqmtt=iqmt
    if (iqmtt.eq.0) iqmtt=1
    ! write prologue to file
    call date_and_time(date=dat, time=tim)
    write(un, *)
    write(un, '("## Date (YYYY-MM-DD): ", a4, "-", a2, "-", a2)') &
	 dat(1:4), dat(5:6), dat(7:8)
    write(un, '("## Time (hh:mm:ss)  : ", a2, ":", a2, ":", a2)') &
	 tim(1:2), tim(3:4), tim(5:6)
    write(un, '("# version	     : ", i1.1, ".", i1.1, ".", i3.3)') version
    write(un, '("# version (xs)      : ", i1.1, ".", i3.3)') versionxs
    write(un, '(a, 2f12.6)') '# efermi (H, eV)	   :', efermi, h2ev*efermi
    write(un, '(a, 3f12.6)') '# vgqlmt		  :', input%xs%qpointset%qpoint(:, iqmtt)
    write(un, '(a, 3i8)') '# ivgmt	       :', ivgmt(:, iqmtt)
    write(un, '(a, 3f12.6)') '# vql		  :', vql(:, iqmtt)
    write(un, '(a, i8)') '# input%xs%tddft%fxctypenumber	   :', input%xs%tddft%fxctypenumber
    write(un, '(a, f12.6)') '# input%xs%tddft%alphalrc		:', input%xs%tddft%alphalrc
    write(un, '(a, f12.6)') '# input%xs%tddft%alphalrcdyn	:', input%xs%tddft%alphalrcdyn
    write(un, '(a, f12.6)') '# input%xs%tddft%betalrcdyn	:', input%xs%tddft%betalrcdyn
    write(un, '(a, l8)') '# input%xs%tddft%intraband	     :', input%xs%tddft%intraband
    write(un, '(a, l8)') '# input%xs%tddft%aresdf	     :', input%xs%tddft%aresdf
    write(un, '(a, l8)') '# input%xs%tddft%aresfxc	     :', input%xs%tddft%aresfxc
    write(un, '(a, l8)') '# input%xs%tddft%torddf	     :', input%xs%tddft%torddf
    write(un, '(a, l8)') '# input%xs%tddft%tordfxc	     :', input%xs%tddft%tordfxc
    write(un, '(a, l8)') '# input%xs%tddft%acont	     :', input%xs%tddft%acont
    write(un, '(a, i8)') '# input%xs%tddft%nwacont	     :', input%xs%tddft%nwacont
    write(un, '(a, 2f12.6)') '# input%xs%broad (H, eV)	    :', input%xs%broad, h2ev*input%xs%broad
    write(un, '(a, 2f12.6)') '# scissor (H, eV)	  :',&
    &input%xs%scissor, h2ev * input%xs%scissor
    write(un, '(a, i8)') '# input%xs%dosWindow%points		   :', input%xs%dosWindow%points
    write(un, '(a, i8)') '# ngq 	      :', ngq(iqmtt)
    write(un, '(a, f12.6)') '# input%xs%gqmax		  :', input%xs%gqmax
    write(un, '(a, f12.6)') '# input%groundstate%gmaxvr 	   :', input%groundstate%gmaxvr
    write(un, '(a, f12.6)') '# input%groundstate%rgkmax 	   :', input%groundstate%rgkmax
    write(un, '(a, f12.6)') '# gkmax		 :', gkmax
    write(un, '(a, 3i8)') '# input%groundstate%ngridk		  :', input%groundstate%ngridk
    write(un, '(a, 3f12.6)') '# input%groundstate%vkloff	    :', input%groundstate%vkloff
    write(un, '(a, l8)') '# input%groundstate%reducek		:', input%groundstate%reducek
    write(un, '(a, i8)') '# nmatmax	      :', nmatmax
    write(un, '(a, i8)') '# ngkmax	      :', ngkmax
    write(un, '(a, i8)') '# nlotot	      :', nlotot
    write(un, '(a, i8)') '# nlomax	      :', nlomax
    write(un, '(a, i8)') '# nst1	      :', nst1
    write(un, '(a, i8)') '# nst2	      :', nst2
    write(un, '(a, 2i8)') '# nstlbse	       :', nbfbse
    write(un, '(a, i8)') '# nstsv	      :', nstsv
    write(un, '(a, 2f12.6)') '# evlmincut (H, eV)  :', evlmincut, h2ev*evlmincut
    write(un, '(a, 2f12.6)') '# evlmaxcut (H, eV)  :', evlmaxcut, h2ev*evlmaxcut
    write(un, '(a, 2f12.6)') '# evlmin (H, eV)	   :', evlmin, h2ev*evlmin
    write(un, '(a, 2f12.6)') '# evlmax (H, eV)	   :', evlmax, h2ev*evlmax
    write(un, '(a, i5, 2f12.6)') '# evlhpo (H, eV)     :', istocc0, evlhpo, h2ev*evlhpo
    write(un, '(a, i5, 2f12.6)') '# evllpu (H, eV)     :', istunocc0, evllpu, &
	 h2ev * evllpu
    write(un, '(a, l8)') '# ksgap	      :', ksgap
    write(un, '(a, i8)') '# input%groundstate%lmaxapw		:', input%groundstate%lmaxapw
    write(un, '(a, i8)') '# input%xs%lmaxapwwf	       :', input%xs%lmaxapwwf
    write(un, '(a, i8)') '# input%groundstate%lmaxmat		:', input%groundstate%lmaxmat
    write(un, '(a, i8)') '# input%groundstate%lmaxvr		:', input%groundstate%lmaxvr
    write(un, '(a, i8)') '# input%groundstate%lmaxinr		:', input%groundstate%lmaxinr
    write(un, '(a, i8)') '# lolmax	      :', lolmax
    write(un, '(a, i8)') '# input%xs%lmaxemat	       :', input%xs%lmaxemat
    write(un, '(a, i8)') '# input%xs%tddft%lmaxalda	     :', input%xs%tddft%lmaxalda
    write(un, '(a, i8)') '# input%groundstate%lradstep		 :', input%groundstate%lradstep
    write(un, '(a, l8)') '# input%xs%tevout	       :', input%xs%tevout
    write(un, '(a, l8)') '# input%xs%fastpmat	       :', input%xs%fastemat
    write(un, '(a, l8)') '# input%xs%fastemat	       :', input%xs%fastemat
    write(un, '(a, l8)') '# input%groundstate%nosym		:', input%groundstate%nosym
    write(un, '(a, l8)') '# input%xs%dfoffdiag	       :', input%xs%dfoffdiag
    write(un, *)
  end subroutine writevars

end module m_writevars
