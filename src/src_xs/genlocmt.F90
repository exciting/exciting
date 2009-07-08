

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine genlocmt(ngp, isti, istf, evecfv, wfcmt)
  use modmain
  implicit none
  ! arguments
  integer, intent(in) :: ngp, isti, istf
  complex(8), intent(in) :: evecfv(nmatmax, nstfv, nspnfv)
  complex(8), intent(out) :: wfcmt(istf-isti+1, nlomax, -lolmax:lolmax, natmtot)
  ! local variables
  integer :: ist, istc, is, ia, ias, i, ilo, l, m, lm
  do istc=isti, istf
     ist=istc-isti+1
     do is=1, nspecies
	do ia=1, natoms(is)
	   ias=idxas(ia, is)
	   do ilo=1, nlorb(is)
	      l=lorbl(ilo, is)
	      do m=-l, l
		 lm=idxlm(l, m)
		 i=idxlo(lm, ilo, ias)
		 wfcmt(ist, ilo, m, ias)=evecfv(ngp+i, istc, 1)
	      end do
	   end do
	end do
     end do
  end do
end subroutine genlocmt
