! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematradoo
! !INTERFACE:
Subroutine ematradoo (iq,ik, igq, integral)
! !USES:
	Use modmain
	Use modinput
	Use modxs
	Use m_getunit 
	Use modxas
! !INPUT/OUTPUT PARAMETERS:
!   iq       : q-point position (in,integer)
!   ik       : k-point position (in,integer)
!   igq      : (q+G)-point position (in,integer)
!   integral : radial planewave integral 
!              (out, complex(lmaxemat+1,nxas,nxas))
! !DESCRIPTION:
!   Calculates the radial integral part $R^{l}_{\mu \mu'}(\mathbf{q}+\mathbf{G})$
!   of the planewave matrix element between two core states. 
!   See Vorwerk's Master thesis for more details. 
!
! !REVISION HISTORY:
!  Based on the subroutine ematqk.F90
!  Created November 2015 (Christian Vorwerk)
!EOP
!BOC      

	Implicit none
	Integer, Intent (In) :: iq, ik, igq
	Complex(8), Intent (Out) :: integral(input%xs%lmaxemat+1,nxas,nxas)
	! local variables
	Integer :: is, ia, ir, nr, lmax2, n1, n2, l2
	Real(8) :: t1
	Real (8), Allocatable :: jl (:, :), jhelp (:)
	Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3,nrmtmax)
	is=input%xs%bse%xasspecies
	ia=input%xs%bse%xasatom
	lmax2 = input%xs%lmaxemat

	Allocate (jl(nrmtmax,0:lmax2))
	Allocate (jhelp(0:lmax2))
    nr = nrmt (is)
    Do ir = 1, nr
		! calculate r^2
		r2 (ir) = spr (ir, is) ** 2
		! calculate spherical Bessel functions of first kind j_l(|G+q|r_a)
		Call sbessel (lmax2, gqc(igq, iq)*spr(ir, is), jhelp)
		jl (ir,:) = jhelp (:)
	End Do  
	Do n1=1,nxas
		Do n2=1,nxas
			Do l2=0, lmax2
				Do ir=1,nr
					t1=ucore(ir,n1+xasstart-1)*ucore(ir,n2+xasstart-1)*r2(ir)
					fr(ir)=t1*jl(ir,l2)
				End Do
				Call fderiv (-1, nr, spr(1, is), fr, gr, cf)
                integral(l2+1,n1,n2) = gr (nr)
			End Do
		End Do
	End Do
End Subroutine ematradoo
! EOC
