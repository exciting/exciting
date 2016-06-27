! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematsumoo
! !INTERFACE:
Subroutine ematsumoo (iq,ik, igq, integral, xi)
! !USES:
	Use modmain
	Use modinput
	Use modxs
	Use modxas
! !INPUT/OUTPUT PARAMETERS:
!   iq       : q-point position (in,integer)
!   ik       : k-point position (in,integer)
!   igq		 : (q+G)-point position (in, integer)
!   integral : radial planewave integral 
!              (in, complex(lmaxemat+1,nxas,nxas))
!	xi       : planewave matrix element (inout,complex(nxas, sta2:sto2, ngq(iq))
! !DESCRIPTION:
! Calculates planewave matrix elements between two core states, using the radial integrals.
! Fore more information, see Christian Vorwerk's Master thesis, Eq. (5.16).
! !REVISION HISTORY:
!  Based on the subroutine ematqkgmt.F90
!  Created November 2015 (Christian Vorwerk)
!EOP
!BOC      
	Implicit none
	Integer, Intent (In) :: iq, ik, igq
	Complex(8), Intent (In) :: integral(input%xs%lmaxemat+1,nxas,nxas)
	Complex(8), Intent (InOut):: xi(nxas, nxas, ngq(iq))
	! local variables
	Integer :: n1, n2, l2, lmax2, m2, lm2, ias,ia, is
	Complex(8) :: prefactor, prefactor2, inter(nxas,nxas)
	Real (8) :: vq (3)
	! Setting xi to zero
	inter(:,:)=zzero
	is=input%xs%bse%xasspecies
	ia=input%xs%bse%xasatom
	ias=idxas(ia,is)
	lmax2=input%xs%lmaxemat
	
	Do n1=1,nxas
		Do n2=1,nxas
			Do l2=0,lmax2
				Do m2=-l2,l2
					lm2=idxlm(l2,m2)
                    inter(n1,n2)=inter(n1,n2)+conjg (zil(l2))*integral(l2+1,n1,n2)*conjg (ylmgq(lm2, igq, iq)) * xsgntoo (n1, lm2, n2)
                    
				End Do
			End Do
		End Do
	End Do
	vq (:) = vql (:, iq)
	prefactor=fourpi*conjg(sfacgq(igq, ias, iq))
	xi(:,:,igq)=prefactor*inter(:,:)
End Subroutine ematsumoo
! EOC
