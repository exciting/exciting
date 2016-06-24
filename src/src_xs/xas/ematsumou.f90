! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematsumou
! !INTERFACE:
Subroutine ematsumou (iq,ik, igq, integral, xi)
! !USES:
	Use modmain
	Use modinput
	Use modxs
	Use m_getunit 
	Use modxas
! !INPUT/OUTPUT PARAMETERS:
!   iq       : q-point position (in,integer)
!   ik       : k-point position (in,integer)
!   igq		 : (q+G)-point position (in, integer)
!   integral : radial planewave integral 
!              (in, complex(lmaxemat+1,lmmaxapw,nxas,sta2:sto2))
!	xi       : planewave matrix element (inout,complex(nxas, sta2:sto2, ngq(iq))
! !DESCRIPTION:
! Calculates planewave matrix elements between a core and a conduction state, using the radial integrals.
! Fore more information, see Christian Vorwerk's Master thesis, Eq. (5.20).
! !REVISION HISTORY:
!  Based on the subroutine ematqkgmt.F90
!  Created November 2015 (Christian Vorwerk)
!EOP
!BOC      

	Implicit none
	Integer, Intent (In) :: iq, ik, igq
	Complex(8), Intent (In) :: integral(input%xs%lmaxemat+1,lmmaxapw,nxas,sta2:sto2)
	Complex(8), Intent (InOut):: xi(nxas, sta2:sto2, ngq(iq))
	! local variables
	Integer :: n1, n2, l2, lmax2, m2, lm2, l3, m3, lm3, ias,ia, is
	Complex(8) :: prefactor
	Real (8) :: vk (3)
	! Setting xioo to zero
	xi(:,:,igq)=zzero
	is=input%xs%bse%xasspecies
	ia=input%xs%bse%xasatom
	ias=idxas(ia,is)
	lmax2=input%xs%lmaxemat

	Do n1=1,nxas
		Do n2=sta2,sto2
			Do l2=0,lmax2
				Do m2=-l2,l2
					lm2=idxlm(l2,m2)
					Do l3=0, input%groundstate%lmaxapw
						Do m3=-l3,l3
							lm3=idxlm(l3,m3)
							xi(n1,n2,igq)=xi(n1,n2,igq)+conjg (zil(l2))*integral(l2+1,lm3,n1,n2)*conjg &
									& (ylmgq(lm2, igq, iq)) * xsgntou &
									& (n1, lm2, lm3)
						End Do
					End Do
				End Do
			End Do
		End Do
	End Do
	vk (:) = vkl (:, ik)
	! Multiply with Structure Factor
	prefactor=fourpi*conjg(sfacgq(igq, ias, iq))
	xi(:,:,igq)=xi(:,:,igq)*prefactor
End Subroutine ematsumou
! EOC
