! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematradou
! !INTERFACE:
Subroutine ematradou (iq,ik, igq, ngp, apwalm, evecfvo, integral)
! !USES:
	Use modmain
	Use modinput
	Use modxs
	Use m_getunit 
	Use modxas
! !INPUT/OUTPUT PARAMETERS:
!   iq       : q-point position (in,integer)
!   ik       : k-point position (in,integer)
!   ngp      : number of G+p-vectors (in,integer)
!   apwalm   : APW matching coefficients
!              (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfvo  : first-variational eigenvector (in,complex(nmatmax))
!   integral : radial planewave integral 
!              (out, complex(lmaxemat+1,lmmaxapw,nxas,sta2:sto2))
! !DESCRIPTION:
!   Calculates the radial integral part $R^{ll'}_{\mu m(\mathbf{k}+\mathbf{q})}(\mathbf{q}+\mathbf{G})$
!   of the planewave matrix element between a 
!   core state and a conduction state. See Vorwerk's Master thesis for more details. 
!
! !REVISION HISTORY:
!  Based on the subroutine ematqk.F90
!  Created November 2015 (Christian Vorwerk)
!EOP
!BOC      

	Implicit none
	Integer, Intent (In) :: iq, ik, igq
	integer, intent(in)     :: ngp
	complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
	complex(8), intent(in)  :: evecfvo(nmatmax,nstfv)
	Complex(8), Intent (Out) :: integral(input%xs%lmaxemat+1,lmmaxapw,nxas,sta2:sto2)
	! local variables
	Integer :: is, ia, ir, nr, lmax2, n1, n2, l2,l3,m3,lm3, irc
	Real(8) :: t1, t2
	Real (8), Allocatable :: jl (:, :), jhelp (:)
	Real (8) :: r2 (nrmtmax), fr2 (nrcmtmax), fr3(nrcmtmax), gr (nrcmtmax), cf (3,nrcmtmax)
	Real (8), allocatable :: fr1(:,:,:)
	complex(8), allocatable :: wfmt(:,:,:)
	lmax2 = input%xs%lmaxemat
	is=input%xs%bse%xasspecies
	ia=input%xs%bse%xasatom

	Allocate (fr1(0:lmax2,nxas,nrcmtmax))
	Allocate (jl(nrmtmax,0:lmax2))
	Allocate (jhelp(0:lmax2))	
	allocate(wfmt(lmmaxapw,nrcmtmax,nstfv))
    nr = nrmt (is)
     irc=0
! Calculate product of core radial wavefunction and Bessel functions    
    Do ir = 1,nrmt(is),input%groundstate%lradstep
		irc=irc+1
		! calculate r^2
		r2 (ir) = spr (ir, is) ** 2
		! calculate spherical Bessel functions of first kind j_l(|G+q|r_a)
		Call sbessel (lmax2, gqc(igq, iq)*spr(ir, is), jhelp)
		jl (ir,:) = jhelp (:)
		do n1=1,nxas
			fr1(:,n1,irc)=r2(ir)*jl(ir,:)*ucore(ir,n1+xasstart-1)
			
		end do
	End Do  
	Do n1=1,nxas
		Do n2=sta2,sto2
		! Obtain radial wavefunction of the conduction state		
			call wavefmt(input%groundstate%lradstep, &
			&  input%groundstate%lmaxapw,input%xs%bse%xasspecies,input%xs%bse%xasatom,ngp,apwalm, &
			&  evecfvo(:,n2+istunocc0-1),lmmaxapw,wfmt(:,:,n2+istunocc0-1))
             
			Do l2=0, lmax2
				Do l3=0,input%xs%lmaxapw
					Do m3=-l3,l3
					lm3=idxlm(l3,m3)
						Do irc=1,nrcmt(input%xs%bse%xasspecies)
							fr2(irc)=fr1(l2,n1,irc)*dble(wfmt(lm3,irc,n2+istunocc0-1))
							fr3(irc)=fr1(l2,n1,irc)*aimag(wfmt(lm3,irc,n2+istunocc0-1))
						End Do
						! Radial integration
						Call fderiv (-1, nrcmt(input%xs%bse%xasspecies), spr(1, is), fr2, gr, cf)
						t1=gr (nrcmt(input%xs%bse%xasspecies))
						Call fderiv (-1, nrcmt(input%xs%bse%xasspecies), spr(1, is), fr3, gr, cf)
						t2=gr (nrcmt(input%xs%bse%xasspecies))						
!						integral(l2+1,lm3,n1,n2) = gr (nrcmt(input%xs%bse%xasspecies))
						integral(l2+1,lm3,n1,n2) = cmplx(t1,t2,8)
					End Do
				End Do
			End Do
		End Do
	End Do
	! Deallocate local variables
	deallocate (fr1)
	deallocate(jl)
	deallocate (jhelp)
	deallocate (wfmt)
End Subroutine ematradou
! EOC
