
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gencore
! !INTERFACE:
subroutine gencore
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the core radial wavefunctions, eigenvalues and densities. The
!   radial Dirac equation is solved in the spherical part of the effective
!   potential to which the atomic potential has been appended for
!   $r>R^{\rm MT}$. In the case of spin-polarised calculations, the effective
!   magnetic field is ignored.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia1,ia2,ias1,ias2,ist,ir
real(8) t1
! automatic arrays
logical done(natmmax)
real(8) vr(spnrmax)
do is=1,nspecies
  done(:)=.false.
  do ia1=1,natoms(is)
    if (.not.done(ia1)) then
      ias1=idxas(ia1,is)
      vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias1)*y00
! append the effective potential from the atomic calculation
      t1=vr(nrmt(is))-spvr(nrmt(is),is)
      do ir=nrmt(is)+1,spnr(is)
        vr(ir)=spvr(ir,is)+t1
      end do
      rhocr(:,ias1)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ir,t1)
!$OMP DO
      do ist=1,spnst(is)
        if (spcore(ist,is)) then
! solve the Dirac equation
          call rdirac(spn(ist,is),spl(ist,is),spk(ist,is),nprad,spnr(is), &
           spr(1,is),vr,evalcr(ist,ias1),rwfcr(1,1,ist,ias1), &
           rwfcr(1,2,ist,ias1))
          t1=spocc(ist,is)
!$OMP CRITICAL
          do ir=1,spnr(is)
            rhocr(ir,ias1)=rhocr(ir,ias1) &
             +t1*(rwfcr(ir,1,ist,ias1)**2+rwfcr(ir,2,ist,ias1)**2)
          end do
!$OMP END CRITICAL
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      do ir=1,spnr(is)
        rhocr(ir,ias1)=rhocr(ir,ias1)/(fourpi*spr(ir,is)**2)
      end do
      done(ia1)=.true.
! copy to equivalent atoms
      do ia2=1,natoms(is)
        if ((.not.done(ia2)).and.(nsymeqat(ia2,ia1,is).gt.0)) then
          ias2=idxas(ia2,is)
          do ist=1,spnst(is)
            if (spcore(ist,is)) then
              evalcr(ist,ias2)=evalcr(ist,ias1)
              rwfcr(1:spnr(is),:,ist,ias2)=rwfcr(1:spnr(is),:,ist,ias1)
            end if
          end do
          rhocr(1:spnr(is),ias2)=rhocr(1:spnr(is),ias1)
          done(ia2)=.true.
        end if
      end do
    end if
  end do
end do
return
end subroutine
!EOC

