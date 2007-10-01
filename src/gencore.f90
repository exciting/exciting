
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
integer is,ia,ja,ias,jas,ist,ir
real(8) t1
! automatic arrays
logical done(natmmax)
real(8) vr(spnrmax)
do is=1,nspecies
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias)*y00
! append the effective potential from the atomic calculation
      t1=vr(nrmt(is))-spvr(nrmt(is),is)
      do ir=nrmt(is)+1,spnr(is)
        vr(ir)=spvr(ir,is)+t1
      end do
      rhocr(:,ias)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ir,t1)
!$OMP DO
      do ist=1,spnst(is)
        if (spcore(ist,is)) then
! solve the Dirac equation
          call rdirac(spn(ist,is),spl(ist,is),spk(ist,is),nprad,spnr(is), &
           spr(1,is),vr,evalcr(ist,ias),rwfcr(1,1,ist,ias),rwfcr(1,2,ist,ias))
          t1=spocc(ist,is)
!$OMP CRITICAL
          do ir=1,spnr(is)
! add to the core density
            rhocr(ir,ias)=rhocr(ir,ias) &
             +t1*(rwfcr(ir,1,ist,ias)**2+rwfcr(ir,2,ist,ias)**2)
          end do
!$OMP END CRITICAL
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      do ir=1,spnr(is)
        rhocr(ir,ias)=rhocr(ir,ias)/(fourpi*spr(ir,is)**2)
      end do
      done(ia)=.true.
! copy to equivalent atoms
      do ja=1,natoms(is)
        if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
          jas=idxas(ja,is)
          do ist=1,spnst(is)
            if (spcore(ist,is)) then
              evalcr(ist,jas)=evalcr(ist,ias)
              rwfcr(1:spnr(is),:,ist,jas)=rwfcr(1:spnr(is),:,ist,ias)
            end if
          end do
          rhocr(1:spnr(is),jas)=rhocr(1:spnr(is),ias)
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
return
end subroutine
!EOC

