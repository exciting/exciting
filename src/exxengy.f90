
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengy
use modmain
implicit none
! local variables
integer is,ia,nr,m1,m2
integer ik,ist1,ist2
real(8) evv,ecv,ecc
complex(8) zpchg,zt1
! allocatable arrays
complex(8), allocatable :: wfcr1(:,:,:)
complex(8), allocatable :: wfcr2(:,:,:)
complex(8), allocatable :: zrhomt(:,:)
complex(8), allocatable :: zvclmt(:,:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(wfcr1(lmmaxvr,nrcmtmax,2))
allocate(wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax))
allocate(zvclmt(lmmaxvr,nrcmtmax))
allocate(zfmt(lmmaxvr,nrcmtmax))
zpchg=0.d0
evv=0.d0
ecv=0.d0
ecc=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL
  write(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call exxengyk(ik,evv,ecv)
end do
!$OMP END DO
!$OMP END PARALLEL
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nr=nrcmt(is)
  do ia=1,natoms(is)
    do ist2=1,spnst(is)
      if (spcore(ist2,is)) then
        do m2=-spk(ist2,is),spk(ist2,is)-1
          call wavefcr(lradstp,is,ia,ist2,m2,nrcmtmax,wfcr2)
          do ist1=1,spnst(is)
            if (spcore(ist1,is)) then
              do m1=-spk(ist1,is),spk(ist1,is)-1
                call wavefcr(lradstp,is,ia,ist1,m1,nrcmtmax,wfcr1)
! calculate the complex overlap density
                call vnlrhomt(is,wfcr1(1,1,1),wfcr2(1,1,1),zrhomt)
                call vnlrhomt(is,wfcr1(1,1,2),wfcr2(1,1,2),zfmt)
                zrhomt(:,1:nr)=zrhomt(:,1:nr)+zfmt(:,1:nr)
! calculate the Coulomb potential
                call zpotclmt(lmaxvr,nr,rcmt(1,is),zpchg,lmmaxvr,zrhomt, &
                 zvclmt)
                zt1=zfmtinp(lmaxvr,nr,rcmt(1,is),lmmaxvr,zrhomt,zvclmt)
                ecc=ecc-0.5d0*dble(zt1)
              end do
! end loop over ist1
            end if
          end do
        end do
! end loop over ist2
      end if
    end do
! end loops over atoms and species
  end do
end do
! total exchange energy
engyx=evv+ecv+ecc
deallocate(wfcr1,wfcr2,zrhomt,zvclmt,zfmt)
return
end subroutine

