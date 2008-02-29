
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmvaryc(sum)
! calculates new evecsv using the gradient of energy w.r.t. evecsv
use modmain
implicit none
! arguments
real(8), intent(out) :: sum
! local variables
integer ik,ist1,ist2
real(8) t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: dedc(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: evecsvt(:)
! external functions
real(8) dznrm2
complex(8) zdotc
external dznrm2,zdotc
allocate(dedc(nstsv,nstsv,nkpt))
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL
  write(*,'("Info(rdmvaryc): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! compute the derivative w.r.t. evecsv
  call rdmdedc(ik,dedc(1,1,ik))
end do
!$OMP END DO
!$OMP END PARALLEL
allocate(evecsv(nstsv,nstsv))
allocate(evecsvt(nstsv))
sum=0.d0
do ik=1,nkpt
! get the eigenvectors from file
  call getevecsv(vkl(1,ik),evecsv)
! calculate new evecsv
  evecsv(:,:)=evecsv(:,:)-taurdmc*dedc(:,:,ik)
! othogonalise evecsv (Gram-Schmidt)
  do ist1=1,nstsv
    evecsvt(:)=evecsv(:,ist1)
    do ist2=1,ist1-1
      zt1=zdotc(nstsv,evecsv(1,ist2),1,evecsv(1,ist1),1)
      evecsvt(:)=evecsvt(:)-zt1*evecsv(:,ist2)
    end do
    t1=dznrm2(nstsv,evecsvt,1)
    t1=1.d0/t1
    evecsv(:,ist1)=t1*evecsvt(:)
  end do
! write new evecsv to file
  call putevecsv(ik,evecsv)
! convergence check
  do ist1=1,nstsv
    do ist2=1,nstsv
      zt1=dedc(ist1,ist2,ik)
      sum=sum+wkpt(ik)*(dble(zt1)**2+aimag(zt1)**2)
    end do
  end do
! end loop over k-points
end do
deallocate(dedc,evecsv,evecsvt)
return
end subroutine

