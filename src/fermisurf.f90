
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fermisurf
use modmain
implicit none
! local variables
integer ik,ist,ist1,ist2,nst
integer ng(3),i1,i2,i3,j1,j2,j3
real(8) prd,v1(3),v2(3)
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
! get eigenvalues from file
do ik=1,nkpt
  call getevalsv(vkl(1,ik),evalsv(1,ik))
end do
! states to include in plot (task=101)
ist=(nstfv-nempty)*nspinor
ist1=max(ist-nstfsp/2,1)
ist2=min(ist+nstfsp/2,nstsv)
nst=min(ist2-ist1+1,30)
! produce Fermi surface plot
open(50,file='FERMISURF.OUT',action='WRITE',form='FORMATTED')
ng(:)=nup3d(:)*ngridk(:)
if (task.eq.100) then
  write(50,'(3I6," : grid size")') ng
else
  write(50,'(4I6," : grid size, number of states")') ng,nst
end if
do i3=0,nup3d(3)-1
  do j3=0,ngridk(3)-1
    do i2=0,nup3d(2)-1
      do j2=0,ngridk(2)-1
        do i1=0,nup3d(1)-1
          do j1=0,ngridk(1)-1
            v1(1)=dble(i1)+(dble(j1)+vkloff(1))/dble(ngridk(1))
            v1(2)=dble(i2)+(dble(j2)+vkloff(2))/dble(ngridk(2))
            v1(3)=dble(i3)+(dble(j3)+vkloff(3))/dble(ngridk(3))
            v2(:)=v1(1)*bvec(:,1)+v1(2)*bvec(:,2)+v1(3)*bvec(:,3)
            ik=ikmap(j1,j2,j3)
            if (task.eq.100) then
! write the product of eigenvalues minus the Fermi energy
              prd=1.d0
              do ist=1,nstsv
                prd=prd*(evalsv(ist,ik)-efermi)
              end do
              write(50,'(4G18.10)') v2,prd
            else
! write the eigenvalues minus the Fermi energy separately
              write(50,'(40F14.8)') v2,(evalsv(ist+ist1-1,ik)-efermi,ist=1,nst)
            end if
          end do
        end do
      end do
    end do
  end do
end do
close(50)
write(*,*)
write(*,'("Info(fermisurf):")')
write(*,'(" 3D Fermi surface data written to FERMISURF.OUT")')
if (task.eq.100) then
  write(*,'(" in terms of the product of eigenvalues minus the Fermi energy")')
else
  write(*,'(" in terms of separate eigenvalues minus the Fermi energy")')
end if
write(*,*)
return
end subroutine
