
! Copyright (C) 2002-2005 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine effmass
use modmain
implicit none
! local variables
integer, parameter :: lwork=10
integer ik,ik0,ist,info
integer i,j,k,l,m,n
integer i1,i2,i3,j1,j2,j3
real(8) em(3,3),v1(3),v2(3)
real(8) w(3),work(lwork)
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: a(:,:)
real(8), allocatable :: b(:,:,:,:)
real(8), allocatable :: c(:,:,:)
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
allocate(ipiv(nkpt))
allocate(a(nkpt,nkpt))
n=2*ndspem+1
allocate(b(0:n-1,0:n-1,0:n-1,nstsv))
allocate(c(0:n-1,0:n-1,0:n-1))
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
ik0=0
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP PRIVATE(i1,i2,i3,j1,j2,j3,ist)
!$OMP DO
do ik=1,nkpt
  allocate(evalfv(nstfv,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  i1=ivk(1,ik); i2=ivk(2,ik); i3=ivk(3,ik)
  if ((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0)) ik0=ik
! solve the first- and second-variational secular equations
  call seceqn(ik,evalfv,evecfv,evecsv)
! copy eigenvalues to new array
  j1=i1+ndspem; j2=i2+ndspem; j3=i3+ndspem
  do ist=1,nstsv
    b(j1,j2,j3,ist)=evalsv(ist,ik)
  end do
  deallocate(evalfv,evecfv,evecsv)
end do
!$OMP END DO
!$OMP END PARALLEL
! set up polynomial matrix
i=0
do i3=-ndspem,ndspem
  do i2=-ndspem,ndspem
    do i1=-ndspem,ndspem
      i=i+1
      v1(1)=dble(i1); v1(2)=dble(i2); v1(3)=dble(i3)
      v1(:)=v1(:)*deltaem
      j=0
      v2(3)=1.d0
      do j3=0,n-1
        v2(2)=1.d0
        do j2=0,n-1
          v2(1)=1.d0
          do j1=0,n-1
            j=j+1
            a(i,j)=v2(1)*v2(2)*v2(3)
            v2(1)=v2(1)*v1(1)
          end do
          v2(2)=v2(2)*v1(2)
        end do
        v2(3)=v2(3)*v1(3)
      end do
    end do
  end do
end do
! solve for the polynomial coefficients
call dgesv(nkpt,nstsv,a,nkpt,ipiv,b,nkpt,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(effmass): could not determine polynomial coefficients")')
  write(*,'(" DGESV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
open(50,file='EFFMASS.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("(effective mass matrices are in Cartesian coordinates)")')
write(50,*)
write(50,'("k-point (lattice coordinates) :")')
write(50,'(3G18.10)') vklem
write(50,*)
write(50,'("k-point (Cartesian coordinates) :")')
call r3mv(bvec,vklem,v1)
write(50,'(3G18.10)') v1
! begin loop over states
do ist=1,nstsv
! compute matrix of derivatives with respect to k-vector
  do k=1,3
    do l=1,3
      c(:,:,:)=b(:,:,:,ist)
      do i=1,2
        if (i.eq.1) then
          m=k
        else
          m=l
        end if
        if (m.eq.1) then
          do j=0,n-2
            c(j,:,:)=dble(j+1)*c(j+1,:,:)
          end do
          c(n-1,:,:)=0.d0
        else if (m.eq.2) then
          do j=0,n-2
            c(:,j,:)=dble(j+1)*c(:,j+1,:)
          end do
          c(:,n-1,:)=0.d0
        else if (m.eq.3) then
          do j=0,n-2
            c(:,:,j)=dble(j+1)*c(:,:,j+1)
          end do
          c(:,:,n-1)=0.d0
        end if
      end do
! derivative evaluated at zero
      em(k,l)=c(0,0,0)
    end do
  end do
  write(50,*)
  write(50,*)
  write(50,'("State, eigenvalue : ",I6,G18.10)') ist,evalsv(ist,ik0)
  write(50,*)
  write(50,'(" matrix of eigenvalue derivatives with respect to k :")')
  do i=1,3
    write(50,'(3G18.10)') (em(i,j),j=1,3)
  end do
  write(50,'(" trace : ",G18.10)') em(1,1)+em(2,2)+em(3,3)
! invert derivative matrix
  call r3minv(em,em)
  write(50,*)
  write(50,'(" effective mass tensor (inverse of derivative matrix) :")')
  do i=1,3
    write(50,'(3G18.10)') (em(i,j),j=1,3)
  end do
  write(50,'(" trace : ",G18.10)') em(1,1)+em(2,2)+em(3,3)
! find the eigenvalues
  call dsyev('N','U',3,em,3,w,work,lwork,info)
  write(50,'(" eigenvalues :")')
  write(50,'(3G18.10)') w
! end loop over states
end do
close(50)
write(*,*)
write(*,'("Info(effmass):")')
write(*,'(" effective mass tensor for each state written to EFFMASS.OUT")')
write(*,'(" for k-point (lattice) ",3G18.10)') vklem
write(*,*)
deallocate(ipiv,a,b,c)
return
end subroutine

