
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynqtor(dynq,dynr)
use modmain
implicit none
! arguments
complex(8), intent(in) :: dynq(3*natmtot,3*natmtot,nqpt)
complex(8), intent(out) :: dynr(3*natmtot,3*natmtot,ngridq(1)*ngridq(2) &
 *ngridq(3))
! local variables
integer ir,iq,i,j,n
integer isym,lspl,iv(3)
integer i1,i2,i3,j1,j2,j3
real(8) v1(3),v2(3),s(3,3),t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: dyns(:,:)
! external functions
real(8) r3taxi
external r3taxi
allocate(dyns(3*natmtot,3*natmtot))
dynr(:,:,:)=0.d0
! loop over q-vectors
do j1=0,ngridq(1)-1
  v1(1)=dble(j1)/dble(ngridq(1))
  do j2=0,ngridq(2)-1
    v1(2)=dble(j2)/dble(ngridq(2))
    do j3=0,ngridq(3)-1
      v1(3)=dble(j3)/dble(ngridq(3))
      iq=iqmap(j1,j2,j3)
! rotate and add the dynamical matrix of the reduced q-point with all symmetries
      n=0
      dyns(:,:)=0.d0
      do isym=1,nsymcrys
        lspl=lsplsymc(isym)
        s(:,:)=dble(symlat(:,:,lspl))
        call r3mtv(s,vql(1,iq),v2)
        call r3frac(epslat,v2,iv)
        if (r3taxi(v1,v2).lt.epslat) then
          call dynsymapp(isym,vql(1,iq),dynq(1,1,iq),dyns)
          n=n+1
        end if
      end do
      if (n.eq.0) then
        write(*,*)
        write(*,'("Error(dynqtor): vector ",3G18.10)') v1
        write(*,'(" cannot be mapped to reduced q-point set")')
        write(*,*)
        stop
      end if
      t1=1.d0/dble(n)
      dyns(:,:)=t1*dyns(:,:)
! loop over R-vectors
      ir=0
      do i3=ngridq(3)/2-ngridq(3)+1,ngridq(3)/2
        do i2=ngridq(2)/2-ngridq(2)+1,ngridq(2)/2
          do i1=ngridq(1)/2-ngridq(1)+1,ngridq(1)/2
            ir=ir+1
            t1=twopi*(v1(1)*dble(i1)+v1(2)*dble(i2)+v1(3)*dble(i3))
            zt1=cmplx(cos(t1),sin(t1),8)
            do i=1,3*natmtot
              do j=1,3*natmtot
                dynr(i,j,ir)=dynr(i,j,ir)+zt1*dyns(i,j)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
t1=1.d0/dble(ngridq(1)*ngridq(2)*ngridq(3))
dynr(:,:,:)=t1*dynr(:,:,:)
deallocate(dyns)
return
end subroutine

