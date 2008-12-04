
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dyndiag(dynp,w,ev)
use modmain
implicit none
! arguments
complex(8), intent(in) :: dynp(3*natmtot,3*natmtot)
real(8), intent(out) :: w(3*natmtot)
complex(8), intent(out) :: ev(3*natmtot,3*natmtot)
! local variables
integer is,ia,ip,js,ja,jp
integer i,j,n,lwork,info
real(8) t1
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
! number of phonon branches
n=3*natmtot
ev(:,:)=0.d0
i=0
do is=1,nspecies
  do ia=1,natoms(is)
    do ip=1,3
      i=i+1
      j=0
      do js=1,nspecies
! mass factor
        if ((spmass(is).le.0.d0).or.(spmass(js).le.0.d0)) then
! infinite mass
          t1=0.d0
        else
          t1=1.d0/sqrt(spmass(is)*spmass(js))
        end if
        do ja=1,natoms(js)
          do jp=1,3
            j=j+1
            if (i.le.j) then
! use Hermitian average of dynamical matrix
              ev(i,j)=0.5d0*t1*(dynp(i,j)+conjg(dynp(j,i)))
            end if
          end do
        end do
      end do
    end do
  end do
end do
allocate(rwork(3*n))
lwork=2*n
allocate(work(lwork))
call zheev('V','U',n,ev,n,w,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(dyndiag): diagonalisation failed")')
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
do i=1,n
  if (w(i).ge.0.d0) then
    w(i)=sqrt(w(i))
  else
    w(i)=-sqrt(abs(w(i)))
  end if
end do
deallocate(rwork,work)
return
end subroutine

