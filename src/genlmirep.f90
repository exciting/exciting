
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genlmirep(lmax,ld,elm,ulm)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(out) :: elm(ld,natmtot)
complex(8), intent(out) :: ulm(ld,ld,natmtot)
! local variables
integer is,ia,ias
integer lmmax,i,j,l,lm,n
integer isym,lspl,info,lwork
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: ulat(:,:,:)
complex(8), allocatable :: a(:,:),b(:,:)
complex(8), allocatable :: h(:,:)
complex(8), allocatable :: work(:)
lmmax=(lmax+1)**2
allocate(rwork(3*lmmax))
allocate(ulat(lmmax,lmmax,nsymlat))
allocate(a(lmmax,lmmax),b(lmmax,lmmax))
allocate(h(lmmax,lmmax))
lwork=2*lmmax
allocate(work(lwork))
! construct (l,m) rotation matrix for each lattice symmetry
a(:,:)=0.d0
do i=1,lmmax
  a(i,i)=1.d0
end do
do isym=1,nsymlat
  call rotzflm(symlatc(:,:,isym),lmax,lmmax,lmmax,a,ulat(:,:,isym))
end do
! set up quasi-random symmetric matrix H
h(:,:)=0.d0
do l=0,lmax
  n=2*l+1
  lm=idxlm(l,-l)
  do i=lm,lm+n-1
    do j=i,lm+n-1
      h(i,j)=dble(i*j)
      h(j,i)=h(i,j)
    end do
  end do
end do
! loop over species and atoms
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! symmetrise H with site symmetries
    b(:,:)=0.d0
    do isym=1,nsymsite(ias)
! spatial rotation element in lattice point group
      lspl=lsplsyms(isym,ias)
! apply lattice symmetry as U*H*conjg(U')
      call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulat(:,:,lspl),lmmax,h,lmmax, &
       zzero,a,lmmax)
      call zgemm('N','C',lmmax,lmmax,lmmax,zone,a,lmmax,ulat(:,:,lspl),lmmax, &
       zone,b,lmmax)
    end do
! block diagonalise symmetrised H
    do l=0,lmax
      n=2*l+1
      lm=idxlm(l,-l)
      call zheev('V','U',n,b(lm,lm),lmmax,elm(lm,ias),work,lwork,rwork,info)
    end do
! the unitary matrix U is the transpose of the eigenvector array
    do i=1,lmmax
      do j=1,lmmax
        ulm(i,j,ias)=b(j,i)
      end do
    end do
  end do
end do
deallocate(rwork,ulat,a,b,h,work)
return
end subroutine

