
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrf(lrstp,rfmt,rfir)
use modmain
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(inout) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(inout) :: rfir(ngrtot)
! local variables
integer is,ia1,ia2,ias1,ias2
integer ir,n,i,isym
real(8) t1
! allocatable arrays
logical, allocatable :: done(:)
real(8), allocatable :: rfmt1(:,:,:)
real(8), allocatable :: rfmt2(:,:)
real(8), allocatable :: rfir1(:)
real(8), allocatable :: rfir2(:)
allocate(done(natmmax))
allocate(rfmt1(lmmaxvr,nrmtmax,natmmax))
allocate(rfmt2(lmmaxvr,nrmtmax))
allocate(rfir1(ngrtot))
allocate(rfir2(ngrtot))
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
  do ia1=1,natoms(is)
    ias1=idxas(ia1,is)
    do ir=1,nrmt(is),lrstp
      rfmt1(:,ir,ia1)=rfmt(:,ir,ias1)
    end do
  end do
  done(:)=.false.
! rotate equivalent atoms and add
  do ia1=1,natoms(is)
    if (.not.done(ia1)) then
      ias1=idxas(ia1,is)
      do ir=1,nrmt(is),lrstp
        rfmt(:,ir,ias1)=0.d0
      end do
      n=0
      do ia2=1,natoms(is)
        do i=1,nsymeqat(ia2,ia1,is)
          isym=symeqat(i,ia2,ia1,is)
          call symrfmt(lrstp,is,symlat(1,1,isym),rfmt1(1,1,ia2),rfmt2)
          do ir=1,nrmt(is),lrstp
            rfmt(:,ir,ias1)=rfmt(:,ir,ias1)+rfmt2(:,ir)
          end do
          n=n+1
        end do
      end do
      t1=1.d0/dble(n)
      do ir=1,nrmt(is),lrstp
        rfmt(:,ir,ias1)=t1*rfmt(:,ir,ias1)
      end do
      done(ia1)=.true.
! rotate and copy to equivalent atoms
      do ia2=1,natoms(is)
        if (.not.done(ia2)) then
          if (nsymeqat(ia1,ia2,is).gt.0) then
            ias2=idxas(ia2,is)
            isym=symeqat(1,ia1,ia2,is)
            call symrfmt(lrstp,is,symlat(1,1,isym),rfmt(1,1,ias1), &
             rfmt(1,1,ias2))
            done(ia2)=.true.
          end if
        end if
      end do
    end if
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
rfir1(:)=0.d0
do isym=1,nsymcrys
  call symrfir(symcrys(1,1,isym),rfir,rfir2)
  rfir1(:)=rfir1(:)+rfir2(:)
end do
t1=1.d0/dble(nsymcrys)
rfir(:)=t1*rfir1(:)
deallocate(done,rfmt1,rfmt2,rfir1,rfir2)
return
end subroutine
