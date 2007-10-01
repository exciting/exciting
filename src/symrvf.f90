
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrvf(tpolar,lrstp,nd,rvfmt,rvfir)
use modmain
implicit none
! arguments
logical, intent(in) :: tpolar
integer, intent(in) :: lrstp
integer, intent(in) :: nd
real(8), intent(inout) :: rvfmt(lmmaxvr,nrmtmax,natmtot,nd)
real(8), intent(inout) :: rvfir(ngrtot,nd)
! local variables
integer is,ia1,ia2,ias1,ias2
integer ir,i,n,md,isym,sym(3,3)
real(8) t1
! allocatable arrays
logical, allocatable :: done(:)
real(8), allocatable :: rvfmt1(:,:,:,:)
real(8), allocatable :: rvfmt2(:,:,:)
real(8), allocatable :: rvfmt3(:,:,:)
real(8), allocatable :: rvfir1(:,:)
real(8), allocatable :: rvfir2(:,:)
! external functions
integer i3mdet
external i3mdet
allocate(done(natmmax))
allocate(rvfmt1(lmmaxvr,nrmtmax,3,natmmax))
allocate(rvfmt2(lmmaxvr,nrmtmax,3))
allocate(rvfmt3(lmmaxvr,nrmtmax,3))
allocate(rvfir1(ngrtot,3))
allocate(rvfir2(ngrtot,3))
if ((nd.ne.1).and.(nd.ne.3)) then
  write(*,*)
  write(*,'("Error(symrvf): nd should be 1 or 3 : ",I8)') nd
  write(*,*)
  stop
end if
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
  do ia1=1,natoms(is)
    ias1=idxas(ia1,is)
    if (nd.eq.3) then
      do ir=1,nrmt(is),lrstp
        rvfmt1(:,ir,:,ia1)=rvfmt(:,ir,ias1,:)
      end do
    else
      do ir=1,nrmt(is),lrstp
        rvfmt1(:,ir,1:2,ia1)=0.d0
        rvfmt1(:,ir,3,ia1)=rvfmt(:,ir,ias1,1)
      end do
    end if
  end do
  done(:)=.false.
! rotate equivalent atoms and add
  do ia1=1,natoms(is)
    if (.not.done(ia1)) then
      ias1=idxas(ia1,is)
      do ir=1,nrmt(is),lrstp
        rvfmt(:,ir,ias1,:)=0.d0
      end do
      n=0
      do ia2=1,natoms(is)
        do i=1,nsymeqat(ia2,ia1,is)
          isym=symeqat(i,ia2,ia1,is)
          if (tpolar) then
            md=1
          else
            md=i3mdet(symlat(1,1,isym))
          end if
          sym(:,:)=md*symlat(:,:,isym)
          call symrvfmt(lrstp,is,sym,rvfmt1(1,1,1,ia2),rvfmt2)
          if (nd.eq.3) then
            do ir=1,nrmt(is),lrstp
              rvfmt(:,ir,ias1,:)=rvfmt(:,ir,ias1,:)+rvfmt2(:,ir,:)
            end do
          else
            do ir=1,nrmt(is),lrstp
              rvfmt(:,ir,ias1,1)=rvfmt(:,ir,ias1,1)+rvfmt2(:,ir,3)
            end do
          end if
          n=n+1
        end do
      end do
      t1=1.d0/dble(n)
      do ir=1,nrmt(is),lrstp
        rvfmt(:,ir,ias1,:)=t1*rvfmt(:,ir,ias1,:)
      end do
      done(ia1)=.true.
! rotate and copy to equivalent atoms
      do ia2=1,natoms(is)
        if (.not.done(ia2)) then
          if (nsymeqat(ia1,ia2,is).gt.0) then
            ias2=idxas(ia2,is)
            isym=symeqat(1,ia1,ia2,is)
            if (tpolar) then
              md=1
            else
              md=i3mdet(symlat(1,1,isym))
            end if
            sym(:,:)=md*symlat(:,:,isym)
            if (nd.eq.3) then
              do ir=1,nrmt(is),lrstp
                rvfmt2(:,ir,:)=rvfmt(:,ir,ias1,:)
              end do
            else
              do ir=1,nrmt(is),lrstp
                rvfmt2(:,ir,1:2)=0.d0
                rvfmt2(:,ir,3)=rvfmt(:,ir,ias1,1)
              end do
            end if
            call symrvfmt(lrstp,is,sym,rvfmt2,rvfmt3)
            if (nd.eq.3) then
              do ir=1,nrmt(is),lrstp
                rvfmt(:,ir,ias2,:)=rvfmt3(:,ir,:)
              end do
            else
              do ir=1,nrmt(is),lrstp
                rvfmt(:,ir,ias2,1)=rvfmt3(:,ir,3)
              end do
            end if
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
if (nd.eq.3) then
  rvfir1(:,:)=rvfir(:,:)
else
  rvfir1(:,1:2)=0.d0
  rvfir1(:,3)=rvfir(:,1)
end if
rvfir(:,:)=0.d0
do isym=1,nsymcrys
  if (tpolar) then
    md=1
  else
    md=i3mdet(symcrys(1,1,isym))
  end if
  sym(:,:)=md*symcrys(:,:,isym)
  call symrvfir(sym,rvfir1,rvfir2)
  if (nd.eq.3) then
    rvfir(:,:)=rvfir(:,:)+rvfir2(:,:)
  else
    rvfir(:,1)=rvfir(:,1)+rvfir2(:,3)
  end if
end do
t1=1.d0/dble(nsymcrys)
rvfir(:,:)=t1*rvfir(:,:)
deallocate(done,rvfmt1,rvfmt2,rvfmt3,rvfir1,rvfir2)
return
end subroutine
!EOC
