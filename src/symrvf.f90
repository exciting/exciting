
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrvf(lrstp,rvfmt,rvfir)
use modmain
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(inout) :: rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag)
real(8), intent(inout) :: rvfir(ngrtot,ndmag)
! local variables
integer is,ia,ja,ias,jas
integer isym,lspl,lspn
integer ir,lm,i,sym(3,3)
real(8) s(3,3),v(3),t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rvfmt1(:,:,:,:),rvfmt2(:,:,:)
allocate(rvfmt1(lmmaxvr,nrmtmax,natmmax,ndmag))
allocate(rvfmt2(lmmaxvr,nrmtmax,ndmag))
t1=1.d0/dble(nsymcrys)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
! make copy of vector field for all atoms of current species
  do i=1,ndmag
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        rvfmt1(:,ir,ia,i)=rvfmt(:,ir,ias,i)
      end do
    end do
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        rvfmt(:,ir,ias,:)=0.d0
      end do
! begin loop over crystal symmetries
      do isym=1,nsymcrys
! equivalent atom
        ja=ieqatom(ia,is,isym)
! parallel transport of vector field
        lspl=lsplsymc(isym)
        do i=1,ndmag
          call symrfmt(lrstp,is,symlat(1,1,lspl),rvfmt1(1,1,ja,i),rvfmt2(1,1,i))
        end do
! global spin rotation matrix in Cartesian coordinates
        lspn=lspnsymc(isym)
        s(:,:)=dble(symlat(:,:,lspn))
        call r3mm(s,ainv,s)
        call r3mm(avec,s,s)
! global spin rotation of vector field
        if (ndmag.eq.3) then
! non-collinear case
          do ir=1,nrmt(is),lrstp
            do lm=1,lmmaxvr
              v(1)=s(1,1)*rvfmt2(lm,ir,1) &
                  +s(1,2)*rvfmt2(lm,ir,2) &
                  +s(1,3)*rvfmt2(lm,ir,3)
              v(2)=s(2,1)*rvfmt2(lm,ir,1) &
                  +s(2,2)*rvfmt2(lm,ir,2) &
                  +s(2,3)*rvfmt2(lm,ir,3)
              v(3)=s(3,1)*rvfmt2(lm,ir,1) &
                  +s(3,2)*rvfmt2(lm,ir,2) &
                  +s(3,3)*rvfmt2(lm,ir,3)
              rvfmt(lm,ir,ias,:)=rvfmt(lm,ir,ias,:)+v(:)
            end do
          end do
        else
! collinear case
          do ir=1,nrmt(is),lrstp
            do lm=1,lmmaxvr
              rvfmt(lm,ir,ias,1)=rvfmt(lm,ir,ias,1)+s(3,3)*rvfmt2(lm,ir,1)
            end do
          end do
        end if
! end loop over crystal symmetries
      end do
! normalise
      do ir=1,nrmt(is),lrstp
        rvfmt(:,ir,ias,:)=t1*rvfmt(:,ir,ias,:)
      end do
! mark atom as done
      done(ia)=.true.
! rotate into equivalent atoms
      do isym=1,nsymcrys
        ja=ieqatom(ia,is,isym)
        if (.not.done(ja)) then
          jas=idxas(ja,is)
! parallel transport of vector field (using operation inverse)
          lspl=lsplsymc(isym)
          call i3minv(symlat(1,1,lspl),sym)
          do i=1,ndmag
            call symrfmt(lrstp,is,sym,rvfmt(1,1,ias,i),rvfmt(1,1,jas,i))
          end do
! inverse of global rotation matrix in Cartesian coordinates
          lspn=lspnsymc(isym)
          call i3minv(symlat(1,1,lspn),sym)
          s(:,:)=dble(sym(:,:))
          call r3mm(s,ainv,s)
          call r3mm(avec,s,s)
! global spin rotation of vector field
          if (ndmag.eq.3) then
! non-collinear case
            do ir=1,nrmt(is),lrstp
              do lm=1,lmmaxvr
                rvfmt(lm,ir,jas,1)=s(1,1)*rvfmt(lm,ir,ias,1) &
                                  +s(1,2)*rvfmt(lm,ir,ias,2) &
                                  +s(1,3)*rvfmt(lm,ir,ias,3)
                rvfmt(lm,ir,jas,2)=s(2,1)*rvfmt(lm,ir,ias,1) &
                                  +s(2,2)*rvfmt(lm,ir,ias,2) &
                                  +s(2,3)*rvfmt(lm,ir,ias,3)
                rvfmt(lm,ir,jas,3)=s(3,1)*rvfmt(lm,ir,ias,1) &
                                  +s(3,2)*rvfmt(lm,ir,ias,2) &
                                  +s(3,3)*rvfmt(lm,ir,ias,3)
              end do
            end do
          else
! collinear case
            do ir=1,nrmt(is),lrstp
              do lm=1,lmmaxvr
                rvfmt(lm,ir,jas,1)=s(3,3)*rvfmt(lm,ir,ias,1)
              end do
            end do
          end if
! mark atom as done
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
call symrvfir(ngvec,rvfir)
deallocate(rvfmt1,rvfmt2)
return
end subroutine

