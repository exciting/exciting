
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvf
! !INTERFACE:
subroutine symrvf(lrstp,rvfmt,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rvfmt : real muffin-tin vector field
!           (in,real(lmmaxvr,nrmtmax,natmtot,ndmag))
!   rvfir : real interstitial vector field
!           (in,real(ngrtot,ndmag))
! !DESCRIPTION:
!   Symmetrises a vector field defined over the entire unit cell using the full
!   set of crystal symmetries. If a particular symmetry involves rotating atom
!   1 into atom 2, then the spatial and spin rotations of that symmetry are
!   applied to the vector field in atom 2 (expressed in spherical harmonic
!   coefficients), which is then added to the field in atom 1. This is repeated
!   for all symmetry operations. The fully symmetrised field in atom 1 is then
!   rotated and copied to atom 2. Symmetrisation of the interstitial part of the
!   field is performed by {\tt symrvfir}. See also {\tt symrfmt} and
!   {\tt findsym}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!   Fixed problem with improper rotations, February 2008 (L. Nordstrom,
!    F. Bultmark and F. Cricchio)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(inout) :: rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag)
real(8), intent(inout) :: rvfir(ngrtot,ndmag)
! local variables
integer is,ia,ja,ias,jas
integer isym,ir,lm,i,md
integer lspl,ilspl,lspn,ilspn
real(8) sc(3,3),v(3),t1
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
          call symrfmt(lrstp,is,symlatc(1,1,lspl),rvfmt1(1,1,ja,i), &
           rvfmt2(1,1,i))
        end do
! global spin proper rotation matrix in Cartesian coordinates
        lspn=lspnsymc(isym)
        md=symlatd(lspn)
        sc(:,:)=dble(md)*symlatc(:,:,lspn)
! global spin rotation of vector field
        if (ndmag.eq.3) then
! non-collinear case
          do ir=1,nrmt(is),lrstp
            do lm=1,lmmaxvr
              v(1)=sc(1,1)*rvfmt2(lm,ir,1) &
                  +sc(1,2)*rvfmt2(lm,ir,2) &
                  +sc(1,3)*rvfmt2(lm,ir,3)
              v(2)=sc(2,1)*rvfmt2(lm,ir,1) &
                  +sc(2,2)*rvfmt2(lm,ir,2) &
                  +sc(2,3)*rvfmt2(lm,ir,3)
              v(3)=sc(3,1)*rvfmt2(lm,ir,1) &
                  +sc(3,2)*rvfmt2(lm,ir,2) &
                  +sc(3,3)*rvfmt2(lm,ir,3)
              rvfmt(lm,ir,ias,:)=rvfmt(lm,ir,ias,:)+v(:)
            end do
          end do
        else
! collinear case
          do ir=1,nrmt(is),lrstp
            do lm=1,lmmaxvr
              rvfmt(lm,ir,ias,1)=rvfmt(lm,ir,ias,1)+sc(3,3)*rvfmt2(lm,ir,1)
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
          ilspl=isymlat(lspl)
          do i=1,ndmag
            call symrfmt(lrstp,is,symlatc(1,1,ilspl),rvfmt(1,1,ias,i), &
             rvfmt(1,1,jas,i))
          end do
! inverse of global proper rotation matrix in Cartesian coordinates
          lspn=lspnsymc(isym)
          ilspn=isymlat(lspn)
          md=symlatd(ilspn)
          sc(:,:)=dble(md)*symlatc(:,:,ilspn)
! global spin rotation of vector field
          if (ndmag.eq.3) then
! non-collinear case
            do ir=1,nrmt(is),lrstp
              do lm=1,lmmaxvr
                v(:)=rvfmt(lm,ir,jas,:)
                rvfmt(lm,ir,jas,1)=sc(1,1)*v(1)+sc(1,2)*v(2)+sc(1,3)*v(3)
                rvfmt(lm,ir,jas,2)=sc(2,1)*v(1)+sc(2,2)*v(2)+sc(2,3)*v(3)
                rvfmt(lm,ir,jas,3)=sc(3,1)*v(1)+sc(3,2)*v(2)+sc(3,3)*v(3)
              end do
            end do
          else
! collinear case
            do ir=1,nrmt(is),lrstp
              do lm=1,lmmaxvr
                rvfmt(lm,ir,jas,1)=sc(3,3)*rvfmt(lm,ir,jas,1)
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
!EOC

