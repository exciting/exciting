
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genstark
! !INTERFACE:
subroutine genstark
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Generates the stars for the ${\bf k}$-point set as reference to crystal
!   symmetries. For a non-zero ${\bf q}$-point the little group of ${\bf q}$
!   is taken instead of the full symmetry group.
!
! !REVISION HISTORY:
!   Created December 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  integer :: i,isym,lspl,ik,iknr,iv(3),nsymcr
  real(8) :: t1,v1(3),v2(3),s(3,3)
  real(8), external :: r3taxi
  ! start with all crystal symmetries
  nsymcr=nsymcrys
  ! consider small/little group of q
  if (iqcu.ne.0) nsymcr=nsymcrysq(iqcu)
  nsymcrysstr(:)=0
  do iknr=1,nkptnr
     iv(:)=ivknr(:,iknr)
     ik=ikmap(iv(1),iv(2),iv(3))
     v1(:)=vkl(:,ik)
     do i=1,nsymcr
        isym=i
        ! consider small/little group of q
        if (iqcu.ne.0) isym=scqmap(i,iqcu)
        lspl=lsplsymc(isym)
        s(:,:)=dble(symlat(:,:,lspl))
        call r3mtv(s,v1,v2)
        call r3frac(epslat,v2,iv)
        t1=r3taxi(vklnr(1,iknr),v2)
        if (t1.lt.epslat) then
           ! symmetry element is in star
           nsymcrysstr(ik)=nsymcrysstr(ik)+1
           scmapstr(nsymcrysstr(ik),ik)=isym
           ikstrmapiknr(nsymcrysstr(ik),ik)=iknr
           strmap(iknr)=ik
           strmapsymc(iknr)=isym
           goto 10
        end if   
     end do
10   continue
  end do
  ! debug output
  if (dbglev.gt.1) then
     write(*,'(a)') 'Debug(genstark):'
     write(*,'(a)') 'ik,i,vkl,isym,lspl,vklnr,iknr'
     do ik=1,nkpt
        v1(:)=vkl(:,ik)
        do i=1,nsymcrysstr(ik)
           isym=scmapstr(i,ik)
           lspl=lsplsymc(isym)
           s(:,:)=dble(symlat(:,:,lspl))
           call r3mtv(s,v1,v2)
           call r3frac(epslat,v2,iv)
           iv(:)=int(v2(:)*dble(ngridk(:)))
           iknr=ikmapnr(iv(1),iv(2),iv(3))
           write(*,'(2i6,3f12.4,2i6,3f12.4,i6)') ik,i,v1,isym,lspl,v2,iknr
        end do
     end do
     write(*,*)
     write(*,'(a)') ' ik,wkpt*nkptnr,nsymcrysstr,scmapstr:'
     do ik=1,nkpt
        write(*,'(i9,f12.4,i9,3x,192i4)') ik,wkpt(ik)*nkptnr,nsymcrysstr(ik),&
             scmapstr(1:nsymcrysstr(ik),ik)
     end do
     write(*,*)
     write(*,'(a)') ' ik,wkpt*nkptnr,nsymcrysstr,lsplmapstr:'
     do ik=1,nkpt
        write(*,'(i9,f12.4,i9,3x,192i4)') ik,wkpt(ik)*nkptnr,nsymcrysstr(ik),&
             lsplsymc(scmapstr(1:nsymcrysstr(ik),ik))
     end do
     write(*,*)
     write(*,'(a)') ' ik,wkpt*nkptnr,nsymcrysstr,ikstrmapiknr:'
     do ik=1,nkpt
        write(*,'(i9,f12.4,i9,3x,192i4)') ik,wkpt(ik)*nkptnr,nsymcrysstr(ik),&
             ikstrmapiknr(1:nsymcrysstr(ik),ik)
     end do
     write(*,*)
  end if
end subroutine genstark
!EOC
