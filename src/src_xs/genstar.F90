
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genstar
! !INTERFACE:
subroutine genstar(tsmgrq,ngridp,nppt,npptnr,vpl,vplnr,ivpnr,ipmap,ipmapnr, &
     nsymcrysst,scmapst,ipstmapipnr,stmap,stmapsymc)
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Generates the stars for the $k$-point set as reference to crystal
!   symmetries. For a non-zero q-point the little group of q is taken
!   instead of the full symmetry group.
!
! !REVISION HISTORY:
!   Created December 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  logical, intent(in) :: tsmgrq
  integer, intent(in) :: ngridp(3)
  integer, intent(in) :: nppt
  integer, intent(in) :: npptnr
  real(8), intent(in) :: vpl(3,nppt)
  real(8), intent(in) :: vplnr(3,npptnr)
  integer, intent(in) :: ivpnr(3,npptnr)
  integer, intent(in) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
  integer, intent(in) :: ipmapnr(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
  integer, intent(out) :: nsymcrysst(nppt)
  integer, intent(out) :: scmapst(nsymcrys,nppt)
  integer, intent(out) :: ipstmapipnr(maxsymcrys,nppt)
  integer, intent(out) :: stmap(npptnr)
  integer, intent(out) :: stmapsymc(npptnr)
  ! local variables
  integer :: i,isym,lspl,ip,ipnr,iv(3),nsymcr
  real(8) :: t1,v1(3),v2(3),s(3,3)
  real(8), external :: r3taxi
  ! start with all crystal symmetries
  nsymcr=nsymcrys
  ! consider small/little group of q
  if (tsmgrq.and.(iqcu.ne.0)) nsymcr=nsymcrysq(iqcu)
  nsymcrysst(:)=0
  do ipnr=1,npptnr
     iv(:)=ivpnr(:,ipnr)
     ip=ipmap(iv(1),iv(2),iv(3))
     v1(:)=vpl(:,ip)
     do i=1,nsymcr
        isym=i
        ! consider small/little group of q
        if (iqcu.ne.0) isym=scqmap(i,iqcu)
        lspl=lsplsymc(isym)
        s(:,:)=dble(symlat(:,:,lspl))
        call r3mtv(s,v1,v2)
        call r3frac(epslat,v2,iv)
        t1=r3taxi(vplnr(1,ipnr),v2)
        if (t1.lt.epslat) then
           ! symmetry element is in star
           nsymcrysst(ip)=nsymcrysst(ip)+1
           scmapst(nsymcrysst(ip),ip)=isym
           ipstmapipnr(nsymcrysst(ip),ip)=ipnr
           stmap(ipnr)=ip
           stmapsymc(ipnr)=isym
           goto 10
        end if   
     end do
10   continue
  end do
  ! debug output
  if (dbglev.gt.1) then
     write(*,'(a)') 'Debug(genstark):'
     write(*,'(a)') 'ip,i,vpl,isym,lspl,vplnr,ipnr'
     do ip=1,nppt
        v1(:)=vpl(:,ip)
        do i=1,nsymcrysst(ip)
           isym=scmapst(i,ip)
           lspl=lsplsymc(isym)
           s(:,:)=dble(symlat(:,:,lspl))
           call r3mtv(s,v1,v2)
           call r3frac(epslat,v2,iv)
           iv(:)=int(v2(:)*dble(ngridp(:)))
           ipnr=ipmapnr(iv(1),iv(2),iv(3))
           write(*,'(2i6,3f12.4,2i6,3f12.4,i6)') ip,i,v1,isym,lspl,v2,ipnr
        end do
     end do
     write(*,*)
     write(*,'(a)') ' ip,wkpt*npptnr,nsymcrysst,scmapst:'
     do ip=1,nppt
        write(*,'(i9,f12.4,i9,3x,192i4)') ip,wkpt(ip)*npptnr,nsymcrysst(ip),&
             scmapst(1:nsymcrysst(ip),ip)
     end do
     write(*,*)
     write(*,'(a)') ' ip,wkpt*npptnr,nsymcrysst,lsplmapst:'
     do ip=1,nppt
        write(*,'(i9,f12.4,i9,3x,192i4)') ip,wkpt(ip)*npptnr,nsymcrysst(ip),&
             lsplsymc(scmapst(1:nsymcrysst(ip),ip))
     end do
     write(*,*)
     write(*,'(a)') ' ip,wkpt*npptnr,nsymcrysst,ipstmapipnr:'
     do ip=1,nppt
        write(*,'(i9,f12.4,i9,3x,192i4)') ip,wkpt(ip)*npptnr,nsymcrysst(ip),&
             ipstmapipnr(1:nsymcrysst(ip),ip)
     end do
     write(*,*)
  end if
end subroutine genstar
!EOC
