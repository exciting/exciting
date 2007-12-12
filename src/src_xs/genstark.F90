
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstark
  use modmain
  use modxs
  implicit none
  ! local variables
  integer :: isym,ik,iknr,iv(3)
  real(8) :: v1(3),v2(3)

  ! number of elements in stars
  nsymcrysstr(:)=nint(wkpt(:)*dble(nkptnr))
  nsymcrysstrmax=maxval(nsymcrysstr)

!!$  do ik=1,nkpt
  ! INSIDE LOOP: add symmetry operation to star
!!$  scmapstr(nint(wppt(ik)/t1),ik)=scimap(isym)
!!$  end do

  do iknr=1,nkptnr
     iv(:)=ivknr(:,iknr)
     ik=ikmap(iv(1),iv(2),iv(3))
     v1(:)=vkl(:,ik)
     do isym=1,nsymcrysq(iqcu)
        lspl=lsplsymc(isym)
        s(:,:)=dble(symlat(:,:,lspl))
        call r3mtv(s,v1,v2)
        call r3frac(epslat,v2,iv)
        t2=r3taxi(vpl(1,jp),v2)
        if (t2.lt.epslat) then
              
     end do
  end do

!!$! *** TEST ***
!!$i1=0
!!$write(*,*) 'writing out: genppts'
!!$write(*,*)
!!$do ip=1,nppt
!!$   do jsym=1,nsymcrysstr(ip)
!!$      i1=i1+1
!!$!      write(*,*) 'TEST: ip,jstar,counter',ip,jsym,i1
!!$      lspl=lsplsymc(scmapstr(jsym,ip))
!!$      v1(:)=vpl(:,ip)
!!$      s(:,:)=dble(symlat(:,:,lspl))
!!$      call r3mtv(s,v1,v2)
!!$      call r3frac(epslat,v2,iv)
!!$      write(*,'(a,4i6,6f12.2)') 'c,ip,jsym,lspl,vpl,s.vpl',&
!!$           i1,ip,jsym,lspl,v1,v2
!!$      write(10,*) v2
!!$   end do
!!$end do

end subroutine genstark
