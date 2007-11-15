
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstark
  use modmain
  use modxs
  implicit none
  ! arguments
  ! local variables


  ! initialize stars
  nsymcrysstr(:)=0
  scmapstr(:,:)=1

!!$  ! INSIDE LOOP: add symmetry operation to star
!!$  scmapstr(nint(wppt(jp)/t1),jp)=scimap(isym)

  ! number of elements in stars
  nsymcrysstr(:)=nint(wkpt(:)*dble(ngridk(1)*ngridk(2)*ngridk(3)))
  nsymcrysstrmax=maxval(nsymcrysstr)

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
