
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phcell(iph,dph,iq,is,ia,ip)
use modmain
implicit none
! arguments
integer, intent(in) :: iph
real(8), intent(in) :: dph
integer, intent(in) :: iq
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ip
! local variables
integer js,ja,na,i,n,iv(3)
integer i1,i2,i3,m(3,3)
real(8) v1(3),v2(3),dmin,t1
if ((iph.ne.0).and.(iph.ne.1)) then
  write(*,*)
  write(*,'("Error(phcell): iph should be 0 or 1 : ",I8)') iph
  write(*,*)
  stop
end if
! check for Gamma-point phonon
if ((ivq(1,iq).eq.0).and.(ivq(2,iq).eq.0).and.(ivq(3,iq).eq.0)) then
  m(:,:)=0
  m(1,1)=1
  m(2,2)=1
  m(3,3)=1
  nphcell=1
  goto 10
end if
! find the first lattice vector
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iq)+dble(i2)*vql(2,iq)+dble(i3)*vql(3,iq)
      if (abs(t1-nint(t1)).lt.epslat) then
        v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
        t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
        if ((t1.lt.dmin).and.(t1.gt.epslat)) then
          m(1,1)=i1
          m(2,1)=i2
          m(3,1)=i3
          dmin=t1
        end if
      end if
    end do
  end do
end do
! find the second lattice vector
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iq)+dble(i2)*vql(2,iq)+dble(i3)*vql(3,iq)
      if (abs(t1-nint(t1)).lt.epslat) then
! area defined by first two lattice vectors
        n=(i2*m(3,1)-i3*m(2,1))**2 &
         +(i3*m(1,1)-i1*m(3,1))**2 &
         +(i1*m(2,1)-i2*m(1,1))**2
        if (n.ne.0) then
          v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
          t1=v1(1)**2+v1(2)**2+v1(3)**2
          if (t1.lt.dmin) then
            m(1,2)=i1
            m(2,2)=i2
            m(3,2)=i3
            dmin=t1
          end if
        end if
      end if
    end do
  end do
end do
! find the third lattice vector
nphcell=0
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iq)+dble(i2)*vql(2,iq)+dble(i3)*vql(3,iq)
      if (abs(t1-nint(t1)).lt.epslat) then
! number of primitive unit cells in supercell
        n=m(1,2)*(i2*m(3,1)-i3*m(2,1)) &
         +m(2,2)*(i3*m(1,1)-i1*m(3,1)) &
         +m(3,2)*(i1*m(2,1)-i2*m(1,1))
        if (n.ne.0) then
          v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
          t1=v1(1)**2+v1(2)**2+v1(3)**2
          if (t1.lt.dmin) then
            nphcell=abs(n)
            m(1,3)=i1
            m(2,3)=i2
            m(3,3)=i3
            dmin=t1
          end if
        end if
      end if
    end do
  end do
end do
if (nphcell.eq.0) then
  write(*,*)
  write(*,'("Error(phcell): unable to generate supercell")')
  write(*,*)
  stop
end if
10 continue
! new lattice vectors
do i=1,3
  avec(:,i)=dble(m(1,i))*avec0(:,1)+dble(m(2,i))*avec0(:,2) &
   +dble(m(3,i))*avec0(:,3)
end do
! inverse of lattice vector matrix
call r3minv(avec,ainv)
! generate offset vectors for each primitive cell in the supercell
n=1
vphcell(:,1)=0.d0
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      if (n.eq.nphcell) goto 30
      v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
      call r3mv(ainv,v1,v2)
      call r3frac(epslat,v2,iv)
      call r3mv(avec,v2,v1)
      do i=1,n
        t1=abs(v1(1)-vphcell(1,i))+abs(v1(2)-vphcell(2,i)) &
         +abs(v1(3)-vphcell(3,i))
        if (t1.lt.epslat) goto 20
      end do
      n=n+1
      vphcell(:,n)=v1(:)
20 continue
    end do
  end do
end do
write(*,*)
write(*,'("Error(phcell): unable to generate supercell")')
write(*,*)
stop
30 continue
! set up new atomic positions
do js=1,nspecies
  na=0
  do ja=1,natoms(js)
    do i=1,nphcell
      na=na+1
      if (na.gt.maxatoms) then
        write(*,*)
        write(*,'("Error(phcell): too many atoms in supercell : ",I8)') na
        write(*,'(" for species ",I4)') js
        write(*,'("Adjust maxatoms in modmain and recompile code")')
        write(*,*)
        stop
      end if
      v1(:)=vphcell(:,i)+atposc0(:,ja,js)
! add small periodic displacement
      if ((is.eq.js).and.(ia.eq.ja)) then
        t1=dot_product(vqc(:,iq),vphcell(:,i))
        if (iph.eq.0) then
          v1(ip)=v1(ip)+dph*cos(t1)
        else
          v1(ip)=v1(ip)+dph*sin(t1)
        end if
      end if
! convert to new lattice coordinates
      call r3mv(ainv,v1,atposl(:,na,js))
      call r3frac(epslat,atposl(:,na,js),iv)
    end do
  end do
  natoms(js)=na
end do
! muffin-tin magnetic fields should be zero
bfcmt(:,:,:)=0.d0
return
end subroutine

