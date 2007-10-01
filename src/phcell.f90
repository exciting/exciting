
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phcell(iph,iq,is,ia,ip)
use modmain
implicit none
! arguments
integer, intent(in) :: iph
integer, intent(in) :: iq
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ip
! local variables
integer js,ja,i,j,n,iv(3)
integer i1,i2,i3,m(3,3)
real(8) a(3,3),ai(3,3),apl(3,maxatoms)
real(8) v(3),vb(3,maxatoms),dmin,t1
! external functions
real(8) r3taxi
external r3taxi
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
        v(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
        t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
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
          v(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
          t1=v(1)**2+v(2)**2+v(3)**2
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
          v(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
          t1=v(1)**2+v(2)**2+v(3)**2
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
  a(:,i)=dble(m(1,i))*avec(:,1)+dble(m(2,i))*avec(:,2)+dble(m(3,i))*avec(:,3)
end do
! inverse of lattice vector matrix
call r3minv(a,ai)
! offset vectors for the atomic basis
n=1
vb(:,1)=0.d0
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      if (n.eq.nphcell) goto 30
      v(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
      call r3mv(ai,v,v)
      call r3frac(epslat,v,iv)
      call r3mv(a,v,v)
      do i=1,n
        if (r3taxi(v,vb(1,i)).lt.epslat) goto 20
      end do
      n=n+1
      vb(:,n)=v(:)
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
  i=0
  do j=1,nphcell
    do ja=1,natoms(js)
      i=i+1
      if (i.gt.maxatoms) then
        write(*,*)
        write(*,'("Error(phcell): too many atoms in supercell : ",I8)') i
        write(*,'(" for species ",I4)') js
        write(*,'("Adjust maxatoms in modmain and recompile code")')
        write(*,*)
        stop
      end if
      v(:)=vb(:,j)+atposc(:,ja,js)
! add small periodic displacement
      if ((is.eq.js).and.(ia.eq.ja)) then
        t1=vb(1,j)*vqc(1,iq)+vb(2,j)*vqc(2,iq)+vb(3,j)*vqc(3,iq)
        if (iph.eq.0) then
          v(ip)=v(ip)+deltaph*cos(t1)
        else
          v(ip)=v(ip)+deltaph*sin(t1)
        end if
      end if
! convert to new lattice coordinates
      call r3mv(ai,v,apl(1,i))
      call r3frac(epslat,apl(1,i),iv)
    end do
  end do
  atposl(:,1:i,js)=apl(:,1:i)
  natoms(js)=i
end do
avec(:,:)=a(:,:)
! muffin-tin magnetic fields should be zero
bfcmt(:,:,:)=0.d0
return
end subroutine

