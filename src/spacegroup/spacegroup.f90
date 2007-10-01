
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program spacegroup
! The authors of the LMGP Suite of programs for the interpretation of X-ray
! Experiments, Jean laugier and Bernard Bochu, are acknowleged for the use of
! their spacegroup data.
implicit none
! local variables
logical primcell
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=1000
! maximum allowed Wyckoff positions
integer, parameter :: maxpos=100
! maximum number of symmetry operations
integer, parameter :: maxop=192
integer nspecies,ncell(3),natmtot
integer is,ia,ip,id(3),ict
integer i1,i2,i3,i,j,nop,iostat
real(8), parameter :: pi=3.1415926535897932385d0
! zero vector tolerance
real(8), parameter :: eps=1.d-6
! conversion from atomic units to Angstroms
real(8), parameter :: au2ang=0.5291772108d0
real(8) a,b,c,ab,ac,bc,abr,acr,bcr
real(8) sab,cab,cac,cbc
real(8) avec(3,3),v1(3),v2(3)
! Hall symbol
character(25) hall
! input string from spacegroup.dat
character(25) str
! crystal type
character(12) ctype
! species symbols
character(3) spsymb(maxspecies)
! number of Wyckoff positions for each species
integer npos(maxspecies)
! number of atoms for each species
integer natoms(maxspecies)
! Wyckoff positions
real(8) pos(3,maxpos,maxspecies)
! atomic positions
real(8) atposl(3,maxatoms,maxspecies)
! magnetic fields
real(8) bfcmt(3,maxatoms,maxspecies)
! species filenames
character(256) spfname(maxspecies)
! symmetry operations
character(25) op(maxop)
! external functions
real(8) r3taxi
external r3taxi
! read in parameters from spacegroup.in
open(50,file='spacegroup.in',action='READ',status='OLD',form='FORMATTED')
read(50,*) hall
hall=adjustl(hall)
read(50,*) a,b,c
read(50,*) ab,ac,bc
read(50,*) ncell
if ((ncell(1).lt.1).or.(ncell(2).lt.1).or.(ncell(3).lt.1)) then
  write(*,*)
  write(*,'("Error(spacegroup): invalid ncell : ",3I8)') ncell
  write(*,*)
  stop
end if
read(50,*) primcell
read(50,*) nspecies
if (nspecies.le.0) then
  write(*,*)
  write(*,'("Error(spacegroup): nspecies <= 0 : ",I8)') nspecies
  write(*,*)
  stop
end if
if (nspecies.gt.maxspecies) then
  write(*,*)
  write(*,'("Error(spacegroup): nspecies too large : ",I8)') nspecies
  write(*,'("Adjust maxspecies and recompile code")')
  write(*,*)
  stop
end if
do is=1,nspecies
  read(50,*) spsymb(is),spfname(is)
  read(50,*) npos(is)
  if (npos(is).le.0) then
    write(*,*)
    write(*,'("Error(spacegroup): npos <=0 : ",I8)') npos(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (npos(is).gt.maxpos) then
    write(*,*)
    write(*,'("Error(spacegroup): npos too large : ",I8)') npos(is)
    write(*,'("Adjust maxpos and recompile code")')
    write(*,*)
    stop
  end if
  do ip=1,npos(is)
    read(50,*) pos(:,ip,is)
  end do
end do
close(50)
! read in the spacegroup data
open(50,file='spacegroup.dat',action='READ',status='OLD',form='FORMATTED')
10 continue
read(50,'(A25)',iostat=iostat) str
str=adjustl(str)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(spacegroup): cannot find Hall symbol : ",A)') trim(hall)
  write(*,*)
  stop
end if
if (trim(str).eq.'TRIC') then
  ict=1
else if (trim(str).eq.'MONO') then
  ict=2
else if (trim(str).eq.'ORTH') then
  ict=3
else if (trim(str).eq.'TETR') then
  ict=4
else if (trim(str).eq.'HEXA') then
  ict=5
else if (trim(str).eq.'RHOM') then
  ict=6
else if (trim(str).eq.'CUBI') then
  ict=7
else if (trim(str).eq.trim(hall)) then
  read(50,*) nop
  do j=1,nop
    read(50,'(A25)') op(j)
  end do
  goto 20
end if
goto 10
20 continue
close(50)
! crystal type
select case(ict)
case(1)
  ctype='triclinic'
case(2)
  ctype='monoclinic'
  bc=90.d0
  ac=90.d0
case(3)
  ctype='orthorhombic'
  ab=90.d0
  ac=90.d0
  bc=90.d0
case(4)
  ctype='tetragonal'
  if (a.ne.b) then
    write(*,*)
    write(*,'("Error(spacegroup): a <> b for tetragonal system :")')
    write(*,'(2G18.10)') a,b
    write(*,*)
    stop
  end if
  ab=90.d0
  ac=90.d0
  bc=90.d0
case(5)
  ctype='hexagonal'
  if (a.ne.b) then
    write(*,*)
    write(*,'("Error(spacegroup): a <> b for hexagonal system :")')
    write(*,'(2G18.10)') a,b
    write(*,*)
    stop
  end if
  ab=120.d0
  ac=90.d0
  bc=90.d0
case(6)
  ctype='rhombohedral'
  if ((a.ne.b).or.(a.ne.c).or.(b.ne.c)) then
    write(*,*)
    write(*,'("Error(spacegroup): lattice constants not equal for &
     &rhombohedral system :")')
    write(*,'(3G18.10)') a,b,c
    write(*,*)
    stop
  end if
  if ((ab.ne.ac).or.(ab.ne.bc).or.(ac.ne.bc)) then
    write(*,*)
    write(*,'("Error(spacegroup): angles not equal for rhombohedral &
     &system :")')
    write(*,'(3G18.10)') ab,ac,bc
    write(*,*)
    stop
  end if
case(7)
  ctype='cubic'
  if ((a.ne.b).or.(a.ne.c).or.(b.ne.c)) then
    write(*,*)
    write(*,'("Error(spacegroup): lattice constants not equal for cubic &
     &system :")')
    write(*,'(3G18.10)') a,b,c
    write(*,*)
    stop
  end if
  ab=90.d0
  ac=90.d0
  bc=90.d0
end select
! convert angles from degrees to radians
abr=ab*(pi/180.d0)
acr=ac*(pi/180.d0)
bcr=bc*(pi/180.d0)
! setup lattice vectors
sab=sin(abr)
if (abs(sab).lt.eps) then
  write(*,*)
  write(*,'("Error(spacegroup): degenerate lattice vectors")')
  write(*,*)
  stop
end if
cab=cos(abr)
cac=cos(acr)
cbc=cos(bcr)
avec(1,1)=a
avec(2,1)=0.d0
avec(3,1)=0.d0
avec(1,2)=b*cab
avec(2,2)=b*sab
avec(3,2)=0.d0
avec(1,3)=c*cac
avec(2,3)=c*(cbc-cab*cac)/sab
avec(3,3)=c*sqrt(sab**2-cac**2+2.d0*cab*cac*cbc-cbc**2)/sab
do i=1,3
  do j=1,3
    if (abs(avec(i,j)).lt.eps) avec(i,j)=0.d0
  end do
end do
! scale lattice vectors by the number of unit cells
do i=1,3
  avec(:,i)=avec(:,i)*dble(ncell(i))
end do
! compute the equivalent atomic positions
do is=1,nspecies
  natoms(is)=0
  do ip=1,npos(is)
    do j=1,nop
      call applyop(op(j),pos(1,ip,is),v1)
      do i1=0,ncell(1)-1
        do i2=0,ncell(2)-1
          do i3=0,ncell(3)-1
            v2(1)=(v1(1)+dble(i1))/dble(ncell(1))
            v2(2)=(v1(2)+dble(i2))/dble(ncell(2))
            v2(3)=(v1(3)+dble(i3))/dble(ncell(3))
            call r3frac(eps,v2,id)
            do ia=1,natoms(is)
              if (r3taxi(v2,atposl(1,ia,is)).lt.eps) goto 30
            end do
            natoms(is)=natoms(is)+1
            if (natoms(is).gt.maxatoms) then
              write(*,*)
              write(*,'("Error(spacegroup): natoms too large")')
              write(*,'(" for species ",I4)') is
              write(*,'("Adjust maxatoms and recompile code")')
              write(*,*)
              stop
            end if
            atposl(:,natoms(is),is)=v2(:)
          end do
        end do
      end do
30 continue
    end do
  end do
  natmtot=natmtot+natoms(is)
end do
! set magnetic fields to zero
bfcmt(:,:,:)=0.d0
! reduce conventional cell to primitive cell if required
if (primcell) call findprim(eps,avec,nspecies,natoms,maxatoms,atposl,bfcmt)
! find the total number of atoms
natmtot=0
do is=1,nspecies
  natmtot=natmtot+natoms(is)
end do
! write atomic positions to GEOMETRY.OUT
open(50,file='GEOMETRY.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("! Atomic positions generated by spacegroup program")')
write(50,'("!  crystal type : ",A12)') ctype
write(50,'("!  spacegroup : ",A)') trim(hall)
write(50,'("!  lattice constants (a,b,c) : ",3G18.10)') a,b,c
write(50,'("!  angles in degrees (ab,ac,bc) : ",3G18.10)') ab,ac,bc
write(50,'("!  number of conventional unit cells : ",3I4)') ncell
write(50,'("!  reduction to primitive cell : ",L1)') primcell
write(50,'("!  Wyckoff positions :")')
do is=1,nspecies
  write(50,'("!   species : ",I4,", ",A)') is,trim(spfname(is))
  do ip=1,npos(is)
    write(50,'("!   ",3G18.10)') pos(:,ip,is)
  end do
end do
write(50,*)
write(50,'("avec")')
write(50,'(3G18.10)') avec(:,1)
write(50,'(3G18.10)') avec(:,2)
write(50,'(3G18.10)') avec(:,3)
write(50,*)
write(50,'("atoms")')
write(50,'(I4,T40," : nspecies")') nspecies
do is=1,nspecies
  write(50,'("''",A,"''",T40," : spfname")') trim(spfname(is))
  write(50,'(I4,T40," : natoms; atposl, bfcmt below")') natoms(is)
  do ia=1,natoms(is)
    write(50,'(3F12.8,"  ",3F12.8)') atposl(:,ia,is),bfcmt(:,ia,is)
  end do
end do
close(50)
! write out the XCrySDen file
open(50,file='crystal.xsf',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("CRYSTAL")')
write(50,*)
write(50,'("PRIMVEC")')
write(50,'(3G18.10)') avec(:,1)*au2ang
write(50,'(3G18.10)') avec(:,2)*au2ang
write(50,'(3G18.10)') avec(:,3)*au2ang
write(50,*)
write(50,'("PRIMCOORD")')
write(50,'(2I8)') natmtot,1
do is=1,nspecies
  do ia=1,natoms(is)
    call r3mv(avec,atposl(1,ia,is),v1)
    write(50,'(A,3G18.10)') spsymb(is),v1*au2ang
  end do
end do
close(50)
write(*,*)
write(*,'("Info(spacegroup):")')
write(*,'(" EXCITING lattice vectors and atomic positions written to &
 &GEOMETRY.OUT")')
write(*,'(" XCrysDen file written to crystal.xsf")')
write(*,*)
stop
end program

