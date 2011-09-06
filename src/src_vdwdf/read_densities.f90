! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!-----------------------------------------------------------------------
subroutine read_densities()

  use param
  
  implicit none
  integer          :: i, j, k
  real*8           :: crossab(3)
 
  character        :: line*120
  integer          :: iostat
  
  real*8           :: sum
  
  open(11,file=trim(xsffile),status='old')
  
  do 
      read(11,'(a120)') line
      line=adjustl(line)

      if ( (line(1:23)=='BEGIN_BLOCK_DATAGRID_3D').or. &
           (line(1:22)=='BEGIN_BLOCK_DATAGRID3D') )&
      then
      	  read(11,*);  read(11,*)
          read(11,*) nx,ny,nz
          ! note the PBC, e.i., density(1,*)=density(nx,*)
      	  allocate(density(nx,ny,nz))
      	  read(11,*) origin
          read(11,*) vectors 
      	  ! to skip empty lines which sometimes precede density data
          iostat=1
          do while (iostat.ne.0)
              read(11,*,IOSTAT=iostat) density
          end do
	  exit
      end if
  
  end do
  close(11)
  
!====================================
! The program internal units:
!   - distance in bohr
!   - density  in eletron/(bohr**3)
!------------------------------------

! convert Angstroms to bohrs  
  vectors = vectors/bohr

! Important:
!   - in Espresso density is electron/(bohr**3)
!   - in VASP density is electron/(angstrom**3)
!

! Convert the input density into program internal units
  if (densunits=='angs') density = density * (bohr**3.d0)

!====================================

! Checking input density for negative values and zeros  
  do i=1,nx
  do j=1,ny
  do k=1,nz
     if ((density(i,j,k).lt.0d0).or.(abs(density(i,j,k)).lt.eps)) then
          density(i,j,k)=1.0d-10
     end if
  end do
  end do
  end do

! calculate unit cell volume  
  crossab(1) = vectors(2,1)*vectors(3,2) - vectors(3,1)*vectors(2,2)
  crossab(2) = vectors(3,1)*vectors(1,2) - vectors(1,1)*vectors(3,2)
  crossab(3) = vectors(1,1)*vectors(2,2) - vectors(2,1)*vectors(1,2)
  volume     = crossab(1)*vectors(1,3) + crossab(2)*vectors(2,3) + crossab(3)*vectors(3,3)

  ! Check the density for correct units
  !sum=0.0
  !do i=1,nx
  !do j=1,ny
  !do k=1,nz
  !    sum=sum+density(i,j,k)
  !end do
  !end do
  !end do
  !sum=volume*sum/(nx*ny*nz)
  !write(*,*) "The total density = ", sum
  !stop

end subroutine read_densities
