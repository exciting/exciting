! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!-----------------------------------------------------------------------
subroutine read_gradients()

  use param
  
  implicit none
  integer          :: i, j, k
  real*8           :: crossab(3)
 
  character        :: line*120
  integer          :: iostat
  integer          :: gnx, gny, gnz
  
  real*8           :: sum
  
  open(11,file=trim(xsfgradients),status='old',ERR=10)

  do                                                                    
      read(11,'(a120)') line                                            
      line=adjustl(line)                                                

      if ( (line(1:23)=='BEGIN_BLOCK_DATAGRID_3D').or. &                
           (line(1:22)=='BEGIN_BLOCK_DATAGRID3D') )&                    
      then                                                              
          read(11,*);  read(11,*)                                       
          read(11,*) gnx,gny,gnz
          ! check for consistence of the gradient grid with the density one
          if ((nx/=gnx).or.(ny/=gny).or.(nz/=gnz)) then
              print*, "ERROR in read_density.f90 : Inconsistent density and density gradients grids"
              stop
          end if
          allocate(graddensity(nx,ny,nz))                            
          read(11,*); read(11,*); read(11,*); read(11,*)
          ! to skip empty lines which sometimes precede density data    
          iostat=1                                                      
          do while (iostat.ne.0)                                        
              read(11,*,IOSTAT=iostat) graddensity
          end do                                                        
          exit                                                          
      end if
  end do                                                                
  close(11)
  
  ! Important:
  !   - in Espresso density is electron/(bohr**3)
  !   - in VASP density is electron/(angstrom**3)

  ! Convert the density gradients into program's internal units
  if (densunits=='angs') graddensity = graddensity * (bohr**4.0d0)
  return

! else if the gradient file does not exist
  10  write(*,*) 'read_density.f90 : Density gradients are calculated using a three-point formula'
  call evalGradients()

end subroutine read_gradients

!-----------------------------------------------------------------------
subroutine evalGradients()
use param
implicit none
integer :: i, j, k
real*8  :: a, b, c, g(3,3), T(3,3)
real*8  :: hx, hy, hz
real*8  :: grad_x, grad_y, grad_z
real*8  :: gr_x, gr_y, gr_z

  allocate(graddensity(nx,ny,nz))

! keep in mind the periodic boundary conditions,
! i.e., density(1,*)=density(nx,*)

! unit cell lengths
  a=dsqrt(vectors(1,1)**2d0+vectors(2,1)**2d0+vectors(3,1)**2d0)
  b=dsqrt(vectors(1,2)**2d0+vectors(2,2)**2d0+vectors(3,2)**2d0)
  c=dsqrt(vectors(1,3)**2d0+vectors(2,3)**2d0+vectors(3,3)**2d0)

! transformational matrix for the basis vectors
  g(1,:)=vectors(:,1)/a
  g(2,:)=vectors(:,2)/b
  g(3,:)=vectors(:,3)/c
  !write(*,*) (g(i,:),i=1,3)
  
! gradient transformation matrix
  call r3minv(g,T)
  !write(*,*) (T(i,:),i=1,3)  

! differentiation step
  hx=a/dble(nx-1)
  hy=b/dble(ny-1)
  hz=c/dble(nz-1)
!  print*, 'hx,hy,hz', hx,hy,hz
  
  do i=1,nx-1
  do j=1,ny-1
  do k=1,nz-1
  
     if (i.gt.1) then
       grad_x=(density(i+1,j,k)-density(i-1,j,k))/(2*hx)     
     else 
       grad_x=(density(2,j,k)-density(nx-1,j,k))/(2*hx)
     end if
   
     if (j.gt.1) then
       grad_y=(density(i,j+1,k)-density(i,j-1,k))/(2*hy)       
     else 
       grad_y=(density(i,2,k)-density(i,ny-1,k))/(2*hy)
     end if
   
     if (k.gt.1) then
       grad_z=(density(i,j,k+1)-density(i,j,k-1))/(2*hz)       
     else 
       grad_z=(density(i,j,2)-density(i,j,nz-1))/(2*hz)
     end if

     ! in case of not rectangular cells one needs to transform
     ! the gradient components
     gr_x=T(1,1)*grad_x+T(2,1)*grad_y+T(3,1)*grad_z
     gr_y=T(1,2)*grad_x+T(2,2)*grad_y+T(3,2)*grad_z
     gr_z=T(1,3)*grad_x+T(2,3)*grad_y+T(3,3)*grad_z
     
     graddensity(i,j,k)=dsqrt(gr_x*gr_x+gr_y*gr_y+gr_z*gr_z)
      
  end do
  end do
  end do

  ! the periodic boundary conditions for gradients
  do i=1,nx
  do j=1,ny
    graddensity(i,j,nz)=graddensity(i,j,1) 
  end do  
  end do
  
  do j=1,ny
  do k=1,nz
    graddensity(nx,j,k)=graddensity(1,j,k) 
  end do  
  end do
  
  do i=1,nx
  do k=1,nz
    graddensity(i,ny,k)=graddensity(i,1,k) 
  end do  
  end do

  !write(99,'(6f13.5)') graddensity
  !stop

return
end
