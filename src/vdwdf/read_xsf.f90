!-----------------------------------------------------------------------
subroutine read_xsf(filename)

  use param
  
  implicit none
  character*80     :: filename
  integer          :: i, j, k
  real*8           :: crossab(3)
 
  character        :: line*120
  integer          :: len, iostat
  integer          :: i_trimleft_white_space
  
  real*8           :: sum
  
  open(11,file=filename,status='old')
  
  do 
      read(11,'(a120)') line
      len = i_trimleft_white_space(line)

      if ( (line(1:23)=='BEGIN_BLOCK_DATAGRID_3D').or. &
           (line(1:22)=='BEGIN_BLOCK_DATAGRID3D') )&
      then
      	  read(11,*);  read(11,*)
          read(11,*) nx,ny,nz
      	  allocate(density(nx,ny,nz))
      	  allocate(graddensity(nx,ny,nz))
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
! Uncomment this line when reading the VASP's xsf file
! density = density * (bohr**3.d0) 

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

! density gradients
  call evalGradient()

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

end subroutine read_xsf

!-----------------------------------------------------------------------
subroutine evalGradient()
use param
implicit none
integer :: i, j, k
real*8  :: a, b, c, g(3,3), T(3,3)
real*8  :: hx, hy, hz
real*8  :: grad_x, grad_y, grad_z
real*8  :: gr_x, gr_y, gr_z

! remove boundary conditions
  nx = nx - 1
  ny = ny - 1
  nz = nz - 1

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
  hx=a/dble(nx)
  hy=b/dble(ny)
  hz=c/dble(nz)
!  print*, 'hx,hy,hz', hx,hy,hz
  
  do i=1,nx
  do j=1,ny
  do k=1,nz
  
     if (i.gt.1) then
       grad_x=(density(i+1,j  ,k  )-density(i-1,j  ,k  ))/(2*hx)     
     else 
       grad_x=(density(2  ,j  ,k  )-density(nx ,j  ,k  ))/(2*hx)
     end if
   
     if (j.gt.1) then
       grad_y=(density(i  ,j+1,k  )-density(i  ,j-1,k  ))/(2*hy)       
     else 
       grad_y=(density(i  ,2  ,k  )-density(i  ,ny ,k  ))/(2*hy)
     end if
   
     if (k.gt.1) then
       grad_z=(density(i  ,j  ,k+1)-density(i  ,j  ,k-1))/(2*hz)       
     else 
       grad_z=(density(i  ,j  ,2  )-density(i  ,j  ,nz ))/(2*hz)
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

  do i=1,nx+1
  do j=1,ny+1
    graddensity(i,j,nz+1)=graddensity(i,j,1) 
  end do  
  end do
  
  do j=1,ny+1
  do k=1,nz+1
    graddensity(nx+1,j,k)=graddensity(1,j,k) 
  end do  
  end do
  
  do i=1,nx+1
  do k=1,nz+1
    graddensity(i,ny+1,k)=graddensity(i,1,k) 
  end do  
  end do

  !write(99,'(6f13.5)') graddensity
  !stop

return
end

!-------------------------------------------------
integer function i_trimleft_white_space(word)
!-------------------------------------------------
! trim left white spaces out of word
!-------------------------------------------------
  character word*(*), auxword*80

  ilen=len(word)
  auxword=word
  do i=1,ilen
     if ( word(i:i) .eq. ' ' ) then
        auxword=word(i+1:ilen)
     else
        exit
     endif
  enddo
  i_trimleft_white_space=len(word)
  word=auxword(1:i_trimleft_white_space)
return
end
