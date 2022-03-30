subroutine reorder(bc_type, zrhoir,ngrid, real_z, imag_z)
implicit none

integer, intent(in) ::ngrid(3), bc_type
complex(8), intent(in) :: zrhoir(ngrid(1),ngrid(2),ngrid(3))
real(8), intent(out) :: real_z(ngrid(1),ngrid(2),ngrid(3)), imag_z(ngrid(1),ngrid(2),ngrid(3))
integer :: lx, ly, lz, lbx, lby, lbz,i, x,y,z, lbyy, lbxx, lbzz,mx,my,mz, j

x = (mod(ngrid(1),2))
y = (mod(ngrid(2),2))
z = (mod(ngrid(3),2))

lx = int(ngrid(1)/2)
ly = int(ngrid(2)/2)
lz = int(ngrid(3)/2)

lbx = int(ngrid(1)/2)+1
lby = int(ngrid(2)/2)+1
lbz = int(ngrid(3)/2)+1


lbxx=int(ngrid(1)/2)+1
lbyy=int(ngrid(1)/2)+1
lbzz=int(ngrid(1)/2)+1

mx = int(ngrid(1)/2)
my = int(ngrid(2)/2)
mz = int(ngrid(3)/2)
if (x.ne.0) then

lbxx=int(ngrid(1)/2)+2


mx = int(ngrid(1)/2)+1
end if
if (y.ne.0) then

lbyy = int(ngrid(2)/2)+2
my = int(ngrid(2)/2)+1
end if

if (z.ne.0) then

lbzz = int(ngrid(3)/2)+2
mz = int(ngrid(3)/2)+1
end if



!lx = ngrid(1)/2
!ly = ngrid(2)/2
!lz = ngrid(3)/2

!lbx = ngrid(1)/2+1
!lby = ngrid(2)/2+1
!lbz = ngrid(3)/2+1
!write(*,*)lx,ly,lz
!write(*,*)lbx











if (bc_type.eq.0) then


real_z(1:lx         , 1:ly        ,1:lz)         =dble(zrhoir(lbxx:ngrid(1),lbyy:ngrid(2),lbzz:ngrid(3)))
real_z(lbx:ngrid(1), 1:ly        ,1:lz)         =dble(zrhoir(1:mx        ,lbyy:ngrid(2),lbzz:ngrid(3)))
real_z(1:lx         ,lby:ngrid(2),1:lz)         =dble(zrhoir(lbxx:ngrid(1),1:my         ,lbzz:ngrid(3)))
real_z(1:lx         ,1:ly         ,lbz:ngrid(3))=dble(zrhoir(lbxx:ngrid(1),lbyy:ngrid(2),1:mz))
real_z(lbx:ngrid(1),lby:ngrid(2),1:lz)         =dble(zrhoir(1:mx         ,1:my         ,lbzz:ngrid(3)))
real_z(lbx:ngrid(1),1:ly         ,lbz:ngrid(3))=dble(zrhoir(1:mx         ,lbyy:ngrid(2),1:mz))
real_z(1:lx         ,lby:ngrid(2),lbz:ngrid(3))=dble(zrhoir(lbxx:ngrid(1),1:my         ,1:mz))
real_z(lbx:ngrid(1),lby:ngrid(2),lbz:ngrid(3))=dble(zrhoir(1:mx         ,1:my         ,1:mz))


imag_z(1:lx         , 1:ly        ,1:lz)         =dimag(zrhoir(lbxx:ngrid(1),lbyy:ngrid(2),lbzz:ngrid(3)))
imag_z(lbx:ngrid(1), 1:ly        ,1:lz)         =dimag(zrhoir(1:mx         ,lbyy:ngrid(2),lbzz:ngrid(3)))
imag_z(1:lx         ,lby:ngrid(2),1:lz)         =dimag(zrhoir(lbxx:ngrid(1),1:my         ,lbzz:ngrid(3)))
imag_z(1:lx         ,1:ly         ,lbz:ngrid(3))=dimag(zrhoir(lbxx:ngrid(1),lbyy:ngrid(2),1:mz))
imag_z(lbx:ngrid(1),lby:ngrid(2),1:lz)         =dimag(zrhoir(1:mx         ,1:my         ,lbzz:ngrid(3)))
imag_z(lbx:ngrid(1),1:ly         ,lbz:ngrid(3))=dimag(zrhoir(1:mx         ,lbyy:ngrid(2),1:mz))
imag_z(1:lx         ,lby:ngrid(2),lbz:ngrid(3))=dimag(zrhoir(lbxx:ngrid(1),1:my         ,1:mz))
imag_z(lbx:ngrid(1),lby:ngrid(2),lbz:ngrid(3))=dimag(zrhoir(1:mx         ,1:my         ,1:mz))







end if


if (bc_type.eq.2) then!surface boundary conditions
!free in y, periodic in x and z

real_z(1:ngrid(1)         , 1:ly        ,1:ngrid(3))         =dble(zrhoir(1:ngrid(1),lbyy:ngrid(2),1:ngrid(3)))!y
imag_z(1:ngrid(1)         , 1:ly        ,1:ngrid(3))         =dimag(zrhoir(1:ngrid(1),lbyy:ngrid(2),1:ngrid(3)))

real_z(1:ngrid(1)         , lby:ngrid(2)        ,1:ngrid(3))         =dble(zrhoir(1:ngrid(1),1:my,1:ngrid(3)))!y
imag_z(1:ngrid(1)         , lby:ngrid(2)        ,1:ngrid(3))         =dimag(zrhoir(1:ngrid(1),1:my,1:ngrid(3)))


write(*,*) "surface"
end if

if (bc_type.eq.1) then!wire boundary conditions
!free in x,y, periodic in z
!----move y axis
real_z(1:ngrid(1)         , 1:ly        ,1:ngrid(3))         =dble(zrhoir(1:ngrid(1),lbyy:ngrid(2),1:ngrid(3)))!y
imag_z(1:ngrid(1)         , 1:ly        ,1:ngrid(3))         =dimag(zrhoir(1:ngrid(1),lbyy:ngrid(2),1:ngrid(3)))

real_z(1:ngrid(1)         , lby:ngrid(2)        ,1:ngrid(3))         =dble(zrhoir(1:ngrid(1),1:my,1:ngrid(3)))!y
imag_z(1:ngrid(1)         , lby:ngrid(2)        ,1:ngrid(3))         =dimag(zrhoir(1:ngrid(1),1:my,1:ngrid(3)))

!---- move x axis
real_z(1:lx         , 1:ngrid(2)        ,1:ngrid(3))         =dble(zrhoir(lbxx:ngrid(1),1:ngrid(2),1:ngrid(3)))!x
imag_z(1:lx         , 1:ngrid(2)        ,1:ngrid(3))         =dimag(zrhoir(lbxx:ngrid(1),1:ngrid(2),1:ngrid(3)))

real_z(lbx:ngrid(1)         , 1:ngrid(2)       ,1:ngrid(3))         =dble(zrhoir(1:mx,1:ngrid(2),1:ngrid(3)))!x
imag_z(lbx:ngrid(1)         , 1:ngrid(2)       ,1:ngrid(3))         =dimag(zrhoir(1:mx,1:ngrid(2),1:ngrid(3)))

!---- move z axis
!real_z(1:ngrid(1)         , 1:ngrid(2)        ,1:lz)         =dble(zrhoir(1:ngrid(1),1:ngrid(2),1:lbzz))!x
!imag_z(1:ngrid(1)         , 1:ngrid(2)        ,1:lz)         =dimag(zrhoir(1:ngrid(1),1:ngrid(2),1:lbzz))

!real_z(1:ngrid(1)         , 1:ngrid(2)       ,lbz:ngrid(3))         =dble(zrhoir(1:ngrid(1),1:ngrid(2),1:mz))!x
!imag_z(1:ngrid(1)         , 1:ngrid(2)       ,lbz:ngrid(3))         =dimag(zrhoir(1:ngrid(1),1:ngrid(2),1:mz))


write(*,*) "wire"

end if




!open (11, file = 'charge.dat', status = 'replace')

! do i = 1, ngrid(1)
!	write(11,*) real_z(i,1,1)
 !end do

 !do i = 1, ngrid(3)
  ! do j = 1, ngrid(2)
   !   write(11,*) real_z(i,j,1)
   !end do
 !end do

 !close(11)















end subroutine
