subroutine reorder(zrhoir,ngrid, real_z, imag_z)
implicit none

integer, intent(in) ::ngrid(3)
complex(8), intent(in) :: zrhoir(ngrid(1),ngrid(2),ngrid(3))
real(8), intent(out) :: real_z(ngrid(1),ngrid(2),ngrid(3)), imag_z(ngrid(1),ngrid(2),ngrid(3))
integer :: lx, ly, lz, lbx, lby, lbz,i, x,y,z, lbyy, lbxx, lbzz,mx,my,mz

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





if (.false.) then!pāris skaitļi
real_z(1:lx         , 1:ly        ,1:lz)         =dble(zrhoir(lbx:ngrid(1),lby:ngrid(2),lbz:ngrid(3)))
real_z(lbx:ngrid(1), 1:ly        ,1:lz)         =dble(zrhoir(1:lx         ,lby:ngrid(2),lbz:ngrid(3)))
real_z(1:lx         ,lby:ngrid(2),1:lz)         =dble(zrhoir(lbx:ngrid(1),1:ly         ,lbz:ngrid(3)))
real_z(1:lx         ,1:ly         ,lbz:ngrid(3))=dble(zrhoir(lbx:ngrid(1),lby:ngrid(2),1:lz))
real_z(lbx:ngrid(1),lby:ngrid(2),1:lz)         =dble(zrhoir(1:lx         ,1:ly         ,lbz:ngrid(3)))
real_z(lbx:ngrid(1),1:ly         ,lbz:ngrid(3))=dble(zrhoir(1:lx         ,lby:ngrid(2),1:lz))
real_z(1:lx         ,lby:ngrid(2),lbz:ngrid(3))=dble(zrhoir(lbx:ngrid(1),1:ly         ,1:lz))
real_z(lbx:ngrid(1),lby:ngrid(2),lbz:ngrid(3))=dble(zrhoir(1:lx         ,1:ly         ,1:lz))





imag_z(1:lx         , 1:ly        ,1:lz)         =dimag(zrhoir(lbx:ngrid(1),lby:ngrid(2),lbz:ngrid(3)))
imag_z(lbx:ngrid(1), 1:ly        ,1:lz)         =dimag(zrhoir(1:lx         ,lby:ngrid(2),lbz:ngrid(3)))
imag_z(1:lx         ,lby:ngrid(2),1:lz)         =dimag(zrhoir(lbx:ngrid(1),1:ly         ,lbz:ngrid(3)))
imag_z(1:lx         ,1:ly         ,lbz:ngrid(3))=dimag(zrhoir(lbx:ngrid(1),lby:ngrid(2),1:lz))
imag_z(lbx:ngrid(1),lby:ngrid(2),1:lz)         =dimag(zrhoir(1:lx         ,1:ly         ,lbz:ngrid(3)))
imag_z(lbx:ngrid(1),1:ly         ,lbz:ngrid(3))=dimag(zrhoir(1:lx         ,lby:ngrid(2),1:lz))
imag_z(1:lx         ,lby:ngrid(2),lbz:ngrid(3))=dimag(zrhoir(lbx:ngrid(1),1:ly         ,1:lz))
imag_z(lbx:ngrid(1),lby:ngrid(2),lbz:ngrid(3))=dimag(zrhoir(1:lx         ,1:ly         ,1:lz))
!
end if










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







if (.false.) then!nepāra
!old 
lx = ngrid(1)/2
ly = ngrid(2)/2
lz = ngrid(3)/2

lbx = ngrid(1)/2+1
lby = ngrid(2)/2+1
lbz = ngrid(3)/2+1
real_z(1:lx         , 1:ly        ,1:lz)         =dble(zrhoir(lbx+1:ngrid(1),lby+1:ngrid(2),lbz+1:ngrid(3)))
real_z(lx+1:ngrid(1), 1:ly        ,1:lz)         =dble(zrhoir(1:lbx         ,lby+1:ngrid(2),lbz+1:ngrid(3)))
real_z(1:lx         ,ly+1:ngrid(2),1:lz)         =dble(zrhoir(lbx+1:ngrid(1),1:lby         ,lbz+1:ngrid(3)))
real_z(1:lx         ,1:ly         ,lz+1:ngrid(3))=dble(zrhoir(lbx+1:ngrid(1),lby+1:ngrid(2),1:lbz))
real_z(lx+1:ngrid(1),ly+1:ngrid(2),1:lz)         =dble(zrhoir(1:lbx         ,1:lby         ,lbz+1:ngrid(3)))
real_z(lx+1:ngrid(1),1:ly         ,lz+1:ngrid(3))=dble(zrhoir(1:lbx         ,lby+1:ngrid(2),1:lbz))
real_z(1:lx         ,ly+1:ngrid(2),lz+1:ngrid(3))=dble(zrhoir(lbx+1:ngrid(1),1:lby         ,1:lbz))
real_z(lx+1:ngrid(1),ly+1:ngrid(2),lz+1:ngrid(3))=dble(zrhoir(1:lbx         ,1:lby         ,1:lbz))


imag_z(1:lx         , 1:ly        ,1:lz)         =dimag(zrhoir(lbx+1:ngrid(1),lby+1:ngrid(2),lbz+1:ngrid(3)))
imag_z(lx+1:ngrid(1), 1:ly        ,1:lz)         =dimag(zrhoir(1:lbx         ,lby+1:ngrid(2),lbz+1:ngrid(3)))
imag_z(1:lx         ,ly+1:ngrid(2),1:lz)         =dimag(zrhoir(lbx+1:ngrid(1),1:lby         ,lbz+1:ngrid(3)))
imag_z(1:lx         ,1:ly         ,lz+1:ngrid(3))=dimag(zrhoir(lbx+1:ngrid(1),lby+1:ngrid(2),1:lbz))
imag_z(lx+1:ngrid(1),ly+1:ngrid(2),1:lz)         =dimag(zrhoir(1:lbx         ,1:lby         ,lbz+1:ngrid(3)))
imag_z(lx+1:ngrid(1),1:ly         ,lz+1:ngrid(3))=dimag(zrhoir(1:lbx         ,lby+1:ngrid(2),1:lbz))
imag_z(1:lx         ,ly+1:ngrid(2),lz+1:ngrid(3))=dimag(zrhoir(lbx+1:ngrid(1),1:lby         ,1:lbz))
imag_z(lx+1:ngrid(1),ly+1:ngrid(2),lz+1:ngrid(3))=dimag(zrhoir(1:lbx         ,1:lby         ,1:lbz))

end if


end subroutine
