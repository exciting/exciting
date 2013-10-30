!
subroutine FINDMIN3(zmin)
!
! compute minima of polynomial describing the potential
! one minimum is expected between xmin and xmax
!
use raman_coeff, only : a1,a2,a3,a4,a5,a6
use raman_inter, only : xmin_r,xmax_r
use modinput
implicit none
double precision :: work(100),d2,da(6),zmin
double precision, allocatable :: cm(:,:),rootr(:),rooti(:),vec(:)
integer :: deg,info,i,j
!
write(66,84)
! determine degree of polynomial
!degree = 0
!if (a6 .eq. 0.d0) then
!   if (a5 .eq. 0.d0) then
!      if (a4 .eq. 0.d0) then
!         if (a3 .eq. 0.d0) then
!            if (a2 .eq. 0.d0) then
!                if (a1 .eq. 0.d0) degree = 1
!            else
!                degree = 2
!            endif
!         else
!            degree = 3
!         endif
!      else
!         degree = 4
!      endif
!   else
!      degree = 5
!   endif
!else
!   degree = 6
!endif
write(66,90) input%properties%raman%degree
! search roots for d/dx V(x) with degree-1
da(1) = a1; da(2) = 2.d0*a2
da(3) = 3.d0*a3; da(4) = 4.d0*a4
da(5) = 5.d0*a5; da(6) = 6.d0*a6
!do while (abs(da(input%properties%raman%degree)) .lt. 1.d-12) 
!   degree = degree - 1
!enddo
deg = input%properties%raman%degree - 1
! construct companion matrix
allocate( cm(deg,deg) )
allocate( rootr(deg) )
allocate( rooti(deg) )
allocate( vec(deg) )
cm = 0.d0
do i = 1,deg
   cm(i,deg) = -da(i)/da(deg+1)
   if (i .gt. 1) cm(i,i-1) = 1.d0
enddo
!
!if (deg .ge. 2) then
   ! call lapack routine DGEEV to compute eigenvalues
   call DGEEV('N','N',deg,cm,deg,rootr,rooti,vec,deg,vec,deg,WORK,100,INFO)
   !
   ! select roots within interal xmin..xmax only
   ! and determine minima
   j = 0
   do i = 1,deg
      if (rooti(i) .eq. 0.d0) then
         if ((rootr(i) .ge. xmin_r) .and. (rootr(i) .le. xmax_r)) then
            d2 = 2.d0*a2          +     6.d0*a3*rootr(i) +  12.d0*a4*rootr(i)**2 + &
               & 2.d1*a5*rootr(i)**3 + 36.d0*a6*rootr(i)**4
            if (d2 .gt. 0.d0) then
               zmin = rootr(i)
               j = j + 1
            endif
         endif
      endif
   enddo
   if (j .ne. 1) then
      write(66,*) 'Warning! ',j,' minima found in potential'
      write(*,*) 'Warning! ',j,' minima found in potential'
      write(*,*) 'Roots are [No., Re, Im]'
      write(*,'(i6,2f12.6)') (i,rootr(i),rooti(i),i=1,deg)
      stop
   endif
!else
!   ! for a quadratic function do simply the following
!   zmin = - a1 / a2
!endif
write(66,94) zmin
!
deallocate( cm,rootr,rooti,vec )
84 format (//,43('*'),'   Find minimum of potential   ',42('*'),/)
90 format (/,' Minima for polynomial of degree',i3,' are computed...',/)
94  format(/,' Minimum at: ',f12.6,/)
return
end
!
!
