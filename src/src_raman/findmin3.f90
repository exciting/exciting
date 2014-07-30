! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
!
subroutine FINDMIN3(zmin)
!
! compute minima of polynomial describing the potential
! one minimum is expected between xmin and xmax
!
use raman_coeff, only : a1,a2,a3,a4,a5,a6
use modinput
implicit none
double precision :: work(100),d2,da(6),zmin
double precision, allocatable :: cm(:,:),rootr(:),rooti(:),vec(:)
integer :: deg,info,i,j
!
write(66,84)
write(66,90) input%properties%raman%degree
! search roots for d/dx V(x) with degree-1
da(1) = a1; da(2) = 2.d0*a2
da(3) = 3.d0*a3; da(4) = 4.d0*a4
da(5) = 5.d0*a5; da(6) = 6.d0*a6
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
   ! call lapack routine DGEEV to compute eigenvalues
   call DGEEV('N','N',deg,cm,deg,rootr,rooti,vec,deg,vec,deg,WORK,100,INFO)
   !
   ! select roots within interal xmin..xmax only
   ! and determine minima
   j = 0
   do i = 1,deg
      if (rooti(i) .eq. 0.d0) then
         d2 = 2.d0*a2          +     6.d0*a3*rootr(i) +  12.d0*a4*rootr(i)**2 + &
            & 2.d1*a5*rootr(i)**3 + 36.d0*a6*rootr(i)**4
         ! adjust limits in case we found a local maximum
         if ((rootr(i) .gt. input%properties%raman%xmin) .and. (rootr(i) .lt. 0.d0) .and. (d2 .lt. 0.d0)) then
            input%properties%raman%xmin = rootr(i)
            cycle
         endif
         if ((rootr(i) .lt. input%properties%raman%xmax) .and. (rootr(i) .gt. 0.d0) .and. (d2 .lt. 0.d0)) then
            input%properties%raman%xmax = rootr(i)
            cycle
         endif
         if ((rootr(i) .ge. input%properties%raman%xmin) .and. &
          &  (rootr(i) .le. input%properties%raman%xmax) .and. (d2 .gt. 0.d0)) then
            zmin = rootr(i)
            j = j + 1
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
write(66, '(/," The oscillator problem for this mode is solved between u = ",f7.3," and ",f7.3," Bohr.")') &
 &   input%properties%raman%xmin, input%properties%raman%xmax
write(66, '(/," Potential minimum at: ",f12.6,/)') zmin
!
deallocate( cm,rootr,rooti,vec )
84 format (//,43('*'),'   Find minimum of potential   ',42('*'),/)
90 format (/,' Minima for polynomial of degree',i3,' are computed...',/)
return
end
!
!
