! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted 2014 for exciting, Stefan Kontur
!
!
subroutine TRANSMAT(ninter,h)
!
      use raman_ew
      use raman_inter, only : xa
      use raman_fij, only : fi
      use modinput
      implicit none
      integer :: ind,ninter,i,j,iint,ii,jj
      double precision :: u(4),v(4),h,ediff
      ind = 0
      write(66,'(/,51("*"),"   TRANSMAT   ",51("*")/)')
      write(66,'("  Transition energy [ cm-1 ]", &
     & "   Transition matrix elements: ",/ 30x, "< q >", 7x, "< q^2 >", &
     & 7x, "< q^3 >", 7x, "< q^4 >",7x, "< q^5 >", 7x, "< q^6 >",/)')
      do i = 1,input%properties%raman%nstate           ! eigen function i
         do j = i-1,1,-1                               ! eigen function j = 1...i-1
            ind = ind + 1
            transme1(ind) = 0.d0
            transme2(ind) = 0.d0
            transme3(ind) = 0.d0
            transme3(ind) = 0.d0
            transme4(ind) = 0.d0
            transme5(ind) = 0.d0
            transme6(ind) = 0.d0
            do iint = 1,ninter                         ! loop over intervals
               u(1) = z1(iint,j)                       ! u refers to j eigenfunction, v to i eigenfunction          
               u(2) = z2(iint,j)
               u(3) = z1(iint+1,j)
               u(4) = z2(iint+1,j)
               v(1) = z1(iint,i)
               v(2) = z2(iint,i)
               v(3) = z1(iint+1,i)
               v(4) = z2(iint+1,i)
               do ii = 1,4                             ! make sum u
                  do jj = 1,4                          ! make sum v
                     ! integral < u | x | v >
                     transme1(ind) = transme1(ind) + h*u(ii)*v(jj)*  &
     &                                       (xa(iint)*FI(ii,jj,0) + &
     &                                               h*FI(ii,jj,1))
                     ! integral < u | x^2 | v >
                     transme2(ind) = transme2(ind) + h*u(ii)*v(jj)*  &
     &                              (xa(iint)*xa(iint)*FI(ii,jj,0) + &
     &                                 2.d0*h*xa(iint)*FI(ii,jj,1) + &
     &                                             h*h*FI(ii,jj,2))
                     ! integral < u | x^3 | v >
                     transme3(ind) = transme3(ind) + h*u(ii)*v(jj)*  &
     &                     (xa(iint)*xa(iint)*xa(iint)*FI(ii,jj,0) + &
     &                        3.d0*xa(iint)*xa(iint)*h*FI(ii,jj,1) + &
     &                               3.d0*xa(iint)*h*h*FI(ii,jj,2) + &
     &                                           h*h*h*FI(ii,jj,3))
                     ! integral < u | x^4 | v >
                     transme4(ind) = transme4(ind) + h*u(ii)*v(jj)*  &
     &            (xa(iint)*xa(iint)*xa(iint)*xa(iint)*FI(ii,jj,0) + &
     &               4.d0*xa(iint)*xa(iint)*xa(iint)*h*FI(ii,jj,1) + &
     &                      6.d0*xa(iint)*xa(iint)*h*h*FI(ii,jj,2) + &
     &                             4.d0*xa(iint)*h*h*h*FI(ii,jj,3) + &
     &                                         h*h*h*h*FI(ii,jj,4))
                     ! integral < u | x^5 | v >
                     transme5(ind) = transme5(ind) + h*u(ii)*v(jj)*  &
     &                                  ((xa(iint)**5)*FI(ii,jj,0) + &
     &                            5.d0*(xa(iint)**4)*h*FI(ii,jj,1) + &
     &                      10.d0*(xa(iint)**3)*(h**2)*FI(ii,jj,2) + &
     &                      10.d0*(xa(iint)**2)*(h**3)*FI(ii,jj,3) + &
     &                            5.d0*xa(iint)*(h**4)*FI(ii,jj,4) + &
     &                                          (h**5)*FI(ii,jj,5))
                     ! integral < u | x^6 | v >
                     transme6(ind) = transme6(ind) + h*u(ii)*v(jj)*  &
     &                                  ((xa(iint)**6)*FI(ii,jj,0) + &
     &                            6.d0*(xa(iint)**5)*h*FI(ii,jj,1) + &
     &                      15.d0*(xa(iint)**4)*(h**2)*FI(ii,jj,2) + &
     &                      20.d0*(xa(iint)**3)*(h**3)*FI(ii,jj,3) + &
     &                      15.d0*(xa(iint)**2)*(h**4)*FI(ii,jj,4) + &
     &                            6.d0*xa(iint)*(h**5)*FI(ii,jj,5) + &
     &                                          (h**6)*FI(ii,jj,6))
                  enddo
               enddo
            enddo
!           transme1(ind) = transme1(ind)/Sqrt(fgew)   ! scale by fgew
!           transme2(ind) = transme2(ind)/fgew
!           transme3(ind) = transme3(ind)/fgew/Sqrt(fgew)
!           transme4(ind) = transme4(ind)/fgew/fgew
!           transme5(ind) = transme5(ind)/fgew/fgew/Sqrt(fgew)
!           transme6(ind) = transme6(ind)/fgew/fgew/fgew
            ediff = eigen(i) - eigen(j)
            ediff = ediff*fhawn
            write(66,'(i3," -> ",i2,f14.2,6f14.8)')         &
     &                             j,i,ediff,transme1(ind), &
     &                             transme2(ind),           &
     &                             transme3(ind),           &
     &                             transme4(ind),           &
     &                             transme5(ind),           &
     &                             transme6(ind)
         enddo
      enddo
      write(66,*)
      return
!
      end
!
!
