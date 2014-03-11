!
subroutine TEMPPOT(temp)
!
!
! comnputes potential coefficients for solids and temperature temp
!
      use modinput
      use raman_ew
      use raman_inter, only : xa,h
      use raman_fij, only : fi
      use raman_input, only : ncell
      use raman_coeff, only : A0,A1,A2,A3,A4,A5,A6
      implicit none
      integer :: ind,i,iint,ii,jj
      real(8) :: u(4),ua(6),ex(input%properties%raman%nstate)
      real(8) :: arg,beta,zsum,temp,sni
      ind = 0
      write(66,95)
      write(66,96)
      do i = 1,input%properties%raman%nstate                                 ! eigen function i
         ind = ind + 1
         transme1(ind) = 0.d0
         transme2(ind) = 0.d0
         transme3(ind) = 0.d0
         transme3(ind) = 0.d0
         transme4(ind) = 0.d0
         transme5(ind) = 0.d0
         transme6(ind) = 0.d0
         do iint = 1,input%properties%raman%ninter                         ! loop over intervals
            u(1) = z1(iint,i)                       ! u refers to j eigenfunction, v to i eigenfunction          
            u(2) = z2(iint,i)
            u(3) = z1(iint+1,i)
            u(4) = z2(iint+1,i)
            do ii = 1,4                             ! make sum u left
               do jj = 1,4                          ! make sum u right
                  ! integral < u | x | u >
                  transme1(ind) = transme1(ind) + h*u(ii)*u(jj)*  &
     &                                    (xa(iint)*FI(ii,jj,0) + &
     &                                            h*FI(ii,jj,1))
                  ! integral < u | x^2 | u >
                  transme2(ind) = transme2(ind) + h*u(ii)*u(jj)*  &
     &                           (xa(iint)*xa(iint)*FI(ii,jj,0) + &
     &                              2.d0*h*xa(iint)*FI(ii,jj,1) + &
     &                                          h*h*FI(ii,jj,2))
                  ! integral < u | x^3 | u >
                  transme3(ind) = transme3(ind) + h*u(ii)*u(jj)*  &
     &                  (xa(iint)*xa(iint)*xa(iint)*FI(ii,jj,0) + &
     &                     3.d0*xa(iint)*xa(iint)*h*FI(ii,jj,1) + &
     &                            3.d0*xa(iint)*h*h*FI(ii,jj,2) + &
     &                                        h*h*h*FI(ii,jj,3))
                  ! integral < u | x^4 | u >
                  transme4(ind) = transme4(ind) + h*u(ii)*u(jj)*  &
     &         (xa(iint)*xa(iint)*xa(iint)*xa(iint)*FI(ii,jj,0) + &
     &            4.d0*xa(iint)*xa(iint)*xa(iint)*h*FI(ii,jj,1) + &
     &                   6.d0*xa(iint)*xa(iint)*h*h*FI(ii,jj,2) + &
     &                          4.d0*xa(iint)*h*h*h*FI(ii,jj,3) + &
     &                                      h*h*h*h*FI(ii,jj,4))
                  ! integral < u | x^5 | u >
                  transme5(ind) = transme5(ind) + h*u(ii)*u(jj)*  &
     &                               ((xa(iint)**5)*FI(ii,jj,0) + &
     &                         5.d0*(xa(iint)**4)*h*FI(ii,jj,1) + &
     &                   10.d0*(xa(iint)**3)*(h**2)*FI(ii,jj,2) + &
     &                   10.d0*(xa(iint)**2)*(h**3)*FI(ii,jj,3) + &
     &                         5.d0*xa(iint)*(h**4)*FI(ii,jj,4) + &
     &                                       (h**5)*FI(ii,jj,5))
                  ! integral < u | x^6 | u >
                  transme6(ind) = transme6(ind) + h*u(ii)*u(jj)*  &
     &                               ((xa(iint)**6)*FI(ii,jj,0) + &
     &                         6.d0*(xa(iint)**5)*h*FI(ii,jj,1) + &
     &                   15.d0*(xa(iint)**4)*(h**2)*FI(ii,jj,2) + &
     &                   20.d0*(xa(iint)**3)*(h**3)*FI(ii,jj,3) + &
     &                   15.d0*(xa(iint)**2)*(h**4)*FI(ii,jj,4) + &
     &                         6.d0*xa(iint)*(h**5)*FI(ii,jj,5) + &
     &                                       (h**6)*FI(ii,jj,6))
               enddo
            enddo
         enddo
         transme1(ind) = transme1(ind)/Sqrt(fgew)   ! scale by fgew
         transme2(ind) = transme2(ind)/fgew
         transme3(ind) = transme3(ind)/fgew/Sqrt(fgew)
         transme4(ind) = transme4(ind)/fgew/fgew
         transme5(ind) = transme5(ind)/fgew/fgew/Sqrt(fgew)
         transme6(ind) = transme6(ind)/fgew/fgew/fgew
         write(66,94) i,transme1(ind), &
     &                  transme2(ind), &
     &                  transme3(ind), &
     &                  transme4(ind), &
     &                  transme5(ind), &
     &                  transme6(ind)
      enddo
      write(66,*)
! compute factors
      beta = 1.d0/fkwn/temp
      do i = 1,input%properties%raman%nstate
         arg = -beta*eigen(i)*frywn
         ex(i) = dexp(arg)
         zsum = zsum + ex(i)                        ! zsum = Sum( exp(-E/kT) )
      enddo
      ua = 0.d0
      sni = 1.d0/dsqrt(dble(ncell))
      do i = 1,input%properties%raman%nstate
         ua(1) = ua(1) + ex(i)/zsum*transme1(i)*(1.d0 - sni)
         ua(2) = ua(2) + ex(i)/zsum*transme2(i)*(1.d0 - sni**2)
         ua(3) = ua(3) + ex(i)/zsum*transme3(i)*(1.d0 - sni**3)
         ua(4) = ua(4) + ex(i)/zsum*transme4(i)*(1.d0 - sni**4)
         ua(5) = ua(5) + ex(i)/zsum*transme5(i)*(1.d0 - sni**5)
         ua(6) = ua(6) + ex(i)/zsum*transme6(i)*(1.d0 - sni**6)
      enddo
      write(66,98) (ua(i),i=1,6)
! compute coefficients incl. contributions of of other phonons
      a0 = a0 +      a1*ua(1) +       a2*ua(2) +       a3*ua(3) +       a4*ua(4) +      a5*ua(5) + a6*ua(6)
      a1 = a1 + 2.d0*a2*ua(1) +  3.d0*a3*ua(2) +  4.d0*a4*ua(3) +  5.d0*a5*ua(4) + 6.d0*a6*ua(5)
      a2 = a2 + 3.d0*a3*ua(1) +  6.d0*a4*ua(2) +  9.d0*a5*ua(3) + 12.d0*a6*ua(4)
      a3 = a3 + 4.d0*a4*ua(1) +  8.d0*a5*ua(2) + 12.d0*a6*ua(3)
      a4 = a4 + 5.d0*a5*ua(1) + 10.d0*a6*ua(2)
      a5 = a5 + 6.d0*a6*ua(1)
      a6 = a6
      write(66,100) a0,a1,a2,a3,a4,a5,a6
      return
 96   format('   Expectation values [ cm-1 ]   ',// 10x, '< q >', 8x, '< q^2 >', &
     & 7x, '< q^3 >', 7x, '< q^4 >',7x, '< q^5 >', 7x, '< q^6 >',/)
 94   format(i3,6f14.8)
 95   format(/,51('*'),'   TEMPPOT   ',51('*')/)
 98   format(/,' Contributions u^1...u^6 :   ',6f15.6,/)
100   format(/' New coefficients are :',7f11.4)
!
end subroutine temppot
!
!
