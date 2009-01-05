
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine genveelu(l,u,j,lmax,vee)
implicit none
! arguments
integer, intent(in) :: l
real(8), intent(in) :: u
real(8), intent(in) :: j
integer, intent(in) :: lmax
real(8), intent(out) :: vee(-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax)
! local variables
integer m1,m2,m3,m4,k,q
real(8), parameter :: fourpi=12.566370614359172954d0
real(8) r1,r2,f(0:6)
real(8) sum1,sum2,t1
! external functions
real(8) gaunt
external gaunt
! Slater integrals F(k) for d and f electrons in Ry, to be converted in Htr
f(:)=0.d0
f(0)=u
select case(l)
case(0)
! s electrons only f(0)=u
case(1)
! p electrons
  f(2)=5.d0*j
case(2)
! d electrons: ratio r1 = F(4)/F(2), see PRB 52, R5467 (1995)
  r1=0.625d0
  f(2)=(14.d0*j)/(1.d0+r1)
  f(4)=f(2)*r1
case(3)
! f electrons: r2 = F(6)/F(2), r1 = F(4)/F(2), see PRB 50, 16861 (1994)
  r1=451.d0/675.d0
  r2=1001.d0/2025.d0
  f(2)=6435.d0*j/(286.d0+195.d0*r1+250.d0*r2)
  f(4)=f(2)*r1
  f(6)=f(2)*r2
case default
  write(*,*)
  write(*,'("Error(genveelu): invalid l : ",I8)') l
  write(*,*)
  stop
end select
do m1=-l,l
  do m2=-l,l
    do m3=-l,l
      do m4=-l,l
        sum1=0.d0
        do k=0,2*l,2
          sum2=0.d0
          do q=-k,k
            t1=gaunt(l,k,l,m1,q,m2)*gaunt(l,k,l,m3,-q,m4)
            if (mod(q,2).eq.0) then
              sum2=sum2+t1
            else
              sum2=sum2-t1
            end if
          end do
          sum1=sum1+f(k)*sum2/dble(2*k+1)
        end do
        vee(m1,m3,m2,m4)=fourpi*sum1
      end do
    end do
  end do
end do
return
end subroutine

