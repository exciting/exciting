subroutine kkinvert(eps2)

! purpose here is to kk invert the Im Chi(2w,w,w) to get real part
! this can be compared with the real part in ChiTotre file. They should be 
! same except in the lower energy limit. WHY????????????????
! have to find out yet!!!+

    use optica
    implicit none
    integer i,j
    complex(8) ii
    real(8) e
    real(8) eps2(emesh),eps1(emesh)
    real(8), allocatable :: ene(:)
    allocate (ene(emesh))
!..................... energy mesh construction..................................
    do i=1,emesh
     e=(i-1)*dw
     ene(i)=e*0.5d0
    enddo

    ene(1)=ene(2)/2.0
    ii=(0.d0,1.d0)

!   call broad(gamma,ene,eps2,eps2_br)
    call kk(ene,eps2,eps1)

    do j=1,emesh
      ene(j)=(j-1)*dw*13.6
      write(9,*) ene(j),eps1(j)
      write(10,*) ene(j),abs(eps1(j)+ii*eps2(j))
    end do
        
return
end subroutine

!........................... the kk inversion routine...............................
! uses the formulae:
! Real= 2/pi P { int_0^inf (w'Im(w')/((w')^2-w^2)) dw' }
!

subroutine kk(ene,eps2,eps1)
    use optica
    implicit none           

    integer i,j
    real(8) delta,sum,pi
    real(8) eps1(emesh),eps2_br(emesh),ene(emesh),f(emesh),eps2(emesh)
    delta=ene(emesh) - ene(emesh-1)

    do j = 2,emesh-1
      do i = 1,emesh
! ------------- calculate integrand --------------------------
        if(i.ne.j) then
          eps2_br(i)=eps2(i)/5.8300348177d-1
          f(i)= (eps2_br(i))*(ene(i))/           & 
                (ene(i)**2 - ene(j)**2)
        endif
      end do
! -------------------- take mean -----------------------------
!     f(j) = (f(j-1) + f(j+1))/2
      sum=0.0
! --------------------- sum calculated ---------------------
      do i = 1, emesh-1
        sum = sum+((f(i)+f(i+1))/2) * delta 
      end do
      eps1(j) = (sum*2.d0*5.8300348177d-1)/pi 
    end do
! -------------- end points --------------------------------
    eps1(1)=eps1(2)
    eps1(emesh)=eps1(emesh-1)

return
end subroutine

!............ broadening if needed..............................
subroutine broad (gamma,ept,ene,eps2,eps2_br)

    integer i,j,ept,p
    real*8 gamma,ene(ept),eps2(ept),eps2_br(ept)
    real*8 f(ept),sum
    real*8 pi2, gamma2

111 Format(f13.6,6e19.8)

    pi2=8*datan(1.0d0)
    gamma2=gamma**2/4.0
      do i=2,ept-1
        do j=1,ept
          if(j.ne.i) then
            f(j)=eps2(j)/(gamma2+(ene(j)-ene(i))**2)
          endif
        enddo
        f(i)=(f(i-1)+f(i+1))/2.0
        sum=0.0
        do p=1,ept-1
          sum=sum+((f(p)+f(p+1))/2.0)*(ene(p+1)-ene(p))
        enddo
        eps2_br(i)=gamma*sum/pi2
      enddo

      eps2_br(1)=eps2_br(2)    
      eps2_br(ept)=eps2_br(ept-1)    

return
end subroutine
	






