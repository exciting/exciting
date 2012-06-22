!BOP
! !ROUTINE: e1xb
!
! !INTERFACE:
        function e1xb(x) result(e1)
!
! !DESCRIPTION:
!
! compute the exponential integral function e1(x) 
!
! !INPUT PARAMETERS:
      implicit none
      
      real(8) :: x
      
!
! !OUTPUT PARAMETERS:
      real(8) :: e1
      
! !LOCAL VARIABLES:

      integer(4) :: k
      integer(4) :: m
      
      real(8)    :: r      
      real(8)    :: ga
      real(8)    :: t0
      real(8)    :: t
            
!
!EOP
!BOC
      if (x.eq.0.0) then
        e1=1.0d+300
      else if (x.le.1.0) then
        e1=1.0d0
        r=1.0d0
        do 10 k=1,25
          r=-r*k*x/(k+1.0d0)**2
          e1=e1+r
          if (dabs(r).le.dabs(e1)*1.0d-15) go to 15
10      continue
15      ga=0.5772156649015328d0
        e1=-ga-dlog(x)+x*e1
      else
        m=20+int(80.0/x)
        t0=0.0d0
        do 20 k=m,1,-1
          t0=k/(1.0d0+k/(x+t0))
20      continue
        t=1.0d0/(x+t0)
        e1=dexp(-x)*t
      endif
      return
      end function e1xb
!EOC
