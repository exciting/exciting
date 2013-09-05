!
subroutine POTENTIAL(maxp)
!
   use modinput
   use raman_inter
   implicit none
   integer :: maxp,i,iidim
   real(8) :: xende
!
   call potx6(maxp)
!
   iidim = 2*input%properties%raman%ninter
   do i = 1,input%properties%raman%ninter       ! change mesh, xa(i+1)-xa(i) contains then 4 knots
      xa(i) = xpot(4*i-3)
   enddo
   h = xa(2) - xa(1)
   xende = xa(input%properties%raman%ninter) + h
   call interpol
   return
end subroutine potential
!
