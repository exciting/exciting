subroutine pickritzvectors(n,segmentsize,evals,firstvec)
      use mod_APW_LO
      implicit none
      integer :: n,firstvec,segmentsize
      real(8) :: evals(n)

      integer :: i
      real(8) :: ta,tb

      call timesec(ta)      
      i=1      
      do while ((i.le.n-segmentsize).and.(evals(i).lt.mine0))
       i=i+1
      enddo
      if (i.le.n-segmentsize+1) then
        firstvec=i
      else
        write(*,*) 'Eigenvalue filtering:'
        write(*,*) 'Too many eigenvalues under the permitted bound of ',mine0
        stop
      endif 
      call timesec(tb)
!write(*,*) 'pickritzvectors', tb-ta

end subroutine pickritzvectors

