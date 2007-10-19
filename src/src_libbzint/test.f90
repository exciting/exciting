program test

integer(4):: i,j,ind

      do i=1,4
       do j=1,4
        ind=mod(j+i-2,4)+1
        write(*,*)i,j,ind
       enddo
      enddo
end
