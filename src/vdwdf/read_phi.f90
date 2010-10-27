subroutine read_phi(filename)

  use param
  
  implicit none
  character*80     :: filename
  integer          :: j
  real*8           :: d
  
  open(11,file=filename,status='old')
  read(11,*) dmin,dmax,dstep
  read(11,*) deltamin,deltamax,deltastep 
  read(11,*) nd, ndelta
  
  allocate(kernel(ndelta,nd))
  
  read(11,*) kernel
  
  do j = 1, nd
      d           = dmin + (j-1)*dstep
      kernel(:,j) = kernel(:,j)/(4.0d0*pi*d*d) 
  end do
  
  close(11)

end subroutine read_phi
