!///////////////////////////////////////////////////////////////////////////////

!!$program test_summations
!!$  use summations
!!$  implicit none
!!$
!!$  integer, parameter :: n1=200,n2=40,m1=2000,m2=2000
!!$  complex(8) :: z2(n1,n2),z(n1,n2),za(m1,n1),za2(n1,m1),zb(m2,n2),zm(m1,m2)
!!$  complex(8) :: zw(m1,n2),zw2(n1,m2),zi,alpha,beta
!!$  real(8) :: rz(n1,n2),ra(m1,n1),rb(m2,n2),rm(m1,m2)
!!$  real(8) :: cpu0,cpu1,cpuzgemm,cpumatmul
!!$
!!$  zi=cmplx(0.d0,1.d0,8)
!!$  alpha=1.d0
!!$  beta=0.d0
!!$
!!$  za=0.d0
!!$  zb=0.d0
!!$  zm=0.d0
!!$  call random_number(ra)
!!$  call random_number(rb)
!!$  call random_number(rm)
!!$  za=za+ra/dble(m1)
!!$  zb=zb+rb/dble(m2)
!!$  zm=zm+rm
!!$  call random_number(ra)
!!$  call random_number(rb)
!!$  call random_number(rm)
!!$  za=za+zi*ra/dble(m1)
!!$  zb=zb+zi*rb/dble(m2)
!!$  zm=zm+zi*rm
!!$
!!$  za2=transpose(za)
!!$
!!$  call cpu_time(cpu0)
!!$  call doublesummation_simple_cz(z,za,zm,zb,alpha,beta)
!!$  call cpu_time(cpu1)
!!$  cpuzgemm=cpu1-cpu0
!!$if (n1.ge.n2) then
!!$   zw=matmul(zm,zb)
!!$   z2=matmul(conjg(transpose(za)),zw)
!!$else
!!$   zw2=matmul(conjg(transpose(za)),zm)
!!$   z2=matmul(zw2,zb)
!!$end if
!!$  call cpu_time(cpu0)
!!$  cpumatmul=cpu0-cpu1
!!$  write(111,*) z
!!$  write(222,*) z2
!!$  write(333,*) z-z2
!!$
!!$write(*,'(a,f12.3)') 'CPU-time (zgemm) :',cpuzgemm
!!$write(*,'(a,f12.3)') 'CPU-time (matmul):',cpumatmul
!!$write(*,*) 'max deviation:',maxval(abs(z-z2))
!!$
!!$end program test_summations

subroutine test_blaswrappers
  ! PASSED TEST November 10, 2008
  use blaswrappers
  implicit none
  integer, parameter :: maxcase=9
  integer, parameter :: nz1=200,nz2=400, na1=200,na2=2000, nb1=2000,nb2=400
  complex(8), parameter :: alpha=(2.5d0,4.7d0),beta=(-1.6d0,3.1d0)
  complex(8), parameter :: zi=(0.d0,1.d0)
  character(1) :: ata(maxcase),atb(maxcase),ta,tb
  integer :: ama1(maxcase),ama2(maxcase),amb1(maxcase),amb2(maxcase)
  integer :: icase,ma1,ma2,mb1,mb2
  real(8) :: cpu0,cpu1,cpuzgemm,cpuzgemmwrap,cpumatmul
  real(8), allocatable :: ra(:,:),rb(:,:)
  complex(8), allocatable :: z(:,:),z2(:,:),z3(:,:),za(:,:),zb(:,:)
  allocate(z(nz1,nz2),z2(nz1,nz2),z3(nz1,nz2))
  ata=(/'n','n','n','t','t','t','c','c','c'/)
  atb=(/'n','t','c','n','t','c','n','t','c'/)
  ! ma1,ma2,mb1,mb2 are the sizes of the arrays A and B
  ! na1,na2,nb1,nb2 are the sizes of the arrays j(A) and j(B)
  ! where j is in {n,t,c} the transformation function
  where (ata.eq.'n')
     ama1=na1
     ama2=na2
  elsewhere
     ama1=na2
     ama2=na1
  end where
  where (atb.eq.'n')
     amb1=nb1
     amb2=nb2
  elsewhere
     amb1=nb2
     amb2=nb1
  end where
  do icase=1,maxcase
     ma1=ama1(icase); ma2=ama2(icase); mb1=amb1(icase); mb2=amb2(icase)
     ta=ata(icase); tb=atb(icase)
     allocate(za(ma1,ma2),zb(mb1,mb2),ra(ma1,ma2),rb(mb1,mb2))
     za=0.d0
     zb=0.d0
     z=0.d0
     z2=0.d0
     z3=0.d0
     call random_number(ra)
     call random_number(rb)
     za=za+ra/dble(ma1)
     zb=zb+rb/dble(mb2)
     call random_number(ra)
     call random_number(rb)
     za=za+zi*ra/dble(ma1)
     zb=zb+zi*rb/dble(mb2)

     call cpu_time(cpu0)
     call zgemm(ta,tb,na1,nb2,na2,alpha,za,ma1,zb,mb1,beta,z,nz1)
     call cpu_time(cpu1)
     cpuzgemm=cpu1-cpu0

     call cpu_time(cpu0)
     call zgemm_wrap(z2,ta,za,tb,zb,alpha,beta)
     call cpu_time(cpu1)
     cpuzgemmwrap=cpu1-cpu0

     call cpu_time(cpu0)
     if ((ta.eq.'n').and.(tb.eq.'n')) then
        z3=alpha*matmul(za,zb)+beta*z3
     else if ((ta.eq.'n').and.(tb.eq.'t')) then
        z3=alpha*matmul(za,transpose(zb))+beta*z3
     else if ((ta.eq.'n').and.(tb.eq.'c')) then
        z3=alpha*matmul(za,conjg(transpose(zb)))+beta*z3
     else if ((ta.eq.'t').and.(tb.eq.'n')) then
        z3=alpha*matmul(transpose(za),zb)+beta*z3
     else if ((ta.eq.'t').and.(tb.eq.'t')) then
        z3=alpha*matmul(transpose(za),transpose(zb))+beta*z3
     else if ((ta.eq.'t').and.(tb.eq.'c')) then
        z3=alpha*matmul(transpose(za),conjg(transpose(zb)))+beta*z3
     else if ((ta.eq.'c').and.(tb.eq.'n')) then
        z3=alpha*matmul(conjg(transpose(za)),zb)+beta*z3
     else if ((ta.eq.'c').and.(tb.eq.'t')) then
        z3=alpha*matmul(conjg(transpose(za)),transpose(zb))+beta*z3
     else if ((ta.eq.'c').and.(tb.eq.'c')) then
        z3=alpha*matmul(conjg(transpose(za)),conjg(transpose(zb)))+beta*z3
     else
        write(*,*) 'Error: see routine summations'
        stop
     end if
     call cpu_time(cpu1)
     cpumatmul=cpu1-cpu0

     write(*,*)
     write(*,'(a,i6)') 'case number: ',icase
     write(*,'(a)') 'cases : '//ta//' , '//tb
     write(*,'(a,f12.3)') 'CPU-time (zgemm)     :',cpuzgemm
     write(*,'(a,f12.3)') 'CPU-time (zgemmwrap) :',cpuzgemmwrap
     write(*,'(a,f12.3)') 'CPU-time (matmul)    :',cpumatmul
     write(*,*) 'max deviation matmul-zgemm     :',maxval(abs(z-z3))
     write(*,*) 'max deviation matmul-zgemmwrap :',maxval(abs(z2-z3))

     deallocate(za,zb,ra,rb)
  end do
  deallocate(z,z2)
end subroutine test_blaswrappers
