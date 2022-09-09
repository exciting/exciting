subroutine darray_assert_equal(name,array1, array2, n,tol,equal)
implicit none
character(256)::name
real(8),intent(in)::array1(*), array2(*)
integer,intent(in)::n
real(8),intent(in)::tol
logical, intent(out)::equal
real(8),parameter::minusone=-1.0
!local:

real(8),allocatable::errarray(:)
real(8)::errmax,erraverage
integer,external:: idamax
real(8),external::dasum
allocate (errarray(n))

call dcopy(n,array1,1,errarray,1)
call daxpy(n,minusone,array2,1,errarray,1)
erraverage=dasum(n,errarray,1)/n
errmax=abs(errarray(IDAMAX(n,errarray,1)))

if (errmax .gt. tol)then
 equal=.false.
 else
 equal=.true.
 endif
call write_arraycompare_statistics(name, erraverage, errmax,tol)
deallocate (errarray)
return
end subroutine

subroutine zarray_assert_equal(name,array1, array2, n,tol,equal)
implicit none
character(256)::name
complex(8),intent(in)::array1(*), array2(*)
integer,intent(in)::n
real(8),intent(in)::tol
logical, intent(out)::equal
complex(8),parameter :: zminusone=(-1,0)
!local:

complex(8),allocatable::errarray(:)
real(8)::errmax,erraverage
integer,external:: izamax
real(8),external::dzasum
allocate (errarray(n))

call zcopy(n,array1,1,errarray,1)
call zaxpy(n,zminusone,array2,1,errarray,1)
erraverage=dzasum(n,errarray,1)/n
errmax=abs(errarray(IZAMAX(n,errarray,1)))
if (errmax .gt. tol)then
 equal=.false.
 else
 equal=.true.
 endif
call write_arraycompare_statistics(name, erraverage, errmax,tol)
deallocate (errarray)
return
end subroutine


subroutine write_arraycompare_statistics(name, erraverage, errmax,tol)
implicit none
character(256),intent(in)::name
real(8),intent(in):: erraverage, errmax,tol
character(256)::tmp,filename

write (tmp,*)"../output/",trim(adjustl (name))
filename=trim(adjustl(tmp))
open(35,file=filename,action="WRITE",form='FORMATTED')
write(35,*) "errmax        erraverage    tol  "
write(35,*) errmax,erraverage,tol
close(35)
end subroutine
