subroutine test_array_assert_equal()
use modreport
implicit none
integer, parameter::n=100
real(8)::darray1(n),darray2(n)
complex(8)::zarray1(n),zarray2(n)
real(8)::tol
logical passed,eq1,eq2
character(256)::name
darray1=0
darray1(2)=1.3e-8
darray2=0
zarray1=0
zarray2=0
zarray2(7)=(1.3e-8,1.3e-8)
tol=1e-7
write(name,*)"darray.err"
call darray_assert_equal(name,darray1, darray2, n,tol,eq1)
write(name,*)"darray.err"
darray1(50)=3e-6
call darray_assert_equal(name,darray1, darray2, n,tol,eq2)
testunitname="darray_assert_equal"
inputf=""
outputf=trim(name)
passed=eq1 .and. (.not. eq2)
call testreport(passed)

write(name,*)"zarray.err"
call zarray_assert_equal(name,zarray1, zarray2, n,tol,passed)
testunitname="zarray_assert_equal"
inputf=""
outputf=trim(name)
call testreport(passed)
end subroutine
