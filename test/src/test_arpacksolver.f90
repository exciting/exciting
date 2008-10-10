
subroutine test_arpacksolver ()
 use modfvsystem
 use modmain
 use sclcontroll
use modreport
implicit none
integer,parameter:: ik=1,ispn=1
type (evsystem)::system
complex(8)::evecfv(nmat(ik,ispn),nstfv) ,evecfvref(nmat(ik,ispn),nstfv)
real(8)::evalfv(nstfv),evalfvref(nstfv)
real(8):: tol=1e-7
logical::equalevec,equaleval,packed
character(256)::name
  packed=.true.


call newsystem(system,packed,nmat(ik,ispn))
!call evSystemRestorefromFile(system)

!call arpacksolve(system,evalfv,evecfv)
call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvref)
call getevalfv(vkl(1,ik),evalfvref)
testunitname="arpackvectors"
inputf="EVECFV.OUT"
outputf="eigenvectors.err"
name=outputf

call zarray_assert_equal(name,evecfv, evecfvref, nmat(ik,ispn)*nstfv,tol,equalevec)
call testreport(equalevec)

testunitname="arpackvalues"
inputf="EIGVEC.OUT"

outputf="eigenvalues.err"
name=outputf

call darray_assert_equal(name,evalfv,evalfvref, nstfv,tol,equaleval)


call testreport(equaleval)
call deleteystem(system)

end subroutine
