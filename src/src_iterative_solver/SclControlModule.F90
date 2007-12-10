! module to switch the different scl solver modes

module sclcontroll

integer iscl !scl index
integer diiscounter !! counter for DIIS iterations
integer iterativetype
integer,parameter:: diismax=5
real lowesteval
real diisthreshould,reps

contains

function calculate_preconditioner()
logical calculate_preconditioner
calculate_preconditioner=.false.
if(diiscounter.eq.1) calculate_preconditioner =.true.
end function

function doDIIScycle()
logical doDIIScycle
doDIIScycle=.false.
if(iterativetype.eq.1) then
if(iscl.ge.2) doDIIScycle=.true.
endif
end function

function doprerotate_preconditioner()
logical doprerotate_preconditioner
doprerotate_preconditioner=.false.
if(diiscounter.ge.3) doprerotate_preconditioner=.true.
end function
function doDavidsoncycle()
logical doDavidsoncycle
doDavidsoncycle=.false.
if(iterativetype.eq.3) then
if(iscl.ge.2) doDavidsoncycle=.true.
endif
end function

function doARPACKiteration()
logical doARPACKiteration
 doARPACKiteration=.false.
if (iterativetype.eq.2)  doARPACKiteration=.true.

end function

function doLAPACKsolver()
logical doLAPACKsolver
doLAPACKsolver=.false.
if (iterativetype.eq.0)  doLAPACKsolver=.true.
end function

function allconverged(n,rnorms)
logical allconverged
integer,intent(in)::n
real(8),intent(in):: rnorms(n)

if (idamax(n,rnorm,1).lt.reps) then 
allconverged=.true.
else
allconverged=.false.
endif

end function
end module