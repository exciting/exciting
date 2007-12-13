! module to switch the different scl solver modes

module sclcontroll
use modmain,only:iscl
implicit none
 !scl index
  integer diiscounter !! counter for DIIS iterations
  integer iterativetype
  integer,parameter:: diismax=5
  real lowesteval
  real diisthreshould,reps
  real idamax
  external idamax
contains

  function calculate_preconditioner()
    logical calculate_preconditioner
    calculate_preconditioner=.false.
    if(diiscounter.eq.1) then
       calculate_preconditioner =.true.
  write(*,*)"precond"
    else
      
    endif
  end function calculate_preconditioner

  function doDIIScycle()
    logical doDIIScycle
    doDIIScycle=.false.
    if(iterativetype.eq.1) then
       !this may get more advanced:
       if(iscl.ge.2) doDIIScycle=.true.
       write(*,*)"DIIS"
    endif
  end function doDIIScycle

  function prediis()
    logical prediis
    prediis=.false.
    if (iterativetype.eq.1.and.(iscl.lt.2))  prediis=.true.
  end function prediis

  function doprerotate_preconditioner()
    logical doprerotate_preconditioner
    doprerotate_preconditioner=.false.
    if(diiscounter.ge.3) doprerotate_preconditioner=.true.
  end function doprerotate_preconditioner

  function doDavidsoncycle()
    logical doDavidsoncycle
    doDavidsoncycle=.false.
    if(iterativetype.eq.3) then
       if(iscl.ge.2) doDavidsoncycle=.true.
    endif
  end function doDavidsoncycle

  function doARPACKiteration()
    logical doARPACKiteration
    doARPACKiteration=.false.
    if (iterativetype.eq.2.or.prediis()) then
       doARPACKiteration=.true.
       write(*,*)"ARPACK"
       diiscounter=1
    endif
 
  end function doARPACKiteration

  function doLAPACKsolver()
    logical doLAPACKsolver
    doLAPACKsolver=.false.
    if ((iterativetype.eq.0).or.prediis()) then
       doLAPACKsolver=.true.
       write(*,*)"LAPACK hevx"
       diiscounter=1
    endif

  end function doLAPACKsolver


  function allconverged(n,rnorms)
    logical allconverged
    integer,intent(in)::n
    real(8),intent(in):: rnorms(n)

    if (idamax(n,rnorms,1).lt.reps) then 
       allconverged=.true.
    else
       allconverged=.false.
    endif

  end function allconverged
end module sclcontroll
