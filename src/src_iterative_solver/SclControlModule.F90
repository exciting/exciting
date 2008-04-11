! module to switch the different scl solver modes

module sclcontroll
use modmain,only:iscl,currentconvergence
implicit none
 !scl index
  integer diiscounter !! counter for DIIS iterations
  integer iterativetype
  logical packedmatrixstorage
  integer,parameter:: diismax=25,diisfirstscl=3
  real(8) lowesteval
  real(8) epsarpack
  real ,parameter::diisthreshould=1,reps=0.5e-7
  integer ,parameter::jacofidavidsonfirstscl=1
integer idamax
  external idamax
 logical  recalculate_preconditioner
contains
 
  function calculate_preconditioner()
    logical calculate_preconditioner
    calculate_preconditioner=.false.
    if(diiscounter.eq.1) then
       calculate_preconditioner =.true.
    
  write(*,*)"precond"
    else if (mod(diiscounter,5).eq.0)then
    if(currentconvergence.gt.1e-3) then
       calculate_preconditioner =.true.
    endif
  write(*,*)"precon"
    endif
  end function calculate_preconditioner

  function doDIIScycle()
    logical doDIIScycle
    doDIIScycle=.false.
    if(iterativetype.eq.1) then
       !this may get more advanced:
       if(iscl.ge.diisfirstscl) doDIIScycle=.true.
       if(currentconvergence.gt.1.0)doDIIScycle=.false.
      if(doDIIScycle) write(*,*)"DIIS"
    endif
  end function doDIIScycle



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
    if (iterativetype.ge.2) then
       doARPACKiteration=.true.
       write(*,*)"ARPACK"
       diiscounter=1
    endif
 
  end function doARPACKiteration 

  function doLAPACKsolver()
    logical doLAPACKsolver
    doLAPACKsolver=.false.
    if ((iterativetype.eq.0.or.iterativetype.eq.1)) then
       doLAPACKsolver=.true.
       write(*,*)"LAPACK hevx"
       diiscounter=1
    endif

  end function doLAPACKsolver


  function allconverged(n,rnorms)
    logical allconverged
    integer,intent(in)::n
    real(8),intent(in):: rnorms(n)

    if (rnorms(idamax(n,rnorms,1)).lt.reps) then 
       allconverged=.true.
         write(*,*)" converged",rnorms(idamax(n,rnorms,1))
    else
       allconverged=.false.
	write(*,*)"not converged",rnorms(idamax(n,rnorms,1)) ,idamax(n,rnorms,1)
    endif

  end function allconverged
  
  function dojacobdavidson()
  logical dojacobdavidson
  dojacobdavidson=.false.
if(iterativetype.eq.3.and.iscl.ge.jacofidavidsonfirstscl) then
dojacobdavidson=.true.
write(*,*)"JDQZ"
endif  
  end function 
end module sclcontroll
