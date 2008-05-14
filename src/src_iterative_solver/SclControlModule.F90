! module to switch the different scl solver modes

module sclcontroll
use modmain,only:iscl,currentconvergence
implicit none
 !scl index
  integer diiscounter !! counter for DIIS iterations
  logical packedmatrixstorage
  logical tarpack,tlapack,tdiis,tjdqz
  integer,parameter:: diismax=35,diisfirstscl=3, maxdiisspace=15
  integer :: iseed(4)=1
  real(8) lowesteval
  real(8) epsarpack
  real(8) epsresid
  
  real(8) ,parameter::diisthreshould=1
  real(8) ::lastresnorm
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
    if(currentconvergence.gt.5.0e-4) then
       calculate_preconditioner =.true.
    endif
  write(*,*)"precon"
    endif
  end function calculate_preconditioner

  function doDIIScycle()
    logical doDIIScycle
    doDIIScycle=.false.
    if(tdiis) then
       !this may get more advanced:
       if(iscl.ge.diisfirstscl) doDIIScycle=.true.
       if(currentconvergence.gt.1.0)doDIIScycle=.false.
      if(doDIIScycle) write(*,*)"DIIS"
    endif
    lastresnorm=1.e10
  end function doDIIScycle



  function doprerotate_preconditioner()
    logical doprerotate_preconditioner
    doprerotate_preconditioner=.false.
    if(diiscounter.ge.3) doprerotate_preconditioner=.true.
  end function doprerotate_preconditioner

 
  function doARPACKiteration()
    logical doARPACKiteration
    doARPACKiteration=.false.
    if (tarpack) then
       doARPACKiteration=.true.
       write(*,*)"ARPACK"
       diiscounter=1
    endif
 
  end function doARPACKiteration 

  function doLAPACKsolver()
    logical doLAPACKsolver
    doLAPACKsolver=.false.
    if ((tlapack)) then
       doLAPACKsolver=.true.
       write(*,*)"LAPACK hevx"
       diiscounter=1
    endif

  end function doLAPACKsolver


  function allconverged(n,rnorms)
    logical allconverged
    integer,intent(in)::n
    real(8),intent(in):: rnorms(n)
	real(8)::rnormmax
	rnormmax=rnorms(idamax(n,rnorms,1))
    if (rnormmax.lt.epsresid) then 
       allconverged=.true.
         write(*,*)" converged",rnorms(idamax(n,rnorms,1))
    else
       allconverged=.false.
	write(*,*)"not converged",rnorms(idamax(n,rnorms,1)) ,idamax(n,rnorms,1)
    endif
    
	if(rnormmax/lastresnorm.gt.1.1)then
 		!allconverged=.true.
        write(*,*)"warning: error is gettig larger again",rnorms(idamax(n,rnorms,1))
        if (rnormmax.gt. .5e-6)then
        allconverged=.false.
        else
        write(*,*)"error:error is gettig larger again"
      !  stop
        endif
        
    endif
    lastresnorm=rnormmax       
  end function allconverged
  
  function dojacobdavidson()
  logical dojacobdavidson
  dojacobdavidson=.false.
if(tjdqz) then
dojacobdavidson=.true.
write(*,*)"JDQZ"
endif  
  end function 
end module sclcontroll
