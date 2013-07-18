!BOP
! !ROUTINE: structureoptimization
! !INTERFACE:
!
!
subroutine structureoptimization
! !USES:
    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!   Created February 2013 (DIN)
!EOP
!BOC
    Implicit None
    integer :: is, ia
    Logical :: exist

! Write first the starting configuration
    if (rank .Eq. 0) then
      if (input%structureoptimization%history) then
        call system ('rm -rf history.*')   
        call writehistory
      end if
    end if

! Use "fromfile" option during the optimization run
    task = 1
    
    if (input%structureoptimization%method=="newton") then

        if (rank .Eq. 0) Then
            write (60,*)
            write (60, '("+---------------------------------------------------------+")')
            write (60, '("| Use Newton-like method for optimizing atomic positions")')
            write (60, '("+---------------------------------------------------------+")')
            !call boxmsg(60,'-','Use Newton-like method for optimizing atomic positions')
            call flushifc(60)
        end if
        
        call newton(input%structureoptimization%epsforce)

    else if (input%structureoptimization%method=="bfgs") then

        If (rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("+---------------------------------------------------------+")')
            Write (60, '("| Use L-BFGS-B method for optimizing atomic positions")')
            Write (60, '("+---------------------------------------------------------+")')
            call flushifc(60)
        End If
        
        call lbfgs_driver

    else if (input%structureoptimization%method=="mixed") then
    
        If (rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("+---------------------------------------------------------+")')
            Write (60, '("| Use mixed Newton/BFGS scheme for optimizing atomic positions")')
            Write (60, '("+---------------------------------------------------------+")')
            call flushifc(60)
        End If
        
        call newton(input%structureoptimization%epsforce0)
        call lbfgs_driver
    
    else
        
        write(*,*)
        write(*,*) 'ERROR(structureoptimization): Unknown method for structure optimization!'
        stop
      
    end if
    
! check force convergence
    if (forcemax .Le. input%structureoptimization%epsforce) Then
        if (rank .Eq. 0) Then
            write(60,*)
            write(60,'("+---------------------------------------------------------+")')
            write(60,'("| Force convergence target achieved")')
            write(60,'("+---------------------------------------------------------+")')
        end if
    else
        if (rank .Eq. 0) Then
          write(60,*)
          write(60,'("+---------------------------------------------------------+")')
          write(60,'(" Required force convergence has not been reached!")')
          write(60,'("+---------------------------------------------------------+")')
          write(60,'(" forcemax=",f12.8," > epsforce=",f12.8)') forcemax, input%structureoptimization%epsforce
          write(60,*)
        end if
    end if

!------------------------------------!
!   SCF cycle for the final structure
!------------------------------------!
    if (rank==0) then
        write(60,*)
        write(60, '("+-----------------------------------------------------------+")')
        write(60, '("| Groundstate module started ")')
        write(60, '("+-----------------------------------------------------------+")')
    end if
    call scf_cycle(input%structureoptimization%outputlevelnumber)
    if (rank==0) then
        write(60,*)
        write(60, '("+-----------------------------------------------------------+")')
        write(60, '("| Groundstate module stopped ")')
        write(60, '("+-----------------------------------------------------------+")')
        write(60,*)
    end if

    return
end subroutine
