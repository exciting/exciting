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
    real(8) :: ts0, ts1
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
!   SCF data
!------------------------------------!
    Call timesec(ts0)      
    If ((input%structureoptimization%outputlevelnumber>0).and.(rank==0)) Then
! output energy components
        call writeengy(60)
        Write (60,*)
        Write (60, '("Density of states at Fermi energy : ", G18.10)') fermidos
        Write (60, '(" (states/Hartree/unit cell)")')
! output charges and moments
        Call writechg (60)
! output effective fields for fixed spin moment calculations
        If (getfixspinnumber() .Ne. 0) Call writefsm (60)
! output forces to INFO.OUT
        call writeforce(60,input%structureoptimization%outputlevelnumber)
        Call scl_iter_xmlout ()
        If (associated(input%groundstate%spin)) Call scl_xml_write_moments()
        Call scl_xml_out_write()
    End If
    Call timesec(ts1)
    timeio=ts1-ts0+timeio

    return
end subroutine
