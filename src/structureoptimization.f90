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
    character*(77) :: string

!_____________________________________________________________
! Write first (if required) the starting configuration on file

    if (rank .Eq. 0) then
      if (input%structureoptimization%history) then
        call system ('rm -rf history.*')   
        call writehistory
      end if
    end if

!__________________________________________________
! Use "fromfile" option during the optimization run

    task = 1
    istep = 0

!_______________________________
! Choose the optimization method

    if (input%structureoptimization%method=="newton") then

        if (rank .Eq. 0) Then
            call printbox(60,"+","Use Newton-like method for optimizing atomic positions")
            call printbox(6, "+","Use Newton-like method for optimizing atomic positions")
            call writeoptminitdata(60,input%structureoptimization%outputlevelnumber)
            call flushifc(60)
            call flushifc(6)
        end if
        
        call newton(input%structureoptimization%epsforce)

    else if (input%structureoptimization%method=="bfgs") then

        If (rank .Eq. 0) Then
            call printbox(60,"+","Use L-BFGS-B method for optimizing atomic positions")
            call printbox(6, "+","Use L-BFGS-B method for optimizing atomic positions")
            call writeoptminitdata(60,input%structureoptimization%outputlevelnumber)
            call flushifc(60)
            call flushifc(6)
        End If
        
        call lbfgs_driver

    else if (input%structureoptimization%method=="mixed") then
    
        If (rank .Eq. 0) Then
            call printbox(60,"+","Use mixed Newton/BFGS scheme for optimizing atomic positions")
            call printbox(6, "+","Use mixed Newton/BFGS scheme for optimizing atomic positions")
            call writeoptminitdata(60,input%structureoptimization%outputlevelnumber)
            call flushifc(60)
            call flushifc(6)
        End If

        write(60,*) input%structureoptimization%epsforce0
        
        call newton(input%structureoptimization%epsforce0)
        call lbfgs_driver
    
    else
         
        call printbox(60,"#","ERROR(structureoptimization): Unknown method for structure optimization!")
        call flushifc(60)
        stop
      
    end if
    
!________________________
! check force convergence

    if (forcemax .Le. input%structureoptimization%epsforce) Then
        if (rank .Eq. 0) Then

            write(6,*)
            call flushifc(6)
            call printbox(60,"+","Force convergence target achieved")
            call flushifc(60)

        end if
    else
        if (rank .Eq. 0) Then

            write(60,*)
            call printline(60,"#")
            call printline(6,"#")
            write(string,'("Required force convergence has not been reached!")') 
            call printtext(60,"#",string)
            call printtext(6,"#",string)
            if (istep+1>input%structureoptimization%maxsteps) then
                write(string,'("Maximum number of optimization steps reached =",I5)') &
               &  input%structureoptimization%maxsteps
                call printtext(60,"#",string)
                call printtext(6,"#",string)
            end if 
            write(string,'("forcemax=",f12.8," > epsforce=",f12.8)') forcemax, input%structureoptimization%epsforce
            call printtext(60,"#",string)
            call printtext(6,"#",string)
            call printline(60,"#")
            call printline(6,"#")
            write(6,*)
            call flushifc(60)
            call flushifc(6)

        end if
    end if

!__________________________________________________
! In any case write summary of the obtained results

    Call timesec(ts0)      

    If (rank==0) Then

!___________
! output SCF

        call writeengy(60)
        Write (60,*)
        Write (60, '(" DOS at Fermi energy (states/Ha/cell)",T45 ": ", F22.12)') fermidos
        Call writechg (60,input%groundstate%outputlevelnumber)
        If (getfixspinnumber() .Ne. 0) Call writefsm (60)

!____________________
! output optimization

        Write (60,*)
        Write (60,'(" Total energy at this optimization step :",F19.9)') engytot
        call writepositions(60,2) 
        call writeforce(60,2)  

        call flushifc(60)

!_____________
! extra output

        Call scl_iter_xmlout ()
        If (associated(input%groundstate%spin)) Call scl_xml_write_moments()
        Call scl_xml_out_write()
    End If

    Call timesec(ts1)
    timeio=ts1-ts0+timeio

    return
end subroutine
!
!-----------------------------------------------------------------------------80
!
subroutine writeoptminitdata(fnum,verbosity)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
    integer, intent(in) :: verbosity
! local variables
    integer :: is, ia

    call printbox(fnum,"-","Optimization step 0: Initialize optimization")
    write(fnum,'(" Maximum force magnitude       (target) : ",F14.8,"    (", F14.8, ")")') &
   &  forcemax, input%structureoptimization%epsforce
    write(fnum,'(" Total energy at this optimization step :",F19.9)') engytot
    call flushifc(fnum)
    if (verbosity < 1) return

    call writepositions(fnum,verbosity) 
    call writeforce(fnum,verbosity)
    
    return
end subroutine
