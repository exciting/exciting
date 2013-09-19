!BOP
! !ROUTINE: relax
! !INTERFACE:
!
!
subroutine relax
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
    integer :: is, ia, idm
    real(8) :: ts0, ts1
    Logical :: exist
    character*(77) :: string

!_________________________________________________________________________________________________
! initialize step sizes, the previous forces, and atomic position in cartesian coordinates

    If (allocated(forcetp)) deallocate (forcetp)
    Allocate (forcetp(3, natmtot))
    If (allocated(tauatm)) deallocate (tauatm)
    Allocate (tauatm(natmtot))
    If (allocated(tauxyz)) deallocate (tauxyz)
    Allocate (tauxyz(3, natmtot))
    If (associated(input%relax)) Then
        tauatm (:) = input%relax%tau0atm
        tauxyz (:, :) = input%relax%tau0atm
    Else
        tauatm (:) = 0
        tauxyz (:, :) = 0
    End If
    forcetp (:, :) = 0.d0
    atposc_1(:,:,:) = atposc(:,:,:)

!_____________________________________________________________
! Write first (if required) the starting configuration on file

    if (rank .Eq. 0) then
      if ((input%relax%history).and. (input%relax%addtohistory)) then
        call writehistory
      elseif (input%relax%history) then
        select case (trim(input%relax%historyformat))
        case('xyz','XYZ')
             call system ('mv history.xyz history.xyz.backup >& /dev/null')
        case('gulp','GULP')
             call system ('mv history.gin history.gin.backup >& /dev/null')
        case default
            write(*,*)'ERROR(relax): Unknown output format'
            stop
        end select
!       call system ('cp history.* backup.history')
        call writehistory
      end if
    end if

!    if (rank .Eq. 0) then
!      if (input%relax%history) then
!        call system ('rm -rf history.*')   
!        call writehistory
!      end if
!    end if

!__________________________________________________
! Use "fromfile" option during the optimization run

    task = 1
    istep = 0
    lstep = .False.
!    call readstate

!_______________________________
! Choose the optimization method

    if (input%relax%method=="newton") then

        if (rank .Eq. 0) Then
            call printbox(60,"+","Use Newton-like method for optimizing atomic positions")
            call writeoptminitdata(60,input%relax%outputlevelnumber)
            call flushifc(60)
        end if
        
        call newton(input%relax%epsforce)

    else if (input%relax%method=="harmonic") then

        If (rank .Eq. 0) Then
            call printbox(60,"+","Use harmonic method for optimizing atomic positions")
            call writeoptminitdata(60,input%relax%outputlevelnumber)
            call flushifc(60)
        End If
        
        call harmonic(input%relax%epsforce)

    else if (input%relax%method=="bfgs") then

        If (rank .Eq. 0) Then
            call printbox(60,"+","Use L-BFGS-B method for optimizing atomic positions")
            call writeoptminitdata(60,input%relax%outputlevelnumber)
            call flushifc(60)
        End If
        
        call lbfgs_driver

    else if (input%relax%method=="mixed") then
    
        If (rank .Eq. 0) Then
            call printbox(60,"+","Use mixed Newton/BFGS scheme for optimizing atomic positions")
            call writeoptminitdata(60,input%relax%outputlevelnumber)
            call flushifc(60)
        End If
        
        call newton(input%relax%epsforce0)
        call lbfgs_driver
    
    else

        call printbox(60,"#","ERROR(relax): Unknown method for structure optimization!")
        call flushifc(60)
        stop
      
    end if
    
!________________________
! check force convergence

    if (forcemax .Le. input%relax%epsforce) Then
        if (rank .Eq. 0) Then

            call printbox(60,"+","Force convergence target achieved")
            call flushifc(60)

        end if
    else
        if (rank .Eq. 0) Then

            write(6,*)
            write(60,*)
            call printline(60,"#")
            call printline(6,"#")
            write(string,'("Required force convergence has not been reached!")') 
            call printtext(60,"#",string)
            call printtext(6,"#",string)
            if (istep+1>input%relax%maxsteps) then
                write(string,'("Maximum number of optimization steps reached =",I5)') &
               &  input%relax%maxsteps
                call printtext(60,"#",string)
                call printtext(6,"#",string)
            end if 
            write(string,'("forcemax=",f12.8," > epsforce=",f12.8)') forcemax, input%relax%epsforce
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

        !Write (60,*)
        !Write (60,'(" Total energy at this optimization step :",F19.9)') engytot
        idm = 0
        write(60,*)
        write(60,'(" Optimized atomic positions (lattice) :")')
        do is = 1, nspecies
            do ia = 1, natoms (is)
                idm = idm+1
                write(60,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
               &  idm, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            end do
        end do
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
   &  forcemax, input%relax%epsforce
    write(fnum,'(" Total energy at this optimization step :",F19.9)') engytot
    call flushifc(fnum)
    if (verbosity < 1) return

    call writepositions(fnum,verbosity) 
    call writeforce(fnum,verbosity)
    if (verbosity > 1) call writechg(fnum,verbosity)
    
    return
end subroutine
