!
!-----------------------------------------------------------------------------80.
! Copyright (C) 2013- exciting team
! This file is distributed under the terms of the GNU General Public License.
! Created       on 15-10-2013 Pasquale Pavone (exciting team)
! Modified      on 15-11-2013 Pasquale Pavone (exciting team)
! Last modified on 14-06-2014 Pasquale Pavone (exciting team)
!-----------------------------------------------------------------------------80
!
    subroutine harmonic(forcetol)

!--------------------------------------------------------------
! Harmonic method for atomic position optimization
!--------------------------------------------------------------

        use modmain
        use modinput
        use modmpi, only: rank
        
        implicit none
        real(8), intent(IN) :: forcetol
        integer :: is, ia, nstep, ias
        logical :: l1, l2, inittime
        character*(77) :: string, rmethod(3), newtonflag

        nstep = 0
        newtonflag = "full newton step"

!_______________________
! start relaxation steps

        do while ((forcemax>forcetol).and.(istep<input%relax%maxsteps))

            nstep = nstep+1
            istep = istep+1

            if (lstep) then 
                istep = istep-1
                lstep = .False.
                inittime = .False.
            else
                if (rank .Eq. 0) then
                    write(string,'("Optimization step ", I4,"    (method = harmonic)")') istep
                    call printbox(60,"-",string)
                    call flushifc(60)
                end if
                inittime = .True.
            end if

            if (inittime) call timesec(tsec1)

!________________________________________________________________
! at the first step initialize atomic possitions and total forces

            if (nstep<2) then
                forcetp(:, :) = forcetot(:, :)
            end if

!________________________
! update atomic positions

            Call updatpos_harmonic

!______________________________________________________
! begin new self-consistent loop with updated positions

            call scf_cycle(-1)

!____________
! output info

            if (rank .Eq. 0) then
                write(60,'(" Number of scf iterations",T45,": ", I9, "  (",A,")")') &
               &  iscl, trim(newtonflag)
                write(60,'(" Maximum force magnitude",T36,"(target) : ",F18.8,"  (", F12.8, ")")') &
               &  forcemax, input%relax%epsforce
                write(60,'(" Total energy at this optimization step",T45,": ",F18.8)') engytot
                if (input%relax%outputlevelnumber>0) then 
                    call writepositions(60,input%relax%outputlevelnumber) 
                    call writeforce(60,input%relax%outputlevelnumber)
                end if
                if (input%relax%outputlevelnumber>1) call writechg(60,input%relax%outputlevelnumber)          
                if (input%relax%printtorque) call writetorque(60)          
                call flushifc(60)

!_____________________________________________________________
! write lattice vectors and optimised atomic positions to file

                Call writehistory
                Call writegeometryxml(.True.)

!__________________________________________________
! write the optimised interatomic distances to file

                Call writeiad(.True.)
            end if

!_______________________________________________
! write the time spent in this optimization step 

            call timesec(tsec2)
            if (rank==0) then
                write(60,*)
                write(60,'(" Time spent in this optimization step",T45,": ",F12.2," seconds")') tsec2-tsec1
                call flushifc(60)
            end if

        end do
!
contains
!    
!-----------------------------------------------------------------------------80
!
    subroutine updatpos_harmonic
        use modinput
        use modmain

!       Updates the current atomic positions according to the force on each atom. 
!       The new configuration is calculated using either a newton-like or a harmonic step.

        Implicit None
        Integer :: ik, ispn, is, ia, ias, i
        Real (8) :: t1, beta, delta(3), zero(3), epsharmonic, tstep, hdelta, xdelta, a,dumdum
        Logical :: checknewton, checkharmonic

        a = 1./input%structure%crystal%scale

        zero(:) = 0.d0
        epsharmonic = 1.0d-8
        checknewton = .False. 
        checkharmonic = .False.

!________________________
! update atomic positions

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)

!_________________________________________________
! for the first step always use a newton-like step

                if (nstep<2) then
                    do i = 1, 3
                        l1 = .not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)
                        l2 = ( abs(forcetot(i,ias)) .gt. min(input%structure%epslat/10.d0,input%relax%epsforce) )
                        if ( l1 .and. l2 ) then
                            xdelta = input%relax%taunewton*forcetot(i,ias)
                            if (xdelta .gt. input%relax%taunewton/2.0) xdelta = input%relax%taunewton/2.0
                            atposc_1(i,ia,is) = atposc(i,ia,is)
                            atposc(i,ia,is) = atposc(i,ia,is) + xdelta   

                   !write(60, '(A,4X,F12.8,3I5)') "n",xdelta,i,ia,is
                   !call flushifc(60)
        
                        end if
                    end do
                    tauatm(ias) = 0.0

!______________________________
! otherwise use harmonic method

                else

                    do i = 1, 3
                        rmethod(i) = "n"
                        l1 = .not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)

!________________________________________________________________
! forces with amplitude below epsharmonic are considered as zeros

                        l2 = ( abs(forcetot(i,ias)) .gt. epsharmonic )
                        if ( l1 .and. l2 ) then 
                            beta = forcetp(i,ias)/forcetot(i,ias)
                            hdelta = atposc_1(i,ia,is)/(1.-beta) - beta*atposc(i,ia,is)/(1.-beta)
                            hdelta = hdelta - atposc(i,ia,is)

!_________________________________________________________________________________________________
! if beta is positive and less then one or if the absolute value of difference between the old and 
! the new cartesian coordinates is larger then 2.0*input%relax%taubfgs, update using a newton-like 
! step, otherwise using a harmonic step

                            if ( (beta .lt. 0.0) .or. (beta .ge. 3) ) then
 
!______________
! harmonic step

                                rmethod(i) = "h"
                                checkharmonic = .True. 
                                xdelta = hdelta                                
                                tauxyz(i,ias) = 0.0

                   !write(60, '(A,4X,F12.8,3I5)') "h",xdelta,i,ia,is
                   !call flushifc(60)

                            else 

!_________________
! newton-like step

                                rmethod(i) = "n"
                                checknewton = .True. 
                                tauxyz(i,ias) = tauxyz(i,ias) + 2.0*input%relax%taunewton
                                xdelta = tauxyz(i,ias)*forcetot(i,ias)

                   !write(60, '(A,4X,F12.8,3I5)') "n",xdelta,i,ia,is
                   !call flushifc(60)
                            
                            end if

!_____________________________
! update old and new positions

                            atposc_1(i,ia,is) = atposc(i,ia,is)
                            atposc(i,ia,is) = atposc(i,ia,is) + xdelta       

                        end if
                    end do
                 end if
             End Do
        End Do

!________________________
! set flags for each step

        newtonflag = "full newton step"
        if (checkharmonic .and. checknewton) then 
            newtonflag = "partial harmonic step"
        else
            if (checkharmonic .and. (.not.checknewton)) then
                newtonflag = "full harmonic step"
            end if
        end if

!________________________________________________________
! compute the lattice coordinates of the atomic positions

        Do is = 1, nspecies
            Do ia = 1, natoms (is)

                Call r3mv(ainv,atposc(:,ia,is),input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))

            End Do
        End Do

!_________________________________________________________
! set the total forces for the previous configurations and
! save the current atomic position in atsave

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)
                forcetp(:,ias) = forcetot(:,ias)
                atsave(:,ia,is) = atposc(:,ia,is)
            End Do
        End Do

!_______________________
! restart initialization

        Call init_relax

!______________________________________________________
! if (tshift .eq. true) shift also old atomic positions 

        if (input%structure%tshift) then 
            do is = 1, nspecies
                do ia = 1, natoms (is)
                   delta(:) = atposc(:,ia,is) - atsave(:,ia,is)
                   atposc_1(:,ia,is) = atposc_1(:,ia,is) + delta(:)
                end do
            end do
        end if

        Return
    End Subroutine updatpos_harmonic
    
end subroutine harmonic
