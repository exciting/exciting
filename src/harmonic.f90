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
        logical :: l1, l2
        character*(77) :: string, rmethod(3), newtonflag

        nstep = 0
        newtonflag = "full newton step"

!_______________________
! start relaxation steps

        do while ((forcemax>forcetol).and.(istep<input%relax%maxsteps))

            nstep = nstep+1
            istep = istep+1
            if (rank .Eq. 0) then

                if (lstep) then 
                    istep = istep-1
                    lstep = .False.
                else
                    write(string,'("Optimization step ", I4,"    (method = harmonic)")') istep
                    call printbox(60,"-",string)
                    call flushifc(60)
                end if

!_____________________________________________________________
! write lattice vectors and optimised atomic positions to file

                Call writehistory
                Call writegeometryxml(.True.)

!__________________________________________________
! write the optimised interatomic distances to file

                Call writeiad(.True.)
            end if

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
                write(60,'(" Number of scf iterations               : ", I5, "    (",A,")")') &
               &  iscl, trim(newtonflag)
                write(60,'(" Maximum force magnitude       (target) : ",F14.8,"    (", F14.8, ")")') &
               &  forcemax, input%relax%epsforce
                write(60,'(" Total energy at this optimization step :",F19.9)') engytot
                if (input%relax%outputlevelnumber>0) then 
                    call writepositions(60,input%relax%outputlevelnumber) 
                    call writeforce(60,input%relax%outputlevelnumber)
                end if
                if (input%relax%outputlevelnumber>1)  then 
                    call writechg (60,input%relax%outputlevelnumber)          
                end if
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

!   Updates the current atomic positions according to the force on each atom. If
!   ${\bf r}_{ij}^m$ is the position and ${\bf F}_{ij}^m$ is the force acting on
!   it for atom $j$ of species $i$ and after time step $m$, then the new
!   position is calculated by an harmonic approximation

        Implicit None
        Integer :: ik, ispn, is, ia, ias, i
        Real (8) :: t1, beta, delta(3), zero(3), epsharmonic, tstep, hdelta, xdelta, a
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

!____________________________________________
! for the first step always use Newton method

                if (nstep<2) then
                    do i = 1, 3
                        l1 = .not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)
                        l2 = ( abs(forcetot(i,ias)) .gt. input%structure%epslat/10.d0 )
                        if ( l1 .and. l2 ) then
                            xdelta = input%relax%tau0atm*forcetot(i,ias)
                            if (xdelta .gt. input%relax%tau0atm/2.0) xdelta = input%relax%tau0atm/2.0
                            atposc_1(i,ia,is) = atposc(i,ia,is)
                            atposc(i,ia,is) = atposc(i,ia,is) + xdelta      

                   !write(60,'(T24,4F14.8)') xdelta*a, xdelta, input%relax%tau0atm/2.0, tauxyz(i,ias)
                   !call flushifc(60)

                        end if
                    end do
                    tauatm(ias) = 0.0

                   !write(60,*) 
                   !write(60,'("0=> 0    ", 6F14.8,I4)') atposc(:,ia,is)*a,   forcetot(:,ias), ia
                   !write(60,'("0=> 1    ", 6F14.8,I4)') atposc_1(:,ia,is)*a, forcetp(:,ias),  ia
                   !call flushifc(60)

!____________________________________________
! otherwise use harmonic method

                else

                   !write(60,*) 
                   !write(60,'("A=> 0    ", 6F14.8,I4)') atposc(:,ia,is)*a,   forcetot(:,ias), ia
                   !write(60,'("A=> 1    ", 6F14.8,I4)') atposc_1(:,ia,is)*a, forcetp(:,ias),  ia
                   !call flushifc(60)

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

                   !write(60,'(T52,2F14.8,2I2)') beta, hdelta*a, ias, i 
                   !call flushifc(60)

!_________________________________________________________________________________________________
! if beta is positive and less then one or if the absolute value of difference between the old and 
! the new cartesian coordinates is larger then 2.0*input%relax%taubfgs, update using newton step 
! otherwise using harmonic step

                            if ( ((beta .ge. 0.0) .and. (beta .lt. 1.0)) .or. &
                           &  (abs(hdelta) .ge. input%relax%taubfgs) ) then
 
                                rmethod(i) = "n"
                                checknewton = .True. 
                                tauxyz(i,ias) = tauxyz(i,ias) + 2.0*input%relax%tau0atm
                                xdelta = tauxyz(i,ias)*forcetot(i,ias)

                   !write(60,'(T24,4F14.8,A)') xdelta, tauxyz(i,ias), forcetot(i,ias), input%relax%tau0atm/2.0, " qui"

                                if (abs(xdelta) .gt. input%relax%tau0atm) then
                                    if (xdelta < 0.0) then 
                                        xdelta = -input%relax%tau0atm
                                    else
                                        xdelta = input%relax%tau0atm
                                    end if
                                end if

                   !write(60,'(T24,4F14.8)') xdelta*a, xdelta, input%relax%tau0atm/2.0, tauxyz(i,ias)
                   !call flushifc(60)

                            else 
                                rmethod(i) = "h"
                                checkharmonic = .True. 
                                xdelta = hdelta
                                tauxyz(i,ias) = 0.0
                   
                   !write(60,'(T24,A,F14.8)') "sono passato da harmonic  ",tauxyz(i,ias)
                   !call flushifc(60)

                            end if

                            atposc_1(i,ia,is) = atposc(i,ia,is)
                            atposc(i,ia,is) = atposc(i,ia,is) + xdelta       

                   !write(60,*) "===> ", ias, i, tauxyz(i,ias), trim(rmethod(i))
                   !call flushifc(60)


                        end if
                    end do
 
                   !write(60,*) 
                   !write(60,'("B=> 0    ", 6F14.8,I4)') atposc(:,ia,is)*a,   forcetot(:,ias), ia
                   !write(60,'("B=> 1    ", 6F14.8,I4)') atposc_1(:,ia,is)*a, forcetp(:,ias),  ia
                   !call flushifc(60)

                 end if
             End Do
        End Do
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

!__________________________________________________________
! if tshift .eq. true calculate the shift of the first atom

        if (input%structure%tshift) delta(:) = atposc(:,1,1) - atposc_1(:,1,1)

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)

                if (input%structure%tshift) then
                    atposc(:,ia,is) = atposc(:,ia,is) - delta(:)
                    atposc_1(:,ia,is) = atposc_1(:,ia,is) - delta(:)
                end if

!_____________________________________________________
! set the total forces for the previous configurations

                forcetp(:,ias) = forcetot(:,ias)

                   !write(60,*) 
                   !write(60,'("delta    ", 3F14.8)') delta(:)*a
                   !write(60,*) 
                   !write(60,'("C=> 0    ", 6F14.8,I4)') atposc(:,ia,is)*a,   forcetot(:,ias), ia
                   !write(60,'("C=> 1    ", 6F14.8,I4)') atposc_1(:,ia,is)*a, forcetp(:,ias),  ia
                   !call flushifc(60)

            End Do
        End Do

!_______________________
! restart initialization

        Call init_relax

        !Do is = 1, nspecies
        !    Do ia = 1, natoms (is)
        !        ias = idxas (ia, is)
        !           write(60,*) 
        !           write(60,'("D=> 0    ", 6F14.8,I4)') atposc(:,ia,is)*a,   forcetot(:,ias), ia
        !           write(60,'("D=> 1    ", 6F14.8,I4)') atposc_1(:,ia,is)*a, forcetp(:,ias),  ia
        !           call flushifc(60)
        !    End Do
        !End Do

        Return
    End Subroutine updatpos_harmonic
    
end subroutine harmonic
