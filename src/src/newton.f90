!
!-----------------------------------------------------------------------------80.
! Copyright (C) 2013- exciting team
! This file is distributed under the terms of the GNU General Public License.
! Created       on 15-10-2013 Dimtrii Nabok (exciting team)
! Modified      on 15-11-2013 Pasquale Pavone (exciting team)
! Last modified on 14-06-2014 Pasquale Pavone (exciting team)
!-----------------------------------------------------------------------------80
!
    subroutine newton(forcetol)

!--------------------------------------------------------------
! Simple (Newton-like) method for atomic position optimization
!--------------------------------------------------------------

        use modmain
        use modinput
        use modmpi, only: rank
        
        implicit none
        real(8), intent(IN) :: forcetol
        integer :: is, ia, nstep, ias
        character*(77) :: string
        logical :: inittime

        nstep = 0

!_______________________
! start relaxation steps

        do while ((forcemax>forcetol).and.(istep<input%relax%maxsteps))

            nstep = nstep+1
            istep = istep+1

            if (lstep) then 
                istep = istep-1
                lstep = .False.
                inittime = .False.
                do is = 1, nspecies
                    do ia = 1, natoms (is)
                        ias = idxas (ia, is)
                        tauatm(ias) = 0.d0
                        forcetp(:, ias) = 0.d0
                    end do
                end do
            else
                if (rank .Eq. 0) then
                    write(string,'("Optimization step ", I4,"    (method = newton)")') istep
                    call printbox(60,"-",string)
                    call flushifc(60)
                end if
                inittime = .True.
            end if

            if (inittime) call timesec(tsec1)

!________________________
! update atomic positions

            Call updatpos

!_______________________
! restart initialization

            Call init_relax

!______________________________________________________
! begin new self-consistent loop with updated positions

            call scf_cycle(-1)

!____________
! output info

            if (rank .Eq. 0) then
                write(60,'(" Number of scf iterations",T45,": ", I9)') iscl
                write(60,'(" Maximum force magnitude",T36,"(target) : ",F18.8,"  (", F12.8, ")")') &
               &  forcemax, input%relax%epsforce
                write(60,'(" Total energy at this optimization step",T45,": ",F18.8)') engytot
                if (input%relax%outputlevelnumber>0) then 
                    call writepositions(60,input%relax%outputlevelnumber) 
                    call writeforce(60,input%relax%outputlevelnumber)
                end if
                if (input%relax%outputlevelnumber>1) call writechg (60,input%relax%outputlevelnumber)          
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
    subroutine updatpos
        use modinput
        use modmain

!   Updates the current atomic positions according to the force on each atom. If
!   ${\bf r}_{ij}^m$ is the position and ${\bf F}_{ij}^m$ is the force acting on
!   it for atom $j$ of species $i$ and after time step $m$, then the new
!   position is calculated by
!   $$ {\bf r}_{ij}^{m+1}={\bf r}_{ij}^m+\tau_{ij}^m\left({\bf F}_{ij}^m
!    +{\bf F}_{ij}^{m-1}\right), $$
!   where $\tau_{ij}^m$ is a parameter governing the size of the displacement.
!   If ${\bf F}_{ij}^m\cdot{\bf F}_{ij}^{m-1}>0$ then $\tau_{ij}^m$ is
!   increased, otherwise it is decreased.

        Implicit None
        Integer :: ik, ispn, is, ia, ias, i
        Real (8) :: t1, xdelta

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)

!_____________________________________________________________________
! compute the dot-product between the current and previous total force

                t1 = dot_product (forcetot(:, ias), forcetp(:, ias))

!________________________________________________________________________
! if the force is in the same direction then increase step size parameter

                If (t1 .Gt. 0.d0) Then
                    tauatm(ias) = tauatm(ias)+input%relax%taunewton
                Else
                    tauatm(ias) = input%relax%taunewton
                End If
                do i = 1, 3
                    if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)) then
                        xdelta = tauatm(ias)*(forcetot(i,ias)+forcetp(i,ias))
                        if (abs(xdelta) .gt. input%relax%taunewton) then
                            if (xdelta < 0.0) then 
                                xdelta = -input%relax%taunewton
                            else
                                xdelta = input%relax%taunewton
                            end if
                        end if
                        atposc(i,ia,is) = atposc(i,ia,is) + xdelta
                    end if 
                end do 
            End Do
        End Do

!________________________________________________________
! compute the lattice coordinates of the atomic positions

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)

                Call r3mv(ainv,atposc(:,ia,is),input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))

!____________________________________________
! set the previous to the current total force

                forcetp(:, ias) = forcetot(:, ias)
            End Do
        End Do
      
        Return
    End Subroutine updatpos
    
end subroutine newton
