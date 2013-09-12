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

        nstep = 0

        do while ((forcemax>forcetol).and.(istep<input%relax%maxsteps))

            nstep = nstep+1
            istep = istep+1

            if (rank .Eq. 0) then

                if (lstep) then 
                    istep = istep-1
                    lstep = .False.
                    do is = 1, nspecies
                        do ia = 1, natoms (is)
                            ias = idxas (ia, is)
                            tauatm(ias) = 0.d0
                            forcetp(:, ias) = 0.d0
                        end do
                    end do
                else
                    write(string,'("Optimization step ", I4,"    (method = newton)")') istep
                    call printbox(60,"-",string)
                    call flushifc(60)
                end if

!_____________________________________________________________
! write lattice vectors and optimised atomic positions to file

!                Call writehistory(istep)
                Call writehistory
                Call writegeometryxml(.True.)

!__________________________________________________
! write the optimised interatomic distances to file

                Call writeiad(.True.)
            end if

!________________________
! update atomic positions

            Call updatpos

!______________________________________________________
! begin new self-consistent loop with updated positions

            call scf_cycle(-1)

!____________
! output info

            if (rank .Eq. 0) then
                write(60,'(" Number of scf iterations               : ", I5)') iscl
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
        Real (8) :: t1

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)

!_____________________________________________________________________
! compute the dot-product between the current and previous total force

                t1 = dot_product (forcetot(:, ias), forcetp(:, ias))

!________________________________________________________________________
! if the force is in the same direction then increase step size parameter

                If (t1 .Gt. 0.d0) Then
                    tauatm(ias) = tauatm(ias)+input%relax%tau0atm
                Else
                    tauatm(ias) = input%relax%tau0atm
                End If
                do i = 1, 3
                    if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)) &
                   &   atposc(i,ia,is) = atposc(i,ia,is) + tauatm(ias)*(forcetot(i,ias)+forcetp(i,ias))
                end do ! i
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

!___________________________________________________________________
! find the crystal symmetries and shift atomic positions if required

        Call findsymcrys

!__________________________________
! check for overlapping muffin-tins

        Call checkmt

!_________________________________________
! generate structure factors for G-vectors

        Call gensfacgp (ngvec, vgc, ngvec, sfacg)

!_____________________________________
! generate the characteristic function

        Call gencfun

!___________________________________________
! generate structure factors for G+k-vectors

        Do ik = 1, nkpt
            Do ispn = 1, nspnfv
                Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), &
               &  ngkmax, sfacgk(:, :, ispn, ik))
            End Do
        End Do

!_________________________________________
! determine the new nuclear-nuclear energy

        Call energynn
      
        Return
    End Subroutine updatpos
    
end subroutine newton
