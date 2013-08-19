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
        character*(77) :: string

        nstep = 0
        atposcp(:,:,:) = atposc(:,:,:)
        atposcd(:,:,:) = atposc(:,:,:)

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
                Do is = 1, nspecies
                   Do ia = 1, natoms (is)
                       ias = idxas (ia, is)
                       forcetp(:, ias) = forcetot(:, ias)
                   End Do
                End Do
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
    subroutine updatpos_harmonic
        use modinput
        use modmain

!   Updates the current atomic positions according to the force on each atom. If
!   ${\bf r}_{ij}^m$ is the position and ${\bf F}_{ij}^m$ is the force acting on
!   it for atom $j$ of species $i$ and after time step $m$, then the new
!   position is calculated by an harmonic approximation


        Implicit None
        Integer :: ik, ispn, is, ia, ias, i
        Real (8) :: t1, beta, delta(3), zero(3), epsharmonic

        zero(:) = 0.d0
        epsharmonic = 1.0d-8
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        write(60,*) istep, "atposcp, forcetp"
!        Do is = 1, nspecies
!            Do ia = 1, natoms (is)
!                ias = idxas (ia, is)
!                write(60,'(6F12.8)') atposcp(:,ia,is), forcetp(:,ias)
!           End Do
!        End Do
!
!        write(60,*) istep, "atposc, forcetot"
!        Do is = 1, nspecies
!            Do ia = 1, natoms (is)
!                ias = idxas (ia, is)
!                write(60,'(6F12.8)') atposc(:,ia,is), forcetot(:,ias)
!           End Do
!        End Do
!
!        call flushifc(60)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!________________________
! update atomic positions


        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)
                tauatm(ias) = input%relax%tau0atm

!____________________________________________
! for the first step always use Newton method

                if (nstep<2) then
                    do i = 1, 3
                        atposcd(i,ia,is) = atposc(i,ia,is)
                        l1 = .not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)
                        l2 = ( abs(forcetot(i,ias)) .gt. input%structure%epslat/10.d0 )
                        if ( l1 .and. l2 ) then
                            atposc(i,ia,is) = atposc(i,ia,is) + tauatm(ias)*forcetot(i,ias)
                            atposcp(i,ia,is) = atposc(i,ia,is)
                        end if
                    end do
!____________________________________________
! otherwise use harmonic method

                else
                    do i = 1, 3
                        atposcd(i,ia,is) = atposc(i,ia,is)
                        l1 = .not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)


! forces with amplitude below epsharmonic are considered as zeros
!________________________________________________________________

                        l2 = ( abs(forcetot(i,ias)) .gt. epsharmonic )
                        if ( l1 .and. l2 ) then 
                            beta = forcetp(i,ias)/forcetot(i,ias)
                            atposc(i,ia,is) = atposcp(i,ia,is)/(1.-beta) - beta*atposc(i,ia,is)/(1.-beta)
                            atposcp(i,ia,is) = atposc(i,ia,is)
                        end if
                    end do
                 end if
             End Do
        End Do

!________________________________________________________
! compute the lattice coordinates of the atomic positions

        Do is = 1, nspecies
            Do ia = 1, natoms (is)

                Call r3mv(ainv,atposc(:,ia,is),input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))

            End Do
        End Do

!___________________________________________________________________
! find the crystal symmetries and shift atomic positions if required

        Call findsymcrys

        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)

!_____________________________________________________
! set the total forces for the previous configurations

                forcetp(:,ias) = forcetot(:,ias)

!__________________________________________________________
! if tshift .eq. true shift back the previous configuration

                if (input%structure%tshift) then
                    delta(:) = atposc(:,ia,is)-atposcp(:,ia,is)
                    atposcp(:,ia,is) = atposcd(:,ia,is)+delta(:)
                else 
                    atposcp(:,ia,is) = atposcd(:,ia,is)
                end if
            End Do
        End Do

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
    End Subroutine updatpos_harmonic
    
end subroutine harmonic
