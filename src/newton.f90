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
        integer :: istep, is, ia
        
        istep = 0
        do while ((forcemax>forcetol).or.(istep>input%structureoptimization%maxsteps))
            
            istep = istep+1
            if (rank .Eq. 0) then
                write(60,*)
                write(60,'("+---------------------------------------------------------+")')
                write(60,'("| Optimization step ", I4,"    (method = newton)")') istep
                write(60,'("+---------------------------------------------------------+")')
! write lattice vectors and optimised atomic positions to file
                Call writehistory
                Call writegeometryxml(.True.)
! write the optimised interatomic distances to file
                Call writeiad(.True.)
            end if

! update atomic positions
            Call updatpos

! begin new self-consistent loop with updated positions
            call scf_cycle(0)

! output info
            if (rank==0) then
                write(60,'("Number of scf iterations : ", I4)') iscl
                write(60,'("Maximum force magnitude (target) : ",G18.10," (", G18.10, ")")') &
               &  forcemax, input%structureoptimization%epsforce
                write(*,'("Maximum force magnitude (target) : ",G18.10," (", G18.10, ")")') &
               &  forcemax, input%structureoptimization%epsforce
                write(60,'("Updated atomic positions :")')
                do is = 1, nspecies
                    do ia = 1, natoms (is)
                        write(60,'(" atom ",I4,4x,A2,T30,": ",3F14.8)') &
                       &  ia, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
                       &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                    end do
                end do
                if (input%structureoptimization%outputlevelnumber>0) call writeforce(60)
                call flushifc(60)
            end if
!
        end do
!
contains
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
! local variables
        Integer :: ik, ispn, is, ia, ias, i
        Real (8) :: t1
        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)
! compute the dot-product between the current and previous total force
                t1 = dot_product (forcetot(:, ias), forcetp(:, ias))
! if the force is in the same direction then increase step size parameter
                If (t1 .Gt. 0.d0) Then
                    tauatm(ias) = tauatm(ias)+input%structureoptimization%tau0atm
                Else
                    tauatm(ias) = input%structureoptimization%tau0atm
                End If
                do i = 1, 3
                    if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)) &
                   &   atposc(i,ia,is) = atposc(i,ia,is) + tauatm(ias)*(forcetot(i,ias)+forcetp(i,ias))
                end do ! i
            End Do
        End Do
        Do is = 1, nspecies
            Do ia = 1, natoms (is)
                ias = idxas (ia, is)
! compute the lattice coordinates of the atomic positions
                Call r3mv(ainv,atposc(:,ia,is),input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
! set the previous to the current total force
                forcetp(:, ias) = forcetot(:, ias)
            End Do
        End Do
! find the crystal symmetries and shift atomic positions if required
        Call findsymcrys
! check for overlapping muffin-tins
        Call checkmt
! generate structure factors for G-vectors
        Call gensfacgp (ngvec, vgc, ngvec, sfacg)
! generate the characteristic function
        Call gencfun
! generate structure factors for G+k-vectors
        Do ik = 1, nkpt
            Do ispn = 1, nspnfv
                Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), &
               &  ngkmax, sfacgk(:, :, ispn, ik))
            End Do
        End Do
! determine the new nuclear-nuclear energy
        Call energynn
      
        Return
    End Subroutine updatpos
    
end subroutine newton
