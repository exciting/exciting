!
! Treatment of the symmetry requires averaging over degenerated states.
! Array n12dgn contains lower and upper indexes of the degenerated states.
! This array is used then in the routines to calculate the self-energies.
!
subroutine sym_state_degeneracy
    use modinput
    use modmain
    use modgw
    implicit none
    integer :: ik, ist, jst, i, n
    real(8) :: e0
    
    if (allocated(n12dgn)) deallocate(n12dgn)
    allocate(n12dgn(2,nstsv,nkpt))

    do ik = 1, nkpt

        if (input%gw%reduceq) then

            ist = 1
            do while (ist <= nstsv)
                e0 = evalsv(ist,ik)
! calculate the number of degenerated states
                n = 1; jst = ist+1
                do while (jst <= nstsv)
                    if (abs(evalsv(jst,ik)-e0) < 1.0d-4) then
                        n = n+1
                        jst = jst+1
                    else
                        exit
                    end if
                end do
! indexation of the degenerated states
                do i = 0, n-1
                    n12dgn(1,ist+i,ik) = ist
                    n12dgn(2,ist+i,ik) = ist+n-1
                end do
                ist = ist+n
            enddo ! ist
   
        else
   
! no symmetry: degeneracy index is set to 1
            do ist = 1, nstsv
                n12dgn(1,ist,ik) = ist
                n12dgn(2,ist,ik) = ist
            end do
   
        end if ! reduceq
   
! debug info
        if (input%gw%debug) then

            write(fdebug,*)
            write(fdebug,*)'Degeneracy state of KS bands: ik = ', ik
            do ist = 1, nstsv
                write(fdebug,*)'ist, n12dgn: ', ist, n12dgn(:,ist,ik)
            end do
            write(fdebug,*)
     
        end if ! debug

    enddo ! ik
    return
end subroutine

