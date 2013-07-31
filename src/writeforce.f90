!
subroutine writeforce(fnum,verbosity)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
    integer, intent(in) :: verbosity
! local variables
    integer :: is, ia, ias, i
    
    if (verbosity < 1) return

    i=0
    write(fnum,*)
    if (input%groundstate%tfibs) then
        write(fnum,'(" Total atomic forces including IBS (lattice) :")')
    else
        write(fnum,'(" Total atomic forces (lattice) :")')
    end if 
    do is = 1, nspecies
        do ia = 1, natoms (is)
            ias = idxas (ia, is)
            i=i+1
            write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
           &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
           &  forcetot(:,ias)
        end do
    end do

   i=0
   if (verbosity > 1) then
        write(fnum,*)
        if (input%groundstate%tfibs) then
            write(fnum,'(" Atomic force components including IBS (lattice) :")')
        else
            write(fnum,'(" Atomic force components (lattice) :")')
        end if
        do is = 1, nspecies
            do ia = 1, natoms (is)
                ias = idxas (ia, is)
                i=i+1
                write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8,"   HF force")') &
               &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  forcehf(:,ias)
                write(fnum,'(                  T18,": ",3F14.8,"   core correction")') &
               &  forcecr(:,ias)
                if (input%groundstate%tfibs) &
               &write(fnum,'(                  T18,": ",3F14.8,"   IBS correction")') &
               &  forceibs(:,ias)
            end do
        end do
    end if
    
    return
end subroutine
