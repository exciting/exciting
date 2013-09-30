!
subroutine writepositions(fnum,verbosity)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
    integer, intent(in) :: verbosity
! local variables
    integer :: is, ia, i

    if (verbosity < 1) return

    i = 0
    write(fnum,*)
    write(fnum,'(" Atomic positions at this step (lattice) :")')
    do is = 1, nspecies
        do ia = 1, natoms (is)
          i = i+1
          write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
         &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
         &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
        end do
    end do
    
    return
end subroutine
